//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file RefinementHandler.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_TURBINEREFINEMENTHANDLER_H
#define TURBINECORE_TURBINEREFINEMENTHANDLER_H

#pragma once

#include <memory>

#include <blockforest/StructuredBlockForest.h>
#include <field/EvaluationFilter.h>

#include "domain/refinement/RefinementType.h"
#include "domain/refinement/TurbineWeightAssignmentFunctor.h"
#include "domain/refinement/TurbineSelectorResetter.h"
#include "domain/refinement/LevelDeterminationFunctors.h"

namespace turbine_core {

    namespace refinement {

        namespace internal {

            template< typename PdfField_T, typename VectorField_T >
            struct VelocityCalculator {

                VelocityCalculator( walberla::StructuredBlockStorage * const storage,
                                    const BlockDataID & pdfFieldID, const BlockDataID & velocityFieldID )
                : storage_(storage), pdfFieldID_(pdfFieldID), velocityFieldID_(velocityFieldID)
                {}

                void operator()(std::vector<std::pair<const walberla::Block *, uint_t> > &,
                                std::vector<const walberla::Block *> &,
                                const walberla::BlockForest &) {

                    for (auto & block : *storage_) {

                        auto * pdfField = block.getData<PdfField_T>(pdfFieldID_);
                        WALBERLA_ASSERT_NOT_NULLPTR(pdfField)

                        auto * velocityField = block.getData<VectorField_T>(velocityFieldID_);
                        WALBERLA_ASSERT_NOT_NULLPTR(velocityField)

                        const auto ci = pdfField->xyzSizeWithGhostLayer();

                        for (const auto & cell : ci) {

                            walberla::Vector3<real_t> velocity{};
                            pdfField->getVelocity(velocity, cell);

                            if constexpr(VectorField_T::F_SIZE == 1) {
                                velocityField->get(cell) = velocity;
                            } else if constexpr(VectorField_T::F_SIZE == 3) {
                                velocityField->get(cell, 0) = velocity[0];
                                velocityField->get(cell, 1) = velocity[1];
                                velocityField->get(cell, 2) = velocity[2];
                            } else {
                                WALBERLA_ABORT("Invalid vector field.")
                            }
                        }
                    }
                }

                const BlockDataID pdfFieldID_{};
                const BlockDataID velocityFieldID_{};

                walberla::StructuredBlockStorage * const storage_{};

            };
        }

        template<typename WindFarm_T, typename Stencil_T, typename FlagField_T, typename VectorField_T, typename PdfField_T>
        class TurbineRefinementHandler {

            using Functor_T = std::function<void()>;
            using FlagFieldEvaluator_T = walberla::field::FlagFieldEvaluationFilter<FlagField_T>;

        public:

            TurbineRefinementHandler( const std::shared_ptr<walberla::StructuredBlockForest> & forest,
                                      WindFarm_T & farm, const uint_t nrGhostLayers,
                                      const walberla::Config::BlockHandle & domainConfig,
                                      StaticRefinementHandler<WindFarm_T> & staticRefinementHandler,
                                      const BlockDataID & pdfFieldID, const BlockDataID & velocityFieldID,
                                      const BlockDataID & flagFieldID, const walberla::FlagUID & fluidFlag )
            : structuredBlockForest_(forest.get()), forest_(&forest->getBlockForest()), farm_(&farm),
              nrGhostLayer_(nrGhostLayers), pdfFieldID_(pdfFieldID), velocityFieldID_(velocityFieldID),
              flagFieldID_(flagFieldID), fluidFlag_(fluidFlag), staticRefinementHandler_(&staticRefinementHandler)
            {
                walberla::Config::BlockHandle refinementConfig = domainConfig;

                if( domainConfig.getNumBlocks("Refinement") ) {
                    refinementConfig = domainConfig.getOneBlock("Refinement");
                }

                finestLevel_ = refinementConfig.getParameter<uint_t>("numLevels", uint_t{0});
                refreshInterval_ = refinementConfig.getParameter<uint_t>("refreshInterval", uint_t{0});

                forest_->recalculateBlockLevelsInRefresh( true );
                forest_->alwaysRebalanceInRefresh( false );
                forest_->reevaluateMinTargetLevelsAfterForcedRefinement( false );
                forest_->allowRefreshChangingDepth( false );

                forest_->checkForEarlyOutInRefresh(true);
                forest_->checkForLateOutInRefresh(true);

                bool curveHilbert = false;
                bool curveAllGather = true;
                forest_->setRefreshPhantomBlockMigrationPreparationFunction( walberla::blockforest::DynamicCurveBalance< walberla::blockforest::NoPhantomData >( curveHilbert, curveAllGather ) );

                setLevelDetermination( domainConfig );

            }

            void addRefreshCallbacks( const Functor_T & pdfCommunication,
                                      const Functor_T & boundarySetter,
                                      const Functor_T & forceFieldResetter ) {

                // if no refinement, don't call refresh callbacks
                if( finestLevel_ == 0 ) {
                    return;
                }

                using RefreshCallbackWrapper = walberla::blockforest::BlockForest::RefreshCallbackWrappper;

                // reset turbine selectors
                TurbineSelectorResetter<WindFarm_T> selectorResetter{farm_, nrGhostLayer_};
                forest_->addRefreshCallbackFunctionAfterBlockDataIsUnpacked(selectorResetter);
                //TODO setRefreshBlockStateDeterminationFunction?!

                forest_->addRefreshCallbackFunctionAfterBlockDataIsUnpacked(
                        RefreshCallbackWrapper(pdfCommunication)
                );

                forest_->addRefreshCallbackFunctionAfterBlockDataIsUnpacked(
                        RefreshCallbackWrapper(boundarySetter)
                );

                forest_->addRefreshCallbackFunctionAfterBlockDataIsUnpacked(
                        RefreshCallbackWrapper(forceFieldResetter)
                );

            }

            void addToTimeloop ( walberla::Timeloop * timeloop ) {
                if( refreshInterval_ != 0 ) {
                    WALBERLA_LOG_WARNING_ON_ROOT( "Dynamic refinement is currently not supported and all corresponding input parameter are ignored.\nThis will be fixed in the future." )
                    return;
                    timeloop->addFuncBeforeTimeStep(forest_->getRefreshFunctor(refreshInterval_),
                                                    "block forest refresh");
                }
            }

            auto getReducedRefreshTiming() {
                return forest_->getRefreshTiming().getReduced();
            }

            auto getRefinementType() const {
                return refinementType_;
            }

        private:

            walberla::StructuredBlockForest * structuredBlockForest_{};
            walberla::BlockForest * forest_{};
            WindFarm_T * farm_{};

            const uint_t nrGhostLayer_{};
            uint_t finestLevel_{};

            uint_t refreshInterval_{};

            const BlockDataID pdfFieldID_{};
            const BlockDataID velocityFieldID_{};

            const BlockDataID flagFieldID_{};
            const walberla::FlagUID fluidFlag_{};

            StaticRefinementHandler<WindFarm_T> * const staticRefinementHandler_{};
            RefinementType refinementType_{RefinementType::NONE};

            void setLevelDetermination( const walberla::Config::BlockHandle & domainConfig ) {

                // if no refinement, don't add level determination functions
                if( finestLevel_ == 0 ) {
                    return;
                }

                walberla::blockforest::MinTargetLevelDeterminationFunctions minTargetLevelDeterminationFunctions;

                /// STATIC REFINEMENT
                minTargetLevelDeterminationFunctions.add(*staticRefinementHandler_);

                /// weight assignment functors
                TurbineWeightAssignmentFunctor<WindFarm_T> weightAssignmentFunctor(farm_);
                forest_->setRefreshPhantomBlockDataAssignmentFunction(weightAssignmentFunctor);

                /// OPTIONAL REFINEMENTS BASED ON PARAMETER FILE

                // no refinement block -> early out
                if( domainConfig.getNumBlocks("Refinement") == 0 ) {
                    return;
                }

                const auto & refinementConfig = domainConfig.getOneBlock("Refinement");

                bool addedVelocityCalculation{false};
                internal::VelocityCalculator<PdfField_T,VectorField_T> velocityCalculator(structuredBlockForest_,pdfFieldID_,velocityFieldID_);

                FlagFieldEvaluator_T flagFieldFilter( flagFieldID_, fluidFlag_ );

                if( refinementConfig.getNumBlocks( "VelocityGradient" ) ) {

                    if( !addedVelocityCalculation ) {
                        minTargetLevelDeterminationFunctions.add( velocityCalculator );
                        addedVelocityCalculation = true;
                    }

                    const auto & velocityGradientConfig = refinementConfig.getOneBlock("VelocityGradient");
                    const real_t lowerLimit = velocityGradientConfig.getParameter<real_t>("lowerLimit");
                    const real_t upperLimit = velocityGradientConfig.getParameter<real_t>("upperLimit");

                    VectorGradientBasedLevelDetermination< Stencil_T, VectorField_T, FlagFieldEvaluator_T > gradientRefinement(
                            velocityFieldID_, flagFieldFilter, structuredBlockForest_, upperLimit, lowerLimit, finestLevel_
                    );

                    minTargetLevelDeterminationFunctions.add( gradientRefinement );

                    refinementType_ |= RefinementType::VELOCITY_GRADIENT;

                }

                if( refinementConfig.getNumBlocks( "Vorticity" ) ) {

                    if( !addedVelocityCalculation ) {
                        minTargetLevelDeterminationFunctions.add( velocityCalculator );
                        addedVelocityCalculation = true;
                    }

                    const auto & vorticityConfig = refinementConfig.getOneBlock("Vorticity");

                    const real_t lowerLimit = vorticityConfig.getParameter<real_t>("lowerLimit");
                    const real_t upperLimit = vorticityConfig.getParameter<real_t>("upperLimit");

                    VorticityBasedLevelDetermination< Stencil_T, VectorField_T, FlagFieldEvaluator_T > vorticityRefinement(
                            velocityFieldID_, flagFieldFilter, upperLimit, lowerLimit, finestLevel_
                    );

                    minTargetLevelDeterminationFunctions.add( vorticityRefinement );

                    refinementType_ |= RefinementType::VORTICITY;

                }

                if( refinementConfig.getNumBlocks( "Curl" ) ) {

                    if( !addedVelocityCalculation ) {
                        minTargetLevelDeterminationFunctions.add( velocityCalculator );
                        addedVelocityCalculation = true;
                    }

                    const auto & vorticityConfig = refinementConfig.getOneBlock("Curl");

                    const real_t lowerLimit = vorticityConfig.getParameter<real_t>("lowerLimit");
                    const real_t upperLimit = vorticityConfig.getParameter<real_t>("upperLimit");

                    CurlBasedLevelDetermination< FlagFieldEvaluator_T, VectorField_T > curlRefinement(
                            velocityFieldID_, structuredBlockForest_, flagFieldFilter, upperLimit, lowerLimit, finestLevel_
                    );

                    minTargetLevelDeterminationFunctions.add( curlRefinement );

                    refinementType_ |= RefinementType::CURL;

                }

                if( refinementConfig.getNumBlocks( "QCriterion" ) ) {

                    if( !addedVelocityCalculation ) {
                        minTargetLevelDeterminationFunctions.add( velocityCalculator );
                        addedVelocityCalculation = true;
                    }

                    const auto & qCriterionConfig = refinementConfig.getOneBlock("QCriterion");

                    const real_t lowerLimit = qCriterionConfig.getParameter<real_t>("lowerLimit");
                    const real_t upperLimit = qCriterionConfig.getParameter<real_t>("upperLimit");

                    QCriterionBasedLevelDetermination<Stencil_T,VectorField_T,FlagFieldEvaluator_T> qCriterionRefinement(
                            velocityFieldID_, flagFieldFilter, structuredBlockForest_,
                            upperLimit, lowerLimit, finestLevel_
                    );

                    minTargetLevelDeterminationFunctions.add( qCriterionRefinement );

                    refinementType_ |= RefinementType::Q_CRITERION;
                }

                forest_->setRefreshMinTargetLevelDeterminationFunction(minTargetLevelDeterminationFunctions);

            }

        };

    }

} // namespace turbine_core

#endif // TURBINECORE_TURBINEREFINEMENTHANDLER_H
