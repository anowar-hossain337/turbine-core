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
//! \file CustomRecursiveTimeStep.h
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#pragma once

#include <blockforest/communication/NonUniformBufferedScheme.h>

#include <lbm/field/PdfField.h>
#include <lbm_generated/communication/NonuniformGeneratedPdfPackInfo.h>

#include <timeloop/SweepTimeloop.h>

#include "wind_turbine_core/WalberlaDataTypes.h"
#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    using walberla::blockforest::communication::NonUniformBufferedScheme;

    namespace refinement {

        /**
         *
         * @tparam LatticeStorageSpecification_T   Generated storage specification
         * @tparam SweepCollection_T LBM SweepCollection (must be able to call stream, collide, streamCollide and streamOnlyNoAdvancement)
         * @tparam BoundaryCollection_T LBM Boundary collection (Functor that runs all boundary kernels at call)
         */
        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        class CustomRecursiveTimeStep {
        public:
            using LatticeStorageSpecification_T = typename PdfField_T::LatticeStorageSpecification;
            using Stencil = typename LatticeStorageSpecification_T::Stencil;
            using CommunicationStencil = typename LatticeStorageSpecification_T::CommunicationStencil;
            using CommScheme = NonUniformBufferedScheme<CommunicationStencil>;
            using PackInfo = walberla::lbm_generated::NonuniformGeneratedPdfPackInfo<PdfField_T>;

            using VoidFunction = std::function<void()>;
            using FunctionTuple = std::tuple<std::string, VoidFunction, uint_t>;

            CustomRecursiveTimeStep(std::shared_ptr<walberla::StructuredBlockForest> &sbfs,
                                    const BlockDataID &pdfFieldId, SweepCollection_T &sweepCollection,
                                    BoundaryCollection_T &boundaryCollection,
                                    std::shared_ptr<CommScheme> &commScheme,
                                    std::shared_ptr<PackInfo> &pdfFieldPackInfo) :
                    sbfs_(sbfs), pdfFieldId_(pdfFieldId), pdfFieldPackInfo_(pdfFieldPackInfo), commScheme_(commScheme),
                    sweepCollection_(sweepCollection), boundaryCollection_(boundaryCollection) {
#ifndef NDEBUG
                for (auto& block : *sbfs)
                   WALBERLA_ASSERT(block.isDataOfType<PdfField_T>(pdfFieldId_), "Template parameter PdfField_T is of different type than BlockDataID pdfFieldId that is provided as constructor argument")
#endif
                maxLevel_ = sbfs->getDepth();

                for (uint_t level = 0; level <= maxLevel_; level++) {
                    std::vector<walberla::Block *> blocks;
                    sbfs->getBlocks(blocks, level);
                    blocks_.push_back(blocks);
                }
            };

            void operator()() { timestep(0); };

            void addRefinementToTimeLoop(walberla::SweepTimeloop &timeloop, uint_t level = 0);

            void addAfterStreamFunction(const VoidFunction &fct, std::string name, const uint_t level) {
                afterStreamFunction_.emplace_back(name, fct, level);
            }
            void addAfterStreamFunction(const VoidFunction &fct, std::string name = "") {
                addAfterStreamFunction(fct, name, maxLevel_);
            }

            void addAfterCollideFunction(const VoidFunction &fct, std::string name, const uint_t level) {
                afterStreamFunction_.emplace_back(name, fct, level);
            }
            void addAfterCollideFunction(const VoidFunction &fct, std::string name = "") {
                addAfterCollideFunction(fct, name, maxLevel_);
            }

            void addAfterBoundaryFunction(const VoidFunction &fct, std::string name, const uint_t level) {
                afterStreamFunction_.emplace_back(name, fct, level);
            }
            void addAfterBoundaryFunction(const VoidFunction &fct, std::string name = "") {
                addAfterBoundaryFunction(fct, name, maxLevel_);
            }

        private:
            void timestep(uint_t level);

            void ghostLayerPropagation(walberla::Block *block);

            std::function<void()> executeStreamOnLevel(uint_t level, bool withGhostLayerPropagation = false);

            std::function<void()> executeCollideOnLevel(uint_t level, bool withGhostLayerPropagation = false);

            std::function<void()> executeBoundaryHandlingOnLevel(uint_t level);

            std::shared_ptr<walberla::StructuredBlockForest> sbfs_;
            uint_t maxLevel_;
            std::vector<std::vector<walberla::Block *>> blocks_;

            std::vector<FunctionTuple> afterStreamFunction_;
            std::vector<FunctionTuple> afterCollideFunction_;
            std::vector<FunctionTuple> afterBoundaryFunction_;

            const BlockDataID pdfFieldId_;
            std::shared_ptr<PackInfo> pdfFieldPackInfo_;
            std::shared_ptr<CommScheme> commScheme_;

            SweepCollection_T &sweepCollection_;
            BoundaryCollection_T &boundaryCollection_;
        };

        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        void CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::timestep(uint_t level) {
            // 1.1 Collision
            for (auto b: blocks_[level]) {
                sweepCollection_.stream(b);
                for (auto &fct: afterStreamFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
                sweepCollection_.collide(b);
                for (auto &fct: afterCollideFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
            }

            // 1.2 Recursive Descent
            if (level < maxLevel_) {
                timestep(level + 1);
            }

            // 1.3 Coarse to Fine Communication, receiving end
            if (level != 0) {
                commScheme_->communicateCoarseToFine(level);
            }

            // 1.4 Equal-Level Communication
            commScheme_->communicateEqualLevel(level);

            // 1.5 Boundary Handling and Coalescence Preparation
            for (auto b: blocks_[level]) {
                boundaryCollection_(b);
                for (auto &fct: afterBoundaryFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
                if (level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
            }

            // 1.6 Fine to Coarse Communication, receiving end
            if (level < maxLevel_) {
                commScheme_->communicateFineToCoarse(level + 1);
            }

            // Stop here if on coarsest level.
            // Otherwise, continue to second subcycle.
            if (level == 0) return;

            // 2.1 Collision and Ghost-Layer Propagation
            for (auto b: blocks_[level]) {
                ghostLayerPropagation(b);  // GL-Propagation first without swapping arrays...
                sweepCollection_.stream(b);                // then Stream-Collide on interior, and swap arrays
                for (auto &fct: afterStreamFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
                sweepCollection_.collide(b);
                for (auto &fct: afterCollideFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
            }

            // 2.2 Recursive Descent
            if (level < maxLevel_) {
                timestep(level + 1);
            }

            // 2.4 Equal-Level Communication
            commScheme_->communicateEqualLevel(level);

            // 2.5 Boundary Handling and Coalescence Preparation
            for (auto b: blocks_[level]) {
                boundaryCollection_(b);
                for (auto &fct: afterBoundaryFunction_) {
                    if(level < std::get<2>(fct) + 1)
                        std::get<1>(fct)();
                }
                if (level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
            }

            // 2.6 Fine to Coarse Communication, receiving end
            if (level < maxLevel_) {
                commScheme_->communicateFineToCoarse(level + 1);
            }
        }


        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        void CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::addRefinementToTimeLoop(
                walberla::SweepTimeloop &timeloop, uint_t level) {
            // 1.1 Collision
            timeloop.addFuncBeforeTimeStep(executeStreamOnLevel(level),
                                           "Refinement Cycle: stream on level " + std::to_string(level));
            for (auto &fct: afterStreamFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }

            timeloop.addFuncBeforeTimeStep(executeCollideOnLevel(level),
                                           "Refinement Cycle: collide on level " + std::to_string(level));
            for (auto &fct: afterCollideFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }

            // 1.2 Recursive Descent
            if (level < maxLevel_) {
                addRefinementToTimeLoop(timeloop, level + 1);
            }

            // 1.3 Coarse to Fine Communication, receiving end
            if (level != 0) {
                timeloop.addFuncBeforeTimeStep(commScheme_->communicateCoarseToFineFunctor(level),
                                               "Refinement Cycle: communicate coarse to fine on level " +
                                               std::to_string(level));
            }

            // 1.4 Equal-Level Communication
            timeloop.addFuncBeforeTimeStep(commScheme_->communicateEqualLevelFunctor(level),
                                           "Refinement Cycle: communicate equal level on level " +
                                           std::to_string(level));


            // 1.5 Boundary Handling and Coalescence Preparation
            timeloop.addFuncBeforeTimeStep(executeBoundaryHandlingOnLevel(level),
                                           "Refinement Cycle: boundary handling on level " + std::to_string(level));
            for (auto &fct: afterBoundaryFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }

            // 1.6 Fine to Coarse Communication, receiving end
            if (level < maxLevel_) {
                timeloop.addFuncBeforeTimeStep(commScheme_->communicateFineToCoarseFunctor(level + 1),
                                               "Refinement Cycle: communicate fine to coarse on level " +
                                               std::to_string(level + 1));
            }

            // Stop here if on coarsest level.
            // Otherwise, continue to second subcycle.
            if (level == 0) return;

            // 2.1 Collision and Ghost-Layer Propagation
            timeloop.addFuncBeforeTimeStep(executeStreamOnLevel(level, true),
                                           "Refinement Cycle: stream with ghost layer propagation on level " +
                                           std::to_string(level));
            for (auto &fct: afterStreamFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }
            timeloop.addFuncBeforeTimeStep(executeCollideOnLevel(level, true),
                                           "Refinement Cycle: collide with ghost layer propagation on level " +
                                           std::to_string(level));
            for (auto &fct: afterCollideFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }

            // 2.2 Recursive Descent
            if (level < maxLevel_)
                addRefinementToTimeLoop(timeloop, level + 1);


            // 2.4 Equal-Level Communication
            timeloop.addFuncBeforeTimeStep(commScheme_->communicateEqualLevelFunctor(level),
                                           "Refinement Cycle: communicate equal level on level " +
                                           std::to_string(level));

            // 2.5 Boundary Handling and Coalescence Preparation
            timeloop.addFuncBeforeTimeStep(executeBoundaryHandlingOnLevel(level),
                                           "Refinement Cycle: boundary handling on level " + std::to_string(level));
            for (auto &fct: afterBoundaryFunction_) {
                if(level < std::get<2>(fct) + 1)
                    timeloop.addFuncBeforeTimeStep(std::get<1>(fct), std::get<0>(fct));
            }

            // 2.6 Fine to Coarse Communication, receiving end
            if (level < maxLevel_)
                timeloop.addFuncBeforeTimeStep(commScheme_->communicateFineToCoarseFunctor(level + 1),
                                               "Refinement Cycle: communicate fine to coarse on level " +
                                               std::to_string(level + 1));

        }


        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        std::function<void()>
        CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::executeStreamOnLevel(uint_t level,
                                                                                                           bool withGhostLayerPropagation) {
            return [level, withGhostLayerPropagation, this]() {
                if (withGhostLayerPropagation) {
                    for (auto b: blocks_[level]) {
                        ghostLayerPropagation(b);
                        sweepCollection_.stream(b);
                    }
                } else {
                    for (auto b: blocks_[level]) {
                        sweepCollection_.stream(b);
                    }
                }
            };
        }

        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        std::function<void()>
        CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::executeCollideOnLevel(
                uint_t level, bool withGhostLayerPropagation) {
            return [level, withGhostLayerPropagation, this]() {
                if (withGhostLayerPropagation) {
                    for (auto b: blocks_[level]) {
                        ghostLayerPropagation(b);
                        sweepCollection_.collide(b);
                    }
                } else {
                    for (auto b: blocks_[level]) {
                        sweepCollection_.collide(b);
                    }
                }
            };
        }


        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        std::function<void()>
        CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::executeBoundaryHandlingOnLevel(
                uint_t level) {
            return [level, this]() {
                for (auto b: blocks_[level]) {
                    boundaryCollection_(b);
                    if (level != maxLevel_) pdfFieldPackInfo_->prepareCoalescence(b);
                }
            };
        }


        template<typename PdfField_T, typename SweepCollection_T, typename BoundaryCollection_T>
        void CustomRecursiveTimeStep<PdfField_T, SweepCollection_T, BoundaryCollection_T>::ghostLayerPropagation(
                walberla::Block *block) {
            auto pdfField = block->getData<PdfField_T>(pdfFieldId_);

            for (auto it = CommunicationStencil::beginNoCenter(); it != CommunicationStencil::end(); ++it) {
                uint_t nSecIdx = walberla::blockforest::getBlockNeighborhoodSectionIndex(*it);
                // Propagate on ghost layers shadowing coarse or no blocks
                if (block->neighborhoodSectionHasLargerBlock(nSecIdx)) {
                    walberla::CellInterval ci;
                    pdfField->getGhostRegion(*it, ci, 1);
                    sweepCollection_.streamOnlyNoAdvancementCellInterval(block, ci);
                }
            }
        }

    } // namespace refinement
} // namespace turbine_core
