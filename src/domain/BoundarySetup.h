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
//! \file BoundarySetup.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_BOUNDARYSETUP_H
#define TURBINECORE_BOUNDARYSETUP_H

#pragma once

#include <core/StringUtility.h>
#include <core/config/Config.h>
#include <core/logging/Logging.h>
#include <field/GhostLayerField.h>
#include <geometry/InitBoundaryHandling.h>
// BEGIN EDIT: add waLBerla mesh headers so boundary setup can initialize STL/OBJ mesh obstacles.
#include <mesh/boundary/BoundarySetup.h>
#include <mesh_common/DistanceComputations.h>
#include <mesh_common/DistanceFunction.h>
#include <mesh_common/MeshIO.h>
#include <mesh_common/MeshOperations.h>
#include <mesh_common/TriangleMeshes.h>
#include <mesh_common/distance_octree/DistanceOctree.h>
// END EDIT: add waLBerla mesh headers so boundary setup can initialize STL/OBJ mesh obstacles.

#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/WalberlaDataTypes.h"
#include "EnvironmentSetup.h"
#include "wind_turbine_core/math/Vector3.h"
#include "BoundarySetter.h"

namespace turbine_core {

    namespace domain {

        template<typename BoundaryHandling_T, typename T>
        static constexpr bool is_flagField_v = std::is_base_of< BoundaryHandling_T, walberla::GhostLayerField<T,1> >::value;

        class BoundarySetup {

        public:
            HOST_PREFIX explicit BoundarySetup( const walberla::Config::BlockHandle & config);

            template< typename BoundaryHandling_T >
            HOST_PREFIX BoundarySetter<BoundaryHandling_T> fillFlagFieldFromSetter( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                            const BlockDataID & boundaryHandlingID,
                                            const uint_t nGhostLayer ) const;

            template< typename FlagField_T >
            HOST_PREFIX void fillFlagFieldFromConfig( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                      const BlockDataID & boundaryHandlingID,
                                                      const walberla::FlagUID & FluidFlagUID,
                                                      const walberla::FlagUID & NoSlipFlagUID,
                                                      const walberla::FlagUID & WFBFlagUID,
                                                      const walberla::FlagUID & SymmetryFlagUID,
                                                      const walberla::FlagUID & UniformInflowFlagUID,
                                                      const walberla::FlagUID & LogLawInflowFlagUID,
                                                      const walberla::FlagUID & OutflowFlagUID ) const;

            friend HOST_PREFIX std::ostream& operator<<(std::ostream& os, const BoundarySetup& setup);

            HOST_PREFIX inline walberla::Vector3<bool> periodicity() const { return periodicity_; }
            HOST_PREFIX inline walberla::Vector3<real_t> inflowVelocity() const { return inflowVelocity_; }

            HOST_PREFIX inline auto inflowType() const { return inflowType_; }
            HOST_PREFIX inline auto outflowType() const { return outflowType_; }
            HOST_PREFIX inline auto wallType() const { return wallType_; }
            HOST_PREFIX inline auto environmentSetup() const { return setup_; }

        private:

            // Tunnel or Open
            EnvironmentSetup::Type setup_;

            walberla::Vector3<bool> periodicity_;

            InflowSetup::Type inflowType_;
            OutflowSetup::Type outflowType_;
            WallSetup::Type wallType_;

            walberla::Vector3<real_t> inflowVelocity_;

            uint_t numberElements_;
            uint_t elementDistance_;
            uint_t elementLength_;

            walberla::Config::BlockHandle config_;

            HOST_PREFIX static void addBoundary( walberla::Config::Block & block, const std::string & dir, const std::string & flag ) {
                auto & border = block.createBlock("Border");
                border.addParameter("direction", dir);
                border.addParameter("walldistance", "-1");
                border.addParameter("flag", flag);
            }

            // BEGIN EDIT: add helpers that preserve existing primitive geometry input while enabling mesh-based obstacles.
            HOST_PREFIX static bool isInitializerGeometryBlock( const std::string & key ) {
                return key == "Border" || key == "CellInterval" || key == "Body" ||
                       key == "VoxelFile" || key == "Image" || key == "GrayScaleImage" ||
                       key == "RGBAImage";
            }

            HOST_PREFIX static void cloneConfigBlockRecursively( const walberla::Config::BlockHandle & source,
                                                                 walberla::Config::Block & target ) {
                for( auto it = source.begin(); it != source.end(); ++it ) {
                    target.setOrAddParameter( it->first, it->second );
                }

                walberla::Config::Blocks childBlocks;
                source.getBlocks( childBlocks );
                for( const auto & childBlock : childBlocks ) {
                    auto & targetChild = target.createBlock( childBlock.getKey() );
                    cloneConfigBlockRecursively( childBlock, targetChild );
                }
            }

            HOST_PREFIX static void copySupportedGeometryBlocks( const walberla::Config::BlockHandle & source,
                                                                 walberla::Config::Block & target ) {
                walberla::Config::Blocks childBlocks;
                source.getBlocks( childBlocks );
                for( const auto & childBlock : childBlocks ) {
                    if( !isInitializerGeometryBlock( childBlock.getKey() ) ) {
                        continue;
                    }

                    auto & targetChild = target.createBlock( childBlock.getKey() );
                    cloneConfigBlockRecursively( childBlock, targetChild );
                }
            }

            template< typename FlagField_T >
            HOST_PREFIX static uint_t getFlagFieldGhostLayers( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                               const BlockDataID & flagFieldID ) {
                for( auto & block : *sbf ) {
                    auto * flagField = block.template getData<FlagField_T>( flagFieldID );
                    WALBERLA_CHECK_NOT_NULLPTR( flagField, "flagFieldID invalid!" )
                    return flagField->nrOfGhostLayers();
                }

                WALBERLA_ABORT( "Cannot determine flag-field ghost layers because the block forest is empty." )
                return uint_t(0);
            }

            template< typename FlagField_T >
            HOST_PREFIX static void registerFlagOnAllBlocks( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                             const BlockDataID & flagFieldID,
                                                             const walberla::FlagUID & flagUID ) {
                for( auto & block : *sbf ) {
                    auto * flagField = block.template getData<FlagField_T>( flagFieldID );
                    WALBERLA_CHECK_NOT_NULLPTR( flagField, "flagFieldID invalid!" )
                    flagField->getOrRegisterFlag( flagUID );
                }
            }

            HOST_PREFIX static bool matchesFlagName( const std::string & requestedFlag,
                                                     const walberla::FlagUID & flagUID,
                                                     const std::string & alias ) {
                return walberla::string_icompare( requestedFlag, alias ) == 0 ||
                       walberla::string_icompare( requestedFlag, flagUID.getIdentifier() ) == 0;
            }

            HOST_PREFIX static const walberla::FlagUID & resolveMeshBoundaryFlag( const std::string & requestedFlag,
                                                                                  const walberla::FlagUID & NoSlipFlagUID,
                                                                                  const walberla::FlagUID & WFBFlagUID,
                                                                                  const walberla::FlagUID & SymmetryFlagUID,
                                                                                  const walberla::FlagUID & UniformInflowFlagUID,
                                                                                  const walberla::FlagUID & LogLawInflowFlagUID,
                                                                                  const walberla::FlagUID & OutflowFlagUID ) {
                if( matchesFlagName( requestedFlag, NoSlipFlagUID, "NoSlip" ) ) return NoSlipFlagUID;
                if( matchesFlagName( requestedFlag, WFBFlagUID, "WFB" ) ) return WFBFlagUID;
                if( matchesFlagName( requestedFlag, SymmetryFlagUID, "Symmetry" ) ) return SymmetryFlagUID;
                if( matchesFlagName( requestedFlag, UniformInflowFlagUID, "UniformInflow" ) ) return UniformInflowFlagUID;
                if( matchesFlagName( requestedFlag, LogLawInflowFlagUID, "LogLawInflow" ) ) return LogLawInflowFlagUID;
                if( matchesFlagName( requestedFlag, OutflowFlagUID, "Outflow" ) ) return OutflowFlagUID;

                WALBERLA_ABORT( "Unsupported mesh boundary flag '" << requestedFlag
                                << "'. Use one of: NoSlip, WFB, Symmetry, UniformInflow, LogLawInflow, Outflow, "
                                << "or the full generated flag identifiers." )
                return NoSlipFlagUID;
            }

            template< typename FlagField_T >
            HOST_PREFIX static void applyMeshObject( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                     const BlockDataID & flagFieldID,
                                                     const walberla::Config::BlockHandle & meshBlock,
                                                     const walberla::FlagUID & NoSlipFlagUID,
                                                     const walberla::FlagUID & WFBFlagUID,
                                                     const walberla::FlagUID & SymmetryFlagUID,
                                                     const walberla::FlagUID & UniformInflowFlagUID,
                                                     const walberla::FlagUID & LogLawInflowFlagUID,
                                                     const walberla::FlagUID & OutflowFlagUID ) {
                using Mesh_T = walberla::mesh::TriangleMesh;

                const std::string meshFile = meshBlock.isDefined( "file" )
                                             ? meshBlock.getParameter<std::string>( "file" )
                                             : meshBlock.getParameter<std::string>( "meshFile" );
                const bool binaryFile = meshBlock.getParameter<bool>( "binary", false );
                const std::string requestedFlag = meshBlock.getParameter<std::string>( "flag", "NoSlip" );
                // BEGIN EDIT: explicitly unwrap config parameters into concrete Vector3 values so NVCC does not route through the wrong Vector3 constructor.
                auto scaleParameter = meshBlock.getParameter<walberla::Vector3<real_t>>(
                        "scale", walberla::Vector3<real_t>( real_t(1) ) );
                const walberla::Vector3<real_t> scale = scaleParameter.operator walberla::Vector3<real_t>();

                auto offsetParameter = meshBlock.isDefined( "offset" )
                                       ? meshBlock.getParameter<walberla::Vector3<real_t>>( "offset" )
                                       : meshBlock.getParameter<walberla::Vector3<real_t>>(
                                               "translation", walberla::Vector3<real_t>( real_t(0) ) );
                const walberla::Vector3<real_t> offset = offsetParameter.operator walberla::Vector3<real_t>();
                // END EDIT: explicitly unwrap config parameters into concrete Vector3 values so NVCC does not route through the wrong Vector3 constructor.

                auto mesh = std::make_shared<Mesh_T>();
                walberla::mesh::readAndBroadcast( meshFile, *mesh, binaryFile );
                walberla::mesh::scale( *mesh, scale );
                walberla::mesh::translate( *mesh, offset );

                const auto meshBounds = walberla::mesh::computeAABB( *mesh );
                const auto & meshFlagUID = resolveMeshBoundaryFlag( requestedFlag, NoSlipFlagUID, WFBFlagUID,
                                                                    SymmetryFlagUID, UniformInflowFlagUID,
                                                                    LogLawInflowFlagUID, OutflowFlagUID );

                WALBERLA_LOG_INFO_ON_ROOT( "Applying mesh obstacle from '" << meshFile << "' with flag '"
                                           << meshFlagUID.getIdentifier() << "', scale " << scale
                                           << ", offset " << offset << ", bounds " << meshBounds )

                auto triDistance = std::make_shared<walberla::mesh::TriangleDistance<Mesh_T>>( mesh );
                auto distanceOctree = std::make_shared<walberla::mesh::DistanceOctree<Mesh_T>>( triDistance );

                registerFlagOnAllBlocks<FlagField_T>( sbf, flagFieldID, meshFlagUID );

                walberla::mesh::BoundarySetup meshBoundarySetup(
                        sbf,
                        walberla::makeMeshDistanceFunction( distanceOctree ),
                        getFlagFieldGhostLayers<FlagField_T>( sbf, flagFieldID ) );
                meshBoundarySetup.template setFlag<FlagField_T>(
                        flagFieldID, meshFlagUID, walberla::mesh::BoundarySetup::INSIDE );
            }
            // END EDIT: add helpers that preserve existing primitive geometry input while enabling mesh-based obstacles.

        };

        HOST_PREFIX BoundarySetup::BoundarySetup(const walberla::Config::BlockHandle &config) {

            setup_ = EnvironmentSetup::toType(config.getParameter<std::string>("setup"));
            config_ = config;

            switch (setup_) {
                case EnvironmentSetup::Periodic    : periodicity_ = walberla::Vector3<bool>(false , true , true ); break;
                case EnvironmentSetup::Open        : periodicity_ = walberla::Vector3<bool>(false, true , false); break;
                case EnvironmentSetup::Tunnel      : periodicity_ = walberla::Vector3<bool>(false, false, false); break;
                case EnvironmentSetup::FreeSlip    : periodicity_ = walberla::Vector3<bool>(false, true, false); break;
                default: WALBERLA_ABORT("Unknown boundary setup.")
            }

            inflowType_ = InflowSetup::toType(config.getParameter<std::string>("inflowType"));
            outflowType_ = OutflowSetup::toType(config.getParameter<std::string>("outflowType"));
            wallType_ = WallSetup::toType(config.getParameter<std::string>("wallType"));

            if((inflowType_ != InflowSetup::Periodic && outflowType_ == OutflowSetup::Periodic) ||
                (inflowType_ == InflowSetup::Periodic && outflowType_ != OutflowSetup::Periodic)) {
                WALBERLA_ABORT("Inflow and outflow type must be either both periodic or none of them is allowed to.");
            }

            if((inflowType_ != InflowSetup::ShiftedPeriodic && outflowType_ == OutflowSetup::ShiftedPeriodic) ||
               (inflowType_ == InflowSetup::ShiftedPeriodic && outflowType_ != OutflowSetup::ShiftedPeriodic)) {
                WALBERLA_ABORT("Inflow and outflow type must be either both shifted periodic or none of them is allowed to.");
            }

            if( inflowType_ == InflowSetup::Periodic ) {
                periodicity_[0] = true;
            } else {
                if( inflowType_ != InflowSetup::Periodic && inflowType_ != InflowSetup::ShiftedPeriodic )
                    inflowVelocity_ = config.getParameter<walberla::Vector3<real_t>>("inflowVelocity");
            }

            numberElements_ = config.getParameter<uint_t>("numberElements", 0);
            if(numberElements_) {
                elementDistance_ = config.getParameter<uint_t>("elementDistance", 0);
                elementDistance_ = uint_t(std::round(real_t(elementDistance_)));
                elementLength_ = config.getParameter<uint_t>("elementLength", 0);
                elementLength_ = uint_t(std::round(real_t(elementLength_)));
            }
        }

        template< typename BoundaryHandling_T >
        HOST_PREFIX BoundarySetter<BoundaryHandling_T> BoundarySetup::fillFlagFieldFromSetter( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                       const BlockDataID & boundaryHandlingID,
                                                       const uint_t nGhostLayer ) const {

            return BoundarySetter<BoundaryHandling_T> (sbf, setup_, inflowType_, outflowType_, boundaryHandlingID, nGhostLayer);
        }

        template< typename FlagField_T >
        HOST_PREFIX void BoundarySetup::fillFlagFieldFromConfig( const std::shared_ptr<walberla::StructuredBlockForest> & sbf,
                                                                 const BlockDataID & flagFieldID,
                                                                 const walberla::FlagUID & FluidFlagUID,
                                                                 const walberla::FlagUID & NoSlipFlagUID,
                                                                 const walberla::FlagUID & WFBFlagUID,
                                                                 const walberla::FlagUID & SymmetryFlagUID,
                                                                 const walberla::FlagUID & UniformInflowFlagUID,
                                                                 const walberla::FlagUID & LogLawInflowFlagUID,
                                                                 const walberla::FlagUID & OutflowFlagUID ) const {


            // BEGIN EDIT: filter primitive geometry through the existing initializer and then overlay optional mesh obstacles.
            walberla::Config::Block boundaryBlock( config_.getKey() );
            copySupportedGeometryBlocks( config_, boundaryBlock );

            if( inflowType_ == InflowSetup::InflowUniform ) {
                addBoundary(boundaryBlock, "W", UniformInflowFlagUID.getIdentifier());
            } else if (inflowType_ == InflowSetup::InflowLogLaw) {
                addBoundary(boundaryBlock, "W", LogLawInflowFlagUID.getIdentifier());
            }

            if( outflowType_ == OutflowSetup::Outflow ) {
                addBoundary(boundaryBlock, "E", OutflowFlagUID.getIdentifier());
            }

            if( setup_ == EnvironmentSetup::Open ) {
                if( wallType_ == WallSetup::NoSlip) {
                    addBoundary(boundaryBlock, "B", NoSlipFlagUID.getIdentifier());
                } else if( wallType_ == WallSetup::WFB ) {
                    addBoundary(boundaryBlock, "B", WFBFlagUID.getIdentifier());
                } else {
                    WALBERLA_ABORT("Invalid wall type in Boundary Configuration")
                }
                addBoundary(boundaryBlock, "T", SymmetryFlagUID.getIdentifier());
            } else if ( setup_ == EnvironmentSetup::Tunnel ) {
                if( wallType_ == WallSetup::NoSlip) {
                    addBoundary(boundaryBlock, "B,T", NoSlipFlagUID.getIdentifier());
                } else if( wallType_ == WallSetup::WFB ) {
                    addBoundary(boundaryBlock, "B,T", WFBFlagUID.getIdentifier());
                } else {
                    WALBERLA_ABORT("Invalid wall type in Boundary Configuration")
                }
            } else if ( setup_ == EnvironmentSetup::FreeSlip ) {
                addBoundary(boundaryBlock, "B", SymmetryFlagUID.getIdentifier());
                addBoundary(boundaryBlock, "T", SymmetryFlagUID.getIdentifier());
            }

            if(numberElements_) {

                const auto domainAABB = sbf->getDomain();
                const auto yMax = int(domainAABB.yMax()) - 3;

                // create all elements lines
                for (uint_t i = 1; i <= numberElements_; ++i) {

                    auto xPos = i * elementDistance_;
                    if (xPos > uint_t(domainAABB.xMax()) - 3)
                        break;

                    // create element segments
                    int yPos = -1 + int(i % 2) * int(real_t(elementLength_) * real_t(0.3));
                    while (yPos < yMax + 1) {
                        auto &ciBlock = boundaryBlock.createBlock("CellInterval");

                        std::string minCI(
                                "<" + std::to_string(xPos) + ", " + std::to_string(std::max(yPos, 3)) + ", 0>");
                        std::string maxCI("<" + std::to_string(xPos + 1) + ", " +
                                          std::to_string(std::min(yPos + int(elementLength_), yMax)) + ", " +
                                          std::to_string(elementLength_ ) + ">");

                        ciBlock.addParameter("min", minCI);
                        ciBlock.addParameter("max", maxCI);
                        ciBlock.addParameter("flag", NoSlipFlagUID.getIdentifier());

                        yPos += int(real_t(1.5) * real_t(elementLength_));
                    }
                }
            }

            walberla::Config::BlockHandle boundariesConfig(&boundaryBlock);
            walberla::geometry::initBoundaryHandling< FlagField_T >( *sbf, flagFieldID, boundariesConfig );

            walberla::Config::Blocks meshBlocks;
            config_.getBlocks( "MeshObject", meshBlocks );

            walberla::Config::Blocks stlBlocks;
            config_.getBlocks( "STLObject", stlBlocks );
            meshBlocks.insert( meshBlocks.end(), stlBlocks.begin(), stlBlocks.end() );

            for( const auto & meshBlock : meshBlocks ) {
                applyMeshObject<FlagField_T>( sbf, flagFieldID, meshBlock, NoSlipFlagUID, WFBFlagUID,
                                              SymmetryFlagUID, UniformInflowFlagUID,
                                              LogLawInflowFlagUID, OutflowFlagUID );
            }

            walberla::geometry::setNonBoundaryCellsToDomain<FlagField_T>(*sbf, flagFieldID, FluidFlagUID);
            // END EDIT: filter primitive geometry through the existing initializer and then overlay optional mesh obstacles.
        }

        HOST_PREFIX std::ostream &operator<<(std::ostream &os, const BoundarySetup &setup) {
            os << "BOUNDARY SETUP INFORMATION:\n";
            os << "\n    - Environmental Setup:         " << EnvironmentSetup::toString(setup.setup_);
            os << "\n    - Periodicity:                 " << setup.periodicity_;
            os << "\n    - Boundary Types:              "
               << "\n        Inflow:                    " << InflowSetup::toString(setup.inflowType_);
            if( setup.inflowType() != InflowSetup::Periodic && setup.inflowType() != InflowSetup::ShiftedPeriodic )
                os << "\n        Inflow Velocity:           " << InflowSetup::toString(setup.inflowType_);

            return os;
        }

    } // domain
} // turbine_core

#endif
