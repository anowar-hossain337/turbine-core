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
//! \file DomainSetup.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_DOMAINSETUP_H
#define TURBINECORE_DOMAINSETUP_H

#pragma once

#include <blockforest/SetupBlockForest.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/config/Config.h>

#include <blockforest/AABBRefinementSelection.h>
#include <blockforest/loadbalancing/all.h>
#include <blockforest/Initialization.h>
#include <core/math/Limits.h>

#include <vtk/VTKOutput.h>

#include "SnapshotHandler.h"
#include "refinement/StaticRefinementHandler.h"

#include "wind_turbine_core/math/Round.h"
#include "wind_turbine_core/math/AABB.h"
#include "wind_turbine_core/Selectors.h"

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace domain {

        static void workloadAndMemoryAssignment( walberla::SetupBlockForest &forest, const walberla::memory_t memoryPerBlock, const std::vector<AABB> & turbineAABBs ) {
            for (auto &block : forest) {

                uint_t level = block.getLevel();
                auto workLoad = real_t(uint_t(1) << level);

                 for (uint_t i = 0; i < turbineAABBs.size(); ++i) {

                     walberla::AABB aabb{walberla::Vector3<real_t>(turbineAABBs[i].min().data()),
                                         walberla::Vector3<real_t>(turbineAABBs[i].max().data())};

                     if (block.getAABB().intersects(aabb)) {
                         block.addState(sets::turbineSelector());
                         block.addState(sets::turbineSelector(i));

                         workLoad *= real_t(1.5);
                         break;
                     }
                 }

                 if( level != uint_t(0) ) {
                     block.addState(sets::refinementSelector());
                 }

                 if( block.getState().empty() ) {
                     block.addState(sets::empty());
                 }

                block.setWorkload(walberla::numeric_cast<walberla::workload_t>(workLoad));
                block.setMemory(memoryPerBlock);
            }
        }

        struct DomainSetup {

            explicit DomainSetup(std::shared_ptr<walberla::Config> & globalConfig,
                                 const walberla::Vector3<bool> & periodicity, const uint_t q = 27);

            friend std::ostream& operator<<(std::ostream& os, const DomainSetup& setup);

            template< typename StaticRefinementHandler_T >
            std::shared_ptr< walberla::blockforest::StructuredBlockForest > createStructuredBlockForest( StaticRefinementHandler_T & turbineAabbRefinementSelector, const uint_t fieldGhostLayers );
            std::shared_ptr< walberla::blockforest::StructuredBlockForest > createUniformBlockForest( const std::vector<AABB> & turbineAABBs );

            void computeDomainDecomposition();
            bool adjustDomainSizeToBlockSize();

            // number of refinement levels
            uint_t numLevels_{0};

            // scaling factor from fine dx to coarse dx, e.g., refinement 1 -> dx_fine = 1.0; dx_coarse = dx_fine * scaling = 1.0 * 2.0
            const uint_t dx_fine_;
            uint_t dx_coarse_;

            // domain size
            walberla::Vector3<bool> periodicity_;
            walberla::Vector3<uint_t> domainSize_;
            walberla::Vector3<uint_t> cellsPerBlock_;
            walberla::Vector3<uint_t> numberOfCoarseBlocksPerDirection_;

            real_t domainAdjustmentThreshold_;
            bool fixedCellsPerBlock_;

            // memory defines
            walberla::memory_t memoryPerCell_;
            walberla::memory_t processMemoryLimit_;

            // IO defines
            bool outputSetupForest_;

            // snapshot
            std::unique_ptr<SnapshotHandler> snapshotHandler_;

        };

        DomainSetup::DomainSetup(std::shared_ptr<walberla::Config> & globalConfig,
                                 const walberla::Vector3<bool> & periodicity, const uint_t q)
                : dx_fine_(1.0), periodicity_(periodicity),
                  numberOfCoarseBlocksPerDirection_(4) {

            const auto config = globalConfig->getOneBlock("DomainSetup");

            domainSize_ = config.getParameter<walberla::Vector3<uint_t>>("domainSize");

            if( config.getNumBlocks("Refinement") ) {
                numLevels_ = config.getOneBlock("Refinement").getParameter<uint_t>("numLevels");
            }

            dx_coarse_ = dx_fine_ << numLevels_;

            // ... in bytes
            memoryPerCell_ = config.getParameter<walberla::memory_t>("memoryPerCell", walberla::memory_t(q * 8 + 1));
            processMemoryLimit_ = config.getParameter<walberla::memory_t>("processMemoryLimit", walberla::math::Limits<walberla::memory_t>::inf());

            if (config.isDefined("cellsPerBlock")) {
                if (config.isDefined("maxCellsPerBlock")) {
                    WALBERLA_ABORT("Please specify either fixed number of cellsPerBlock or maximum cellsPerBlock, not both!")
                }

                cellsPerBlock_ = config.getParameter<walberla::Vector3<uint_t>>("cellsPerBlock");
                fixedCellsPerBlock_ = true;
            } else {
                cellsPerBlock_ = config.getParameter<walberla::Vector3<uint_t>>("maxCellsPerBlock", walberla::Vector3<uint_t>(48, 32, 32));
                fixedCellsPerBlock_ = false;
            }

            domainAdjustmentThreshold_ = config.getParameter<real_t>("domainAdjustmentThreshold", real_t(0.5));

            outputSetupForest_ = config.getParameter<bool>("outputSetupForest", false);

            snapshotHandler_ = std::make_unique<SnapshotHandler>(globalConfig);

        }

        /*!
        * @brief Stream domain setup information to std::ostream.
        * @param os
        * @param setup
        * @return
        */
        std::ostream& operator<<(std::ostream& os, const DomainSetup& setup) {
            os << "DOMAIN SETUP INFORMATION:\n";
            os << "\n    - domain size:                 " << setup.domainSize_
               << "\n    - cells per block:             " << setup.cellsPerBlock_
               << "\n    - number of refinement levels: " << setup.numLevels_
               << "\n        -> dx_coarse:              " << setup.dx_coarse_;

            if(setup.numLevels_)
                os << "\n        -> dx_fine:                " << setup.dx_fine_;

            return os;
        }

        /*!
        * @brief Based on a given block size (number of cells per block), this function increases or decreases the domain size.
        * This is needed to have complete blocks in every spatial direction.
        * @note The decision on whether to increase or decrease the domain size is based on the parameter `domainAdjustmentThreshold_`.
        * For very small values, the function will most probably increase the domain size, for values that tend to 1, the domain size will be decreased.
        */
        bool DomainSetup::adjustDomainSizeToBlockSize() {

            bool changedSize{false};

            for(uint_t d = 0; d < 3; ++d) {

                const uint_t factor = cellsPerBlock_[d] * dx_coarse_;
                numberOfCoarseBlocksPerDirection_[d] = std::max(uint_t( math::round(real_t(domainSize_[d]) / real_t(factor), domainAdjustmentThreshold_)), uint_t(1) );

                const uint_t nFineCells = numberOfCoarseBlocksPerDirection_[d] * factor;

                if(nFineCells != uint_t(domainSize_[d])) {
                    domainSize_[d] = nFineCells;
                    changedSize = true;
                }

            }

            return changedSize;

        }

        /*!
        * @brief Computes the domain decomposition for either a fixed block size or a given maximum block size.
        */
        void DomainSetup::computeDomainDecomposition() {

            bool changedSize;

            if(fixedCellsPerBlock_) {
                changedSize = adjustDomainSizeToBlockSize();
            } else {
                //TODO do something clever to avoid too elongated blocks -> some threshold for side ratio?
                const walberla::Vector3<uint_t> maxCellsPerBlock = cellsPerBlock_;

                // obey maximum number of cells per block
                for(uint_t d = 0; d < 3; ++d) {

                    if( domainSize_[d] < cellsPerBlock_[d] ) {
                        WALBERLA_ABORT("Your domain is smaller than the block size. Please use a fixed block size in this case.")
                    } else {
                        const real_t nCoarseCells = real_t(domainSize_[d]) / real_t(dx_coarse_);

                        cellsPerBlock_[d] = uint_t(math::roundToMultiple<2>(nCoarseCells / real_t(numberOfCoarseBlocksPerDirection_[d])));
                        while (cellsPerBlock_[d] > maxCellsPerBlock[d]) {
                            numberOfCoarseBlocksPerDirection_[d] += uint_t(1);
                            cellsPerBlock_[d] = uint_t(math::roundToMultiple<2>(nCoarseCells / real_t(numberOfCoarseBlocksPerDirection_[d])));
                        }
                    }

                    if( cellsPerBlock_[d] < 16 ) {
                        WALBERLA_LOG_WARNING_ON_ROOT("You're using only " << cellsPerBlock_[d] << " cells per block in direction " << d << ". Are you sure about this?!")
                    }
                }

                changedSize = adjustDomainSizeToBlockSize();
            }

            if( changedSize ) {
                WALBERLA_LOG_WARNING_ON_ROOT("Domain size was changed to fit the number of cells per Block. New domain size : " << domainSize_)
            }

        } // function computeDomainDecomposition

        ////////////////////////////
        /// BLOCKFOREST CREATION ///
        ////////////////////////////

        template< typename StaticRefinementHandler_T >
        std::shared_ptr< walberla::blockforest::StructuredBlockForest >
        DomainSetup::createStructuredBlockForest( StaticRefinementHandler_T & staticRefinementHandler, const uint_t fieldGhostLayers ) {

            std::string blockforestFileName;
            if(snapshotHandler_->loadSnapshot || snapshotHandler_->storeSnapshot) {
                blockforestFileName = std::string(walberla::filesystem::path(snapshotHandler_->baseFolder) / walberla::filesystem::path(snapshotHandler_->blockforestFile).c_str());
            }

            if(snapshotHandler_->loadSnapshot) {

                WALBERLA_LOG_INFO_ON_ROOT("Reading the block structure from file " << snapshotHandler_->blockforestFile << "...")

                // load block forest from file
                walberla::MPIManager::instance()->useWorldComm();

                auto bf = std::make_shared< walberla::BlockForest >(
                        uint_t(walberla::MPIManager::instance()->rank()),
                        blockforestFileName.c_str(),
                        true, false);

                // save data from blockforest for domain decomposition calculation
                const auto & aabb = bf->getDomain();
                numLevels_ = bf->getNumberOfLevels();
                for(uint_t d = 0; d < 3; ++d) {
                    periodicity_[d] = bf->isPeriodic(d);
                    domainSize_[d] = uint_t(aabb.max(d));
                    numberOfCoarseBlocksPerDirection_[d] = bf->getSize(d);
                    cellsPerBlock_[d] = domainSize_[d] / numberOfCoarseBlocksPerDirection_[d];
                }

                fixedCellsPerBlock_ = true;

            }

            WALBERLA_LOG_INFO_ON_ROOT("Creating the block structure ...")

            /// SetupBlockForest
            std::shared_ptr< walberla::SetupBlockForest > forest = std::make_shared<walberla::SetupBlockForest>();

            computeDomainDecomposition();

            /// workload
            const walberla::memory_t memoryPerBlock = walberla::numeric_cast<walberla::memory_t>(
                    (cellsPerBlock_[0] + uint_t(2) * fieldGhostLayers)  *
                    (cellsPerBlock_[1] + uint_t(2) * fieldGhostLayers)  *
                    (cellsPerBlock_[2] + uint_t(2) * fieldGhostLayers)) * memoryPerCell_;

            const auto turbineAABBs = staticRefinementHandler.getTurbineAABBs();

            forest->addWorkloadMemorySUIDAssignmentFunction(
                    std::bind(workloadAndMemoryAssignment, std::placeholders::_1, memoryPerBlock, turbineAABBs));

            /// refinement selection
            walberla::blockforest::RefinementSelectionFunctions refinementSelectionFunctions;
            refinementSelectionFunctions.add(staticRefinementHandler);
            forest->addRefinementSelectionFunction(refinementSelectionFunctions);

            auto domainAABB = walberla::AABB(real_t(0), real_t(0), real_t(0),
                                             real_t(domainSize_[0]), real_t(domainSize_[1]), real_t(domainSize_[2]));

            forest->init(domainAABB,
                         numberOfCoarseBlocksPerDirection_[0], numberOfCoarseBlocksPerDirection_[1],
                         numberOfCoarseBlocksPerDirection_[2],
                         periodicity_[0], periodicity_[1], periodicity_[2]);

            if( !walberla::MPIManager::instance()->rankValid() )
                walberla::MPIManager::instance()->useWorldComm();

            auto numberOfProcesses = uint_t(walberla::MPIManager::instance()->numProcesses());

            WALBERLA_CHECK(forest->getNumberOfBlocks(forest->getDepth()) >= numberOfProcesses,
                           "You have to have at least as many fine blocks as MPI processes in mesh refinement.")


            forest->balanceLoad(walberla::blockforest::StaticLevelwiseCurveBalanceWeighted(true), numberOfProcesses, real_t(0),
                                processMemoryLimit_, false);

            WALBERLA_ROOT_SECTION() {
                if (outputSetupForest_) {
                    forest->writeVTKOutput("initialDomainDecomposition");
                    forest->writeCSV("process_distribution");
                }
            }

            WALBERLA_CHECK_EQUAL(forest->getNumberOfBufferProcesses(), 0, "waLBerla-wind does not support buffer processes for the moment. Please don't use more processes than you have blocks.")
            WALBERLA_LOG_INFO_ON_ROOT("SetupBlockForest created successfully:\n" << *forest)

            auto bf = std::make_shared<walberla::blockforest::BlockForest>(uint_t(walberla::MPIManager::instance()->rank()), *forest, true);

            if(snapshotHandler_->storeSnapshot && !snapshotHandler_->loadSnapshot) {
                bf->saveToFile(blockforestFileName);
            }

            auto sbf = std::make_shared<walberla::blockforest::StructuredBlockForest>(bf, cellsPerBlock_[0],
                                                                                      cellsPerBlock_[1],
                                                                                      cellsPerBlock_[2]);
            sbf->createCellBoundingBoxes();

            return sbf;

        } // function calculateStructuredBlockForest

        std::shared_ptr< walberla::blockforest::StructuredBlockForest >
        DomainSetup::createUniformBlockForest( const std::vector<AABB> & turbineAABBs ) {

            std::string blockforestFileName;
            if(snapshotHandler_->loadSnapshot || snapshotHandler_->storeSnapshot) {
                blockforestFileName = std::string(walberla::filesystem::path(snapshotHandler_->baseFolder) / walberla::filesystem::path(snapshotHandler_->blockforestFile).c_str());
            }

            if(snapshotHandler_->loadSnapshot) {

                WALBERLA_LOG_INFO_ON_ROOT("Reading the block structure from file " << snapshotHandler_->blockforestFile << "...")

                // load block forest from file
                walberla::MPIManager::instance()->useWorldComm();

                auto bf = std::make_shared< walberla::BlockForest >(
                        uint_t(walberla::MPIManager::instance()->rank()),
                        blockforestFileName.c_str(),
                        true, false);

                // save data from blockforest for domain decomposition calculation
                const auto & aabb = bf->getDomain();
                for(uint_t d = 0; d < 3; ++d) {
                    periodicity_[d] = bf->isPeriodic(d);
                    domainSize_[d] = uint_t(aabb.max(d));
                    numberOfCoarseBlocksPerDirection_[d] = bf->getSize(d);
                    cellsPerBlock_[d] = domainSize_[d] / numberOfCoarseBlocksPerDirection_[d];
                }

                fixedCellsPerBlock_ = true;

            }

            WALBERLA_CHECK( numLevels_ == uint_t(0), "Cannot use numLevel != 0 for uniform blockforest." )

            auto nrOfProcesses = uint_t( walberla::MPIManager::instance()->numProcesses() );

            if( !fixedCellsPerBlock_ ) {
                walberla::blockforest::calculateCellDistribution(domainSize_, nrOfProcesses,
                                                                 numberOfCoarseBlocksPerDirection_, cellsPerBlock_);
                for( uint_t d = 0; d < 3; ++d ) {
                    domainSize_[d] = numberOfCoarseBlocksPerDirection_[d] * cellsPerBlock_[d];
                }
            } else {
                const auto changedSize = adjustDomainSizeToBlockSize();

                if( changedSize ) {
                    WALBERLA_LOG_WARNING_ON_ROOT("Domain size was changed to fit the number of cells per Block. New domain size : " << domainSize_)
                }
            }

            auto domainAABB = walberla::AABB(real_t(0), real_t(0), real_t(0),
                                             real_t(domainSize_[0]), real_t(domainSize_[1]), real_t(domainSize_[2]));

            walberla::SetupBlockForest sforest;

            const walberla::memory_t memoryPerBlock = walberla::numeric_cast<walberla::memory_t>(
                    (cellsPerBlock_[0] + uint_t(2) * uint_t(1))  *
                    (cellsPerBlock_[1] + uint_t(2) * uint_t(1))  *
                    (cellsPerBlock_[2] + uint_t(2) * uint_t(1))) * memoryPerCell_;

            sforest.addWorkloadMemorySUIDAssignmentFunction( std::bind(workloadAndMemoryAssignment, std::placeholders::_1, memoryPerBlock, turbineAABBs ) );

            sforest.init( domainAABB, numberOfCoarseBlocksPerDirection_[0], numberOfCoarseBlocksPerDirection_[1], numberOfCoarseBlocksPerDirection_[2],
                          periodicity_[0], periodicity_[1], periodicity_[2] );

            sforest.balanceLoad(walberla::blockforest::StaticLevelwiseCurveBalanceWeighted(true), walberla::uint_c( walberla::MPIManager::instance()->numProcesses() ), real_t(0),
                                processMemoryLimit_, false);

            if( !walberla::MPIManager::instance()->rankValid() )
                walberla::MPIManager::instance()->useWorldComm();

            WALBERLA_ROOT_SECTION() {
                if (outputSetupForest_) {
                    sforest.writeVTKOutput("initialDomainDecomposition");
                    sforest.writeCSV("process_distribution");
                }
            }

            WALBERLA_CHECK_EQUAL(sforest.getNumberOfBufferProcesses(), 0, "waLBerla-wind does not support buffer processes for the moment. Please don't use more processes than you have blocks.")
            WALBERLA_LOG_INFO_ON_ROOT("SetupBlockForest created successfully:\n" << sforest)

            // create StructuredBlockForest (encapsulates a newly created BlockForest)

            auto bf = std::make_shared< walberla::BlockForest >( walberla::uint_c( walberla::MPIManager::instance()->rank() ), sforest, true );

            if(snapshotHandler_->storeSnapshot && !snapshotHandler_->loadSnapshot) {
                bf->saveToFile(blockforestFileName);
            }

            auto forest = std::make_shared< walberla::StructuredBlockForest >( bf, cellsPerBlock_[0], cellsPerBlock_[1], cellsPerBlock_[2] );
            forest->createCellBoundingBoxes();

            return forest;

        } // function createUniformBlockForest

    } // domain

} // turbine_core 

#endif // TURBINECORE_DOMAINSETUP_H
