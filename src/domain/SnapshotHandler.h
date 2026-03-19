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
//! \file SnapshotHandler.h
//! \author Helen SCHOTTENHAMML <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#ifndef TURBINECORE_SNAPSHOTHANDLER_H
#define TURBINECORE_SNAPSHOTHANDLER_H

#pragma once

#include <blockforest/SetupBlockForest.h>
#include <blockforest/StructuredBlockForest.h>
#include <core/config/Config.h>
#include <core/mpi/Broadcast.h>

#include <vtk/VTKOutput.h>

#include "wind_turbine_core/WalberlaDataTypes.h"

namespace turbine_core {

    namespace domain {

        struct SnapshotHandler {

            explicit SnapshotHandler(std::shared_ptr<walberla::Config> & globalConfig) {

                const auto config = globalConfig->getOneBlock("DomainSetup");

                // load parameters
                loadSnapshot           = config.getParameter< bool >("loadSnapshot", false);
                storeSnapshot          = config.getParameter< bool >("storeSnapshot", false);
                snapshotFrequency      = config.getParameter< uint_t >("snapshotFrequency", uint_t(1));
                if(storeSnapshot) {
                    WALBERLA_CHECK(snapshotFrequency > uint_t(0), "Snapshot saving frequency \"snapshotFrequency\" must be greater than 0, "
                                                                  "please check your input file.")
                }

                baseFolder             = config.getParameter< std::string >("snapshotBaseFolder", "snapshot");
                configFile             = config.getParameter< std::string >("snapshotConfigFile", "");
                welfordFile            = config.getParameter< std::string >("welfordFile", "");
                blockforestFile        = config.getParameter< std::string >("snapshotBlockforestFile", "");
                pdfFile                = config.getParameter< std::string >("snapshotPdfFile", "");
                forceFile              = config.getParameter< std::string >("snapshotForceFile", "");
                meanVelFile            = config.getParameter< std::string >("snapshotMeanVelFile", "");
                sosFile                = config.getParameter< std::string >("snapshotSOSFile", "");
                socFile                = config.getParameter< std::string >("snapshotSOCFile", "");

                // early out
                if(!loadSnapshot && !storeSnapshot)
                    return;

                if(storeSnapshot) {
                    // create directories if they do not exist already
                    WALBERLA_ROOT_SECTION() {
                        walberla::filesystem::path snapshotPath(baseFolder);
                        if (!walberla::filesystem::exists(snapshotPath)) {
                            walberla::filesystem::create_directory(snapshotPath);
                        }
                    }
                }

                if(loadSnapshot) {
                    // read snapshot config file
                    WALBERLA_ROOT_SECTION()
                    {
                        auto filePath = walberla::filesystem::path(baseFolder) / walberla::filesystem::path(configFile);
                        {
                            std::ifstream file(filePath);
                            if (!file) {
                                WALBERLA_ABORT("Error: " << configFile << " could not be opened!")
                            } else {
                                file >> startingTimestep;
                            }
                        }

                        filePath = walberla::filesystem::path(baseFolder) / walberla::filesystem::path(welfordFile);
                        {
                            std::ifstream file(filePath);
                            if (!file) {
                                WALBERLA_ABORT("Error: " << welfordFile << " could not be opened!")
                            } else {
                                file >> welfordCounter;
                            }
                        }
                    }
                    walberla::mpi::broadcastObject(startingTimestep);
                    walberla::mpi::broadcastObject(welfordCounter);

                    WALBERLA_LOG_INFO_ON_ROOT("Successfully read config parameters from snapshot config file:")
                    WALBERLA_LOG_INFO_ON_ROOT(" - startingTimestep = " << startingTimestep)

                    WALBERLA_LOG_INFO_ON_ROOT("Successfully read welford counter from snapshot config file:")
                    WALBERLA_LOG_INFO_ON_ROOT(" - welfordCounter = " << welfordCounter)

                    // modify VTK output to start with output number xxx instead of 0
                    std::vector< walberla::config::Config::Block* > configOutputBlock;
                    globalConfig->getWritableGlobalBlock().getWritableBlocks("Output", configOutputBlock, 1, 1);
                    std::vector< walberla::config::Config::Block* > configVTKBlock;
                    configOutputBlock[0]->getWritableBlocks("VTK", configVTKBlock, 1, 1);
                    configVTKBlock[0]->setOrAddParameter("initialExecutionCount", std::to_string(startingTimestep));

                }

            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID) {
                if(storeSnapshot) {
                    if(t == 0){
                        if(pdfFile.empty()) WALBERLA_ABORT("No pdf filename given for snapshot creation.")
                        if(forceFile.empty()) WALBERLA_ABORT("No force filename given for snapshot creation.")
                    }

                    const auto pdfFileName = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(pdfFile));
                    const auto forceFileName = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(forceFile));

                    if (t % snapshotFrequency == 0 && t > 0) {
                        blocks->saveBlockData(pdfFileName, pdfFieldID);
                        blocks->saveBlockData(forceFileName, forceFieldID);
                    }

                    WALBERLA_ROOT_SECTION()
                    {
                        const auto checkpointConfigFile = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(configFile));
                        std::ofstream file;
                        file.open(checkpointConfigFile.c_str());

                        if(file.is_open()) {
                            file << std::setprecision(16);
                            file << t << "\n";
                            file.close();
                        }
                    }

                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID ) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID);

                if(storeSnapshot) {
                    if(t == 0 && meanVelFile.empty()){
                        WALBERLA_ABORT("No mean velocity filename given for snapshot creation.")
                    }
                    const auto meanVelocityFileName = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(meanVelFile));

                    if (t % snapshotFrequency == 0 && t > 0) {
                        blocks->saveBlockData(meanVelocityFileName, meanVelocityFieldID);
                    }

                    WALBERLA_ROOT_SECTION()
                    {
                        std::ofstream file;

                        const auto welfordConfigFile = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(welfordFile));
                        file.open(welfordConfigFile.c_str());

                        if(file.is_open()) {
                            file << std::setprecision(16);
                            file << welfordCounter-2 << "\n";
                            file.close();
                        }
                    }

                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID,
                                const walberla::BlockDataID & sosFieldID) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID, welfordCounter, meanVelocityFieldID);

                if(storeSnapshot) {
                    // check for valid file names
                    if(t == 0 && sosFile.empty()){
                        WALBERLA_ABORT("No SoS filename given for snapshot creation.")
                    }

                    const auto sosFileName = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(sosFile));

                    if (t % snapshotFrequency == 0 && t > 0) {
                        blocks->saveBlockData(sosFileName, sosFieldID);
                    }
                }
            }

            void createSnapshot(const uint_t t, const std::shared_ptr<walberla::StructuredBlockForest> & blocks,
                                const walberla::BlockDataID & pdfFieldID, const walberla::BlockDataID & forceFieldID,
                                const uint_t welfordCounter,
                                const walberla::BlockDataID & meanVelocityFieldID,
                                const walberla::BlockDataID & sosFieldID,
                                const walberla::BlockDataID & socFieldID) {
                createSnapshot(t, blocks, pdfFieldID, forceFieldID, welfordCounter, meanVelocityFieldID, sosFieldID);

                if(storeSnapshot) {
                    // check for valid file names
                    if(t == 0 && socFile.empty()){
                        WALBERLA_ABORT("No SoC filename given for snapshot creation.")
                    }

                    const auto socFileName = std::string(walberla::filesystem::path(baseFolder) / walberla::filesystem::path(socFile));

                    if (t % snapshotFrequency == 0 && t > 0) {
                        blocks->saveBlockData(socFileName, socFieldID);
                    }
                }
            }

            bool loadSnapshot;
            bool storeSnapshot;

            uint_t snapshotFrequency;
            uint_t startingTimestep;
            uint_t welfordCounter = 0;

            std::string baseFolder;
            std::string configFile;
            std::string blockforestFile;
            std::string pdfFile;
            std::string forceFile;
            std::string meanVelFile;
            std::string sosFile;
            std::string socFile;
            std::string welfordFile;
        };

    } // domain

} // turbine_core 

#endif // TURBINECORE_SNAPSHOTHANDLER_H
