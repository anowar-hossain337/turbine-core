
#pragma once

#ifndef TURBINECORE_FILEPARSER_H
#define TURBINECORE_FILEPARSER_H

#include <string>
#include <memory>
#include <vector>

#include <fstream>

#include <core/math/Constants.h>

#include "data_interpolator/DataInterpolator.h"
#include "aerodynamics/AirfoilPolar.h"

#include "conversion/Conversion.h"
#include "wind_turbine_core/ProjectDefines.h"
#include "wind_turbine_core/math/Vector3.h"

namespace turbine_core {

    namespace creator {

        template<bool CUBIC>
        HOST_PREFIX void parseDescriptionFile( const std::string & filename, const uint_t & nPoints,
                                               const std::shared_ptr<Vector3DataInterpolator<CUBIC>> & positionInterpolator,
                                               const std::shared_ptr<Vector3DataInterpolator<CUBIC>> & angleInterpolator,
                                               const std::shared_ptr<ScalarDataInterpolator<CUBIC>> & widthInterpolator,
                                               std::vector<aerodynamics::AirfoilPolar> & airfoilPolars ) {

            std::vector<Vector3<real_t>> positions, angles;
            std::vector<real_t> width;

            std::string curLine;
            std::ifstream inputFile(filename.c_str());

            if (!inputFile.is_open()) {
                std::cout << "Error opening file " << filename
                          << ", please check file existence." << std::endl;
            } else {
                std::getline(inputFile, curLine); // skip header

                // get first line
                std::getline(inputFile, curLine);
                std::stringstream firstIn(curLine);

                real_t tmpX, tmpY, tmpZ, tmpXA, tmpYA, tmpZA, tmpW;
                std::string tmpAirfoil;

                firstIn >> tmpX >> tmpY >> tmpZ >> tmpXA >> tmpYA >> tmpZA >> tmpW >> tmpAirfoil;

                positions.emplace_back(Vector3<real_t>(tmpX, tmpY, tmpZ) / Conversion::C_l());

                angles.emplace_back(Vector3<real_t>(tmpXA, tmpYA, tmpZA) / real_t(180.0) * walberla::math::pi);

                width.push_back(tmpW / Conversion::C_l());

                airfoilPolars.emplace_back(tmpAirfoil);

                std::string lastLine;
                while( std::getline(inputFile, curLine) ) {
                    bool is_empty = true;
                    for (char ch : curLine) {
                        is_empty = is_empty && isspace(ch);
                    }
                    if (!is_empty) {
                        lastLine = curLine;
                        if (nPoints) {
                            std::stringstream inStream(curLine);

                            inStream >> tmpX >> tmpY >> tmpZ >> tmpXA >> tmpYA >> tmpZA >> tmpW >> tmpAirfoil;

                            positions.emplace_back(Vector3<real_t>(tmpX, tmpY, tmpZ) / Conversion::C_l());

                            auto tmpAngles = Vector3<real_t>(tmpXA, tmpYA, tmpZA) / real_t(180.0) * walberla::math::pi;
                            angles.push_back(tmpAngles);

                            width.push_back(tmpW / Conversion::C_l());

                            airfoilPolars.emplace_back(tmpAirfoil);
                        }
                    }

                } // end while

                if(!nPoints) {
                    std::stringstream lastIn(lastLine);

                    lastIn >> tmpX >> tmpY >> tmpZ >> tmpXA >> tmpYA >> tmpZA >> tmpW >> tmpAirfoil;

                    positions.emplace_back(Vector3<real_t>(tmpX, tmpY, tmpZ) / Conversion::C_l());

                    auto tmpAngles = Vector3<real_t>(tmpXA, tmpYA, tmpZA) / real_t(180.0) * walberla::math::pi;
                    angles.push_back(tmpAngles);

                    width.push_back(tmpW / Conversion::C_l());

                    airfoilPolars.emplace_back(tmpAirfoil);
                }
            }

            const uint_t interpolationPoints_ = positions.size();

            std::vector<real_t> span;
            span.reserve(interpolationPoints_);

            span.push_back(real_t(0.0));

            for(uint_t i = 1; i < interpolationPoints_; ++i) {
                span.push_back( (positions[i] - positions[i-1]).length() + span.back() );
            }

            positionInterpolator->setPoints(interpolationPoints_, span.data(), positions.data());
            angleInterpolator->setPoints(interpolationPoints_, span.data(), angles.data());
            widthInterpolator->setPoints(interpolationPoints_, span.data(), width.data());

        }

    }

}

#endif //TURBINECORE_FILEPARSER_H
