/*
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

Header file for writeVTKDataOutput.cpp function.
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>
#include <cstddef>
#include <string>
#include <vector>

#include "Node.hpp"
#include "Triangle.hpp"
// #include "Settings.hpp"

long long int writeVTKDataOutput(const std::vector<Node> &nodes, const std::vector<Triangle> &triangles, const int &stepcount,
                                 const double &time, const double &currDialInFactor, const size_t &progTensorSequenceCounter,
                                 const std::vector<double> &gaussCurvatures, const std::vector<double> &meanCurvatures,
                                 const std::vector<double> &angleDeficits, const std::vector<double> &interiorNodeAngleDeficits,
                                 const std::vector<double> &boundaryNodeAngleDeficits,
                                 const std::vector<double> &stretchEnergies, const std::vector<double> &bendEnergies,
                                 const std::vector<double> &kineticEnergies, const std::vector<double> &strainMeasures,
                                 const std::vector<Eigen::Vector2d> &cauchyStressEigenvals,
                                 const std::vector<Eigen::Matrix<double, 3, 2> > &cauchyStressEigenvecs,
                                 const SettingsNew &settings, const std::string &outputDirName);
