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

#include "../Eigen/Dense"
#include <cstddef>
#include <string>
#include <vector>

#include "Node.hpp"
#include "Triangle.hpp"
#include "SettingsStruct.hpp"

void writeVTKDataOutput(
        const std::vector<Node> &,
        const std::vector<Triangle> &,
        const int &,
        const double &,
        const double &,
        const size_t &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<double> &,
        const std::vector<Eigen::Vector2d> &,
        const std::vector<Eigen::Matrix<double, 3, 2> > &,
        const SettingsStruct &,
        const std::string &);
