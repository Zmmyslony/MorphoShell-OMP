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

Header file for calcTriangleGeometries_and_DialledProgTensors.cpp function.
*/

#include <cstddef>
#include <vector>
#include "Eigen/Dense"

#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../SimulationStatus.hpp"
#include "../Settings.hpp"

void calcTriangleGeometries_and_DialledProgTensors(
        const std::vector<Node> &,
        std::vector<Triangle> &,
        const SimulationStatus &,
        const double &,
        const size_t &,
        const std::vector<std::vector<Eigen::Vector3d>> &,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2>>> &,
        const std::vector<std::vector<double>> &,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2>>> &,
        const Settings &);
