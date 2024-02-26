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

Function to loop over triangular elements, and calculate vectors
describing two of the current sides, the current normal to the face, and the
current area. The current 'dialled in' programmed metric and second fundamental
form for each triangle are also calculated.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <cstddef>
#include <Eigen/Dense>
#include <vector>
#include <cmath>

#include "calcTriangleGeometries_and_DialledProgTensors.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../settings_new.h"

void
updateTriangleProperties(const std::vector<Node> &nodes, std::vector<Triangle> &triangles, const SimulationStatus &status,
                         const double &dial_in_factor, const size_t &stage_counter, const SettingsNew &settings) {

    // Will use some functions of dial-in factor, pre-calculated here.
    const double dial_in_factor_root = sqrt(dial_in_factor);
    bool is_LCE_metric_used = settings.getCore().isLceModeEnabled() && !settings.getCore().isAnsatzMetricUsed();


#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        if (status == Dialling) {
            triangles[i].updateProgrammedQuantities(stage_counter, dial_in_factor, dial_in_factor_root, is_LCE_metric_used);
        }
        triangles[i].updateGeometricProperties(nodes);
    }

}