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

Function to calculate a 3x6 matrix for each triangle that can be stored and
used repeatedly to obtain either the first or second spatial derivatives of the
patch for each triangle during the simulation. See main() for more details. To
do this, for each triangle t, the labels of the 3 out of 6 patch-defining nodes
that are not vertices of t are also calculated and stored, and are also used
later.
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <cstddef>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <string>

#include "calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"

void validatePatchSearch(const std::vector<double> &patch_values, double patch_threshold) {
    auto max_iterator = std::max_element(patch_values.begin(), patch_values.end());
    if (*max_iterator < patch_threshold) { return ;}

    int wrong_triangle_count = 0;
    for (auto patch_value : patch_values) {
        if (patch_value >= patch_threshold) { wrong_triangle_count++ ;}
    }
    long worst_offender_index = std::distance(patch_values.begin(), max_iterator);

    throw std::runtime_error(
            "Patch search failed. \n" +
            std::to_string(wrong_triangle_count) + " triangles exceed the maximum threshold value of " + std::to_string(patch_threshold) + "\n" +
            "with maximal value of  " + std::to_string(*max_iterator) + " for triangle " + std::to_string(worst_offender_index)
            );

}

void createNodePatches(const std::vector<Node> &nodes, std::vector<Triangle> &triangles, double patch_threshold) {
    std::vector<double> patch_values(triangles.size());
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        patch_values[i] = triangles[i].updateMatForPatchDerivs(triangles, nodes);
    }
    validatePatchSearch(patch_values, patch_threshold);
}