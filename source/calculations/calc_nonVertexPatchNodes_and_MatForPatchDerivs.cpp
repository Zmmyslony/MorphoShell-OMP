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


void calc_nonVertexPatchNodes_and_MatForPatchDerivs(const std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                                                    double patch_threshold) {

    // Number of boundary triangles that tried multiple candidates for patches
    int invalid_non_boundary_triangle_count = 0;

#pragma omp parallel for reduction (+ : invalid_non_boundary_triangle_count)
    for (int i = 0; i < triangles.size(); i++) {
        invalid_non_boundary_triangle_count += triangles[i].updateMatForPatchDerivs(triangles, nodes, patch_threshold);

    }

    /*Print also the number of non-boundary triangles that had to search through
    multiple possible patch options to find one satisfying the determinant
    condition. It seems unlikely that the determinant will be a small for
    triangles in the interior of a reasonable mesh, so if this number is large,
    that suggests something suspicious. One explanation might be that
    settings.PatchMatrixDimensionlessConditioningThreshold has been set to too
    low a value*/
//    logStream.open();
    std::cout << "Number of non-boundary triangles that had to search through \n" <<
              "multiple possible patch options to find one \nsatisfying the condition number " <<
              "criterion was " << invalid_non_boundary_triangle_count <<
              ", \nwhich should not be a large proportion of the mesh's triangles." << std::endl;
//    logStream.close();
}