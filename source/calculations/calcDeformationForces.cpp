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

Function to calculate the current force due to strain and bending on each
*node*.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <vector>

#include "calcDeformationForces.hpp"


void assignForceLocationsToNodes(std::vector<Triangle> &triangles, std::vector<Node> &nodes) {
    for (int i = 0; i < triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {
            nodes[triangles[i].vertexLabels(j)].addForcePoints(triangles[i].getNodeForce(j));
            nodes[triangles[i].nonVertexPatchNodesLabels(j)].addForcePoints(triangles[i].getNodeForce(j + 3));
        }
    }
}