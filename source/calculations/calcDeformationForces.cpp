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

#include <Eigen/Dense>
#include <vector>

#include "calcDeformationForces.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../functions/zeroForces.hpp"
#include "../configuration/core_config.h"

std::vector<std::vector<std::pair<int, int>>>
getCorrespondingTrianglesForNodes(const std::vector<Triangle> &triangles, const std::vector<Node> &nodes) {
    std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes(nodes.size());
    for (int i = 0; i < triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {
            correspondingTrianglesForNodes[triangles[i].vertexLabels(j)].emplace_back(i, j);
            correspondingTrianglesForNodes[triangles[i].nonVertexPatchNodesLabels(j)].emplace_back(i, j + 3);
        }
    }
    return correspondingTrianglesForNodes;
}
//
//void update_elastic_forces(
//        std::vector<Node> &nodes,
//        std::vector<Triangle> &triangles,
//        const CoreConfig &core_config,
//        const std::vector<std::vector<std::pair<int, int>>> &correspondingTrianglesForNodes) {
//
//    double stretchingPreFac = 0.5 * core_config.getThickness() * core_config.getShearModulus();
//    std::vector<Eigen::Vector3d> forcesForEachTriangle(6 * triangles.size());
//#pragma omp parallel
//    {
//    #pragma omp for
//        for (int i = 0; i < triangles.size(); i++) {
//            triangles[i].updateHalfPK1Stress(stretchingPreFac);
//            Eigen::Matrix<double, 3, 3> stretchForces = triangles[i].getStretchingForces();
//
//            Eigen::Matrix<double, 3, 3> triangleEdgeNormals = triangles[i].getTriangleEdgeNormals();
//            Eigen::Matrix<double, 3, 3> normalDerivPiece =
//                    0.5 * triangles[i].currAreaInv * (triangles[i].patchSecDerivs.transpose() * triangleEdgeNormals);
//
//            for (int n = 0; n < 3; ++n) {
//                forcesForEachTriangle[6 * i + n] = triangles[i].getBendingForce(normalDerivPiece, n) + stretchForces.col(n);
//            }
//            for (int n = 3; n < 6; ++n) {
//                forcesForEachTriangle[6 * i + n] = triangles[i].getBendingForce(normalDerivPiece, n);
//            }
//        }
//
//    #pragma omp for
//        for (int i = 0; i < nodes.size(); i++) {
//            for (auto &trianglesForNode: correspondingTrianglesForNodes[i]) {
//                int index = 6 * trianglesForNode.first + trianglesForNode.second;
//                nodes[i].force += forcesForEachTriangle[index];
//            }
//        }
//    }
//}