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

#include "../../Eigen/Dense"
#include <vector>
#include <omp.h>

#include "calcDeformationForces.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../SettingsStruct.hpp"
#include "../functions/zeroForces.hpp"

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

void calcDeformationForces(
        std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings,
        const std::vector<std::vector<std::pair<int, int>>> &correspondingTrianglesForNodes) {

    double stretchingPreFac = 0.5 * settings.SheetThickness * settings.ShearModulus;
    std::vector<std::vector<Eigen::Vector3d>>>
    omp_set_num_threads(8);
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].updateHalfPK1Stress(stretchingPreFac);
        Eigen::Matrix<double, 3, 3> stretchForces = triangles[i].getStretchingForces();

        for (int v = 0; v < 3; ++v) {
            nodes[triangles[i].vertexLabels(v)].force += stretchForces.col(v);
        }

        Eigen::Matrix<double, 3, 3> outwardTriNormals = triangles[i].getOutwardTriangleNormals();
        Eigen::Matrix<double, 3, 3> normalDerivPiece =
                0.5 * triangles[i].invCurrArea * (triangles[i].patchSecDerivs.transpose() * outwardTriNormals);

        for (int n = 0; n < 3; ++n) {
            nodes[triangles[i].vertexLabels(n)].force += triangles[i].getBendingForce(normalDerivPiece, n);
            nodes[triangles[i].nonVertexPatchNodesLabels(n)].force += triangles[i].getBendingForce(normalDerivPiece, n + 3);
        }
    }
}

void calcDeformationForces(
        std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings) {

    double stretchingPreFac = 0.5 * settings.SheetThickness * settings.ShearModulus;

    for (auto &triangle: triangles) {
        triangle.updateHalfPK1Stress(stretchingPreFac);
        Eigen::Matrix<double, 3, 3> stretchForces = triangle.getStretchingForces();

        for (int v = 0; v < 3; ++v) {
            nodes[triangle.vertexLabels(v)].force += stretchForces.col(v);
        }

        Eigen::Matrix<double, 3, 3> outwardTriNormals = triangle.getOutwardTriangleNormals();
        Eigen::Matrix<double, 3, 3> normalDerivPiece =
                0.5 * triangle.invCurrArea * (triangle.patchSecDerivs.transpose() * outwardTriNormals);

        for (int n = 0; n < 3; ++n) {
            nodes[triangle.vertexLabels(n)].force += triangle.getBendingForce(normalDerivPiece, n);
            nodes[triangle.nonVertexPatchNodesLabels(n)].force += triangle.getBendingForce(normalDerivPiece, n + 3);
        }
    }
}