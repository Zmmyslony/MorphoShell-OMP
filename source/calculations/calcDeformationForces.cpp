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

void calcDeformationForces(
        std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings) {

    // Energy prefactors that are the same for each triangle.
    double stretchingPreFac = 0.5 * settings.SheetThickness * settings.ShearModulus;

    std::vector<std::vector<Eigen::Vector3d>> listOfStretchingForces(triangles.size(),
                                                                     std::vector<Eigen::Vector3d>(3));
    std::vector<std::vector<Eigen::Vector3d>> listOfBendingForces(triangles.size(),
                                                                  std::vector<Eigen::Vector3d>(
                                                                          6));
    // Loop over triangles.
    omp_set_num_threads(8);
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); ++i) {
        triangles[i].updateHalfPK1Stress(stretchingPreFac);
        Eigen::Matrix<double, 3, 3> stretchForces = triangles[i].getStretchingForces();

        for (int v = 0; v < 3; ++v) {
            listOfStretchingForces[i][v] = stretchForces.col(v);
        }

        Eigen::Matrix<double, 3, 3> outwardTriNormals = triangles[i].getOutwardTriangleNormals();
        Eigen::Matrix<double, 3, 3> normalDerivPiece =
                0.5 * triangles[i].invCurrArea * (triangles[i].patchSecDerivs.transpose() * outwardTriNormals);

        for (int n = 0; n < 6; ++n) {
            listOfBendingForces[i][n] = triangles[i].getBendingForce(normalDerivPiece, n);
        }
    }


    for (int i = 0; i < triangles.size(); i++) {
        for (int n = 0; n < 6; n++) {
            if (n < 3) {
                nodes[triangles[i].vertexLabels(n)].force += listOfStretchingForces[i][n];
                nodes[triangles[i].vertexLabels(n)].force += listOfBendingForces[i][n];
            }
            else {
                nodes[triangles[i].nonVertexPatchNodesLabels(n - 3)].force += listOfBendingForces[i][n];
            }
        }
    }
}