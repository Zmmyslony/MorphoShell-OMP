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
#include "../../Eigen/Dense"
#include <vector>
#include <cmath>
#include <omp.h>

#include "calcTriangleGeometries_and_DialledProgTensors.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../StatusEnum.hpp"
#include "../SettingsStruct.hpp"

void calcTriangleGeometries_and_DialledProgTensors(
        const std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const StatusEnum &status,
        const double &currDialInFactor,
        const size_t &progTensorSequenceCounter,
        const std::vector<std::vector<Eigen::Vector3d> > &sequenceOf_ProgMetricInfo,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_InvProgMetrics,
        const std::vector<std::vector<double> > &sequenceOf_ProgTaus,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_ProgSecFFs,
        const SettingsStruct &settings) {

    // Will use some functions of dial-in factor, pre-calculated here.
    const double rootCurrDialInFactor = sqrt(currDialInFactor);

    omp_set_num_threads(8);
#pragma omp parallel for
    for (int i = 0; i < settings.NumTriangles; ++i) {

        /*set the first triangle side (col(0)) to be the vector
        from the first vertex to the second*/
        triangles[i].currSides.col(0) =
                nodes[triangles[i].vertexLabels(1)].pos - nodes[triangles[i].vertexLabels(0)].pos;

        /*set the second triangle side to be the vector
        from the first vertex to the third*/
        triangles[i].currSides.col(1) =
                nodes[triangles[i].vertexLabels(2)].pos - nodes[triangles[i].vertexLabels(0)].pos;


        /* Calculate deformation gradient (3x2 matrix).
        Since we model deformation as constant across each triangle, the
        deformation gradient of each triangle is given by multiplying the
        triangleSides matrix by a sequence of deformation gradient matrices.
        Thus:
        triangleSides(current) = defGradient * triangleSides_init
        Rearranging:
        defGradient = triangleSides(current) * inv(triangleSides_init.*/
        triangles[i].defGradient = triangles[i].currSides * triangles[i].invInitSidesMat;


        // Calculate corresponding metric, and its det and inverse.
        triangles[i].metric = triangles[i].defGradient.transpose() * triangles[i].defGradient;
        triangles[i].detInvMetric = 1.0 / (triangles[i].metric(0, 0) * triangles[i].metric(1, 1) -
                                           triangles[i].metric(0, 1) * triangles[i].metric(0, 1));
        Eigen::Matrix<double, 2, 2> metricAdjMatrix;
        metricAdjMatrix << triangles[i].metric(1, 1), -triangles[i].metric(0, 1), -triangles[i].metric(0,1), triangles[i].metric(0, 0);
        triangles[i].invMetric = triangles[i].detInvMetric * metricAdjMatrix;

        /* If dialling in is in progress (as opposed to waiting for equilibrium),
        update the programmed tensors for each triangle.
        The square root is just to make the dialling in a bit more even - using
        just currDialInFactor leads to slow growth at the start which accelerates.
        This is to be expected since the programmed metric is akin to the square
        of the 'programmed' (up to rotations) deformation gradient. */
        if (status == dialling) {

            // Usual LCE (inverse) metric (not possible to do this if dialling
            // from an ansatz).
            if (settings.isLCEModeEnabled && !settings.isDialingFromAnsatzEnabled) {
                double dirAng = (1.0 - currDialInFactor) * sequenceOf_ProgMetricInfo[progTensorSequenceCounter][i](0) +
                                currDialInFactor * sequenceOf_ProgMetricInfo[progTensorSequenceCounter + 1][i](0);
                double cosDirAng = cos(dirAng);
                double sinDirAng = sin(dirAng);
                double lambda = (1.0 - currDialInFactor) * sequenceOf_ProgMetricInfo[progTensorSequenceCounter][i](1) +
                                currDialInFactor * sequenceOf_ProgMetricInfo[progTensorSequenceCounter + 1][i](1);
                double nu = (1.0 - currDialInFactor) * sequenceOf_ProgMetricInfo[progTensorSequenceCounter][i](2) +
                            currDialInFactor * sequenceOf_ProgMetricInfo[progTensorSequenceCounter + 1][i](2);
                double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
                double lambdaToThe2Nu = pow(lambda, 2.0 * nu);

                triangles[i].dialledInvProgMetric(0, 0) =
                        lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
                triangles[i].dialledInvProgMetric(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
                triangles[i].dialledInvProgMetric(1, 0) = triangles[i].dialledInvProgMetric(0, 1);
                triangles[i].dialledInvProgMetric(1, 1) =
                        lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
            }
                // Dialling in of more general metric.
            else {
                triangles[i].dialledInvProgMetric =
                        (1.0 - currDialInFactor) * sequenceOf_InvProgMetrics[progTensorSequenceCounter][i] +
                        currDialInFactor * sequenceOf_InvProgMetrics[progTensorSequenceCounter + 1][i];
            }


            triangles[i].dialledProgTau = (1.0 - currDialInFactor) * sequenceOf_ProgTaus[progTensorSequenceCounter][i] +
                                          currDialInFactor * sequenceOf_ProgTaus[progTensorSequenceCounter + 1][i];
            triangles[i].dialledProgSecFF =
                    (1.0 - rootCurrDialInFactor) * sequenceOf_ProgSecFFs[progTensorSequenceCounter][i] +
                    rootCurrDialInFactor * sequenceOf_ProgSecFFs[progTensorSequenceCounter + 1][i];
            triangles[i].detDialledInvProgMetric = triangles[i].dialledInvProgMetric.determinant();
        }


        Eigen::Matrix<double, 3, 6> matrixOfPatchNodeCoords;
        // Loop over patch nodes for this triangle and get their positions.
        for (int n = 0; n < 6; ++n) {
            if (n < 3) {
                matrixOfPatchNodeCoords.col(n) = nodes[triangles[i].vertexLabels(n)].pos;
            } else {
                matrixOfPatchNodeCoords.col(n) = nodes[triangles[i].nonVertexPatchNodesLabels(n - 3)].pos;
            }
        }

        /* Calculate the second spatial derivatives of the patch and some other
        quantities.*/

        triangles[i].patchSecDerivs = matrixOfPatchNodeCoords * triangles[i].matForPatchSecDerivs;

        /* Triangle face normal, using cross product of first and second sides
        from above. The vertex labelling convention should mean these face
        normals point the same way in the flat reference state. The current face
        area is also calculated. These are both only needed if the simpler secFF
        estimate is used.*/
        triangles[i].faceNormal = triangles[i].currSides.col(0).cross(triangles[i].currSides.col(1));
        triangles[i].invCurrArea = 2.0 / triangles[i].faceNormal.norm();
        triangles[i].faceNormal = 0.5 * triangles[i].faceNormal * triangles[i].invCurrArea; // Normalising
    }
}
