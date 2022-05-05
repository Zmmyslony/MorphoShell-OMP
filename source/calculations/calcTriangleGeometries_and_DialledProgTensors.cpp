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
#include <omp.h>

#include "calcTriangleGeometries_and_DialledProgTensors.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../SimulationStatus.hpp"
#include "../SettingsStruct.hpp"


void updateMetrics(std::vector<Triangle> &triangles, const std::vector<Node> &nodes) {
    omp_set_num_threads(8);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < triangles.size(); ++i) {
        triangles[i].updateMetric(nodes);
    }
}


void updateDialledInverseProgrammedMetricForATriangleWithoutAnsatz(int index, std::vector<Triangle> &triangles,
                                                                   double currDialInFactor,
                                                                   const Eigen::Vector3d &metricPrevious,
                                                                   const Eigen::Vector3d &metricNew) {
    double dirAng = (1.0 - currDialInFactor) * metricPrevious(0) + currDialInFactor * metricNew(0);
    double cosDirAng = cos(dirAng);
    double sinDirAng = sin(dirAng);
    double lambda = (1.0 - currDialInFactor) * metricPrevious(1) + currDialInFactor * metricNew(1);
    double nu = (1.0 - currDialInFactor) * metricPrevious(2) + currDialInFactor * metricNew(2);
    double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
    double lambdaToThe2Nu = pow(lambda, 2.0 * nu);

    triangles[index].dialledInvProgMetric(0, 0) =
            lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
    triangles[index].dialledInvProgMetric(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
    triangles[index].dialledInvProgMetric(1, 0) = triangles[index].dialledInvProgMetric(0, 1);
    triangles[index].dialledInvProgMetric(1, 1) =
            lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
}


void updateDialledInverseProgrammedMetricForATriangleWithoutAnsatz(int index, std::vector<Triangle> &triangles,
                                                                   double currDialInFactor,
                                                                   const size_t &progTensorSequenceCounter,
                                                                   const std::vector<std::vector<Eigen::Vector3d>> &programmedMetricInfoSequence) {
    Eigen::Vector3d metricPrevious = programmedMetricInfoSequence[progTensorSequenceCounter][index];
    Eigen::Vector3d metricNew = programmedMetricInfoSequence[progTensorSequenceCounter + 1][index];

    updateDialledInverseProgrammedMetricForATriangleWithoutAnsatz(index, triangles, currDialInFactor, metricPrevious,
                                                                  metricNew);
}

void updateDialledInverseProgrammedMetricFromAnsatz(int index, std::vector<Triangle> &triangles,
                                                    double currDialInFactor,
                                                    const size_t &progTensorSequenceCounter,
                                                    const std::vector<std::vector<Eigen::Matrix<double, 2, 2>>> &programmedMetricSequence) {
    triangles[index].dialledInvProgMetric =
            (1.0 - currDialInFactor) * programmedMetricSequence[progTensorSequenceCounter][index] +
            currDialInFactor * programmedMetricSequence[progTensorSequenceCounter + 1][index];

}


void updateDialledSecondFundamantalForm(int index, std::vector<Triangle> &triangles, double currDialInFactor,
                                        double rootCurrDialInFactor, const size_t &progTensorSequenceCounter,
                                        const std::vector<std::vector<double> > &programmedTausSequence,
                                        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedSecondFundamentalFormSequence) {
    triangles[index].dialledProgTau =
            (1.0 - currDialInFactor) * programmedTausSequence[progTensorSequenceCounter][index] +
            currDialInFactor * programmedTausSequence[progTensorSequenceCounter + 1][index];
    triangles[index].dialledProgSecFF =
            (1.0 - rootCurrDialInFactor) *
            programmedSecondFundamentalFormSequence[progTensorSequenceCounter][index] +
            rootCurrDialInFactor * programmedSecondFundamentalFormSequence[progTensorSequenceCounter + 1][index];

}

void
updateDeterminantOfDialledInverseProgrammedMetric(int index, std::vector<Triangle> &triangles) {
    triangles[index].detDialledInvProgMetric = triangles[index].dialledInvProgMetric.determinant();
}

void updateGeometricPropertiesOfAllTriangles(std::vector<Triangle> &triangles, const std::vector<Node> &nodes) {
    omp_set_num_threads(8);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].updateGeometricProperties(nodes);
    }
}


void calcTriangleGeometries_and_DialledProgTensors(
        const std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SimulationStatus &status,
        const double &currDialInFactor,
        const size_t &progTensorSequenceCounter,
        const std::vector<std::vector<Eigen::Vector3d> > &programmedMetricInfoSequence,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedMetricSequence,
        const std::vector<std::vector<double> > &programmedTausSequence,
        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedSecondFundamentalFormSequence,
        const SettingsStruct &settings) {

    // Will use some functions of dial-in factor, pre-calculated here.
    const double rootCurrDialInFactor = sqrt(currDialInFactor);

    updateMetrics(triangles, nodes);

    omp_set_num_threads(8);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < triangles.size(); i++) {
        if (status == dialling) {
            if (settings.isLCEModeEnabled && !settings.isDialingFromAnsatzEnabled) {
                updateDialledInverseProgrammedMetricForATriangleWithoutAnsatz(i, triangles, currDialInFactor,
                                                                                 progTensorSequenceCounter,
                                                                                 programmedMetricInfoSequence);

            } else {
                updateDialledInverseProgrammedMetricFromAnsatz(i, triangles, currDialInFactor,
                                                               progTensorSequenceCounter,
                                                               programmedMetricSequence);
            }
            updateDialledSecondFundamantalForm(i, triangles, currDialInFactor, rootCurrDialInFactor,
                                               progTensorSequenceCounter, programmedTausSequence,
                                               programmedSecondFundamentalFormSequence);
            updateDeterminantOfDialledInverseProgrammedMetric(i, triangles);
        }
    }
    updateGeometricPropertiesOfAllTriangles(triangles, nodes);
}




//    omp_set_num_threads(8);
//#pragma omp parallel for schedule(static)
//    for (int i = 0; i < settings.NumTriangles; ++i) {
//        /* If dialling in is in progress (as opposed to waiting for equilibrium),
//        update the programmed tensors for each triangle.
//        The square root is just to make the dialling in a bit more even - using
//        just currDialInFactor leads to slow growth at the start which accelerates.
//        This is to be expected since the programmed metric is akin to the square
//        of the 'programmed' (up to rotations) deformation gradient. */
//        if (status == dialling) {
//
//            // Usual LCE (inverse) metric (not possible to do this if dialling
//            // from an ansatz).
//            if (settings.isLCEModeEnabled && !settings.isDialingFromAnsatzEnabled) {
//                double dirAng =
//                        (1.0 - currDialInFactor) * programmedMetricInfoSequence[progTensorSequenceCounter][i](0) +
//                        currDialInFactor * programmedMetricInfoSequence[progTensorSequenceCounter + 1][i](0);
//                double cosDirAng = cos(dirAng);
//                double sinDirAng = sin(dirAng);
//                double lambda =
//                        (1.0 - currDialInFactor) * programmedMetricInfoSequence[progTensorSequenceCounter][i](1) +
//                        currDialInFactor * programmedMetricInfoSequence[progTensorSequenceCounter + 1][i](1);
//                double nu = (1.0 - currDialInFactor) * programmedMetricInfoSequence[progTensorSequenceCounter][i](2) +
//                            currDialInFactor * programmedMetricInfoSequence[progTensorSequenceCounter + 1][i](2);
//                double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
//                double lambdaToThe2Nu = pow(lambda, 2.0 * nu);
//
//                triangles[i].dialledInvProgMetric(0, 0) =
//                        lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
//                triangles[i].dialledInvProgMetric(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
//                triangles[i].dialledInvProgMetric(1, 0) = triangles[i].dialledInvProgMetric(0, 1);
//                triangles[i].dialledInvProgMetric(1, 1) =
//                        lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
//            }
//                // Dialling in of more general metric.
//            else {
//                triangles[i].dialledInvProgMetric =
//                        (1.0 - currDialInFactor) * programmedMetricSequence[progTensorSequenceCounter][i] +
//                        currDialInFactor * programmedMetricSequence[progTensorSequenceCounter + 1][i];
//            }
//
//
//            triangles[i].dialledProgTau =
//                    (1.0 - currDialInFactor) * programmedTausSequence[progTensorSequenceCounter][i] +
//                    currDialInFactor * programmedTausSequence[progTensorSequenceCounter + 1][i];
//            triangles[i].dialledProgSecFF =
//                    (1.0 - rootCurrDialInFactor) *
//                    programmedSecondFundamentalFormSequence[progTensorSequenceCounter][i] +
//                    rootCurrDialInFactor * programmedSecondFundamentalFormSequence[progTensorSequenceCounter + 1][i];
//            triangles[i].detDialledInvProgMetric = triangles[i].dialledInvProgMetric.determinant();
//        }
//
//
//        triangles[i].updateGeometricProperties(nodes);
//    }
//}
