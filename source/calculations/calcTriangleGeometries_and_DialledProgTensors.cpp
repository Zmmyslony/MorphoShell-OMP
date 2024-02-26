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


void updateMetrics(std::vector<Triangle> &triangles, const std::vector<Node> &nodes) {
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

    triangles[index].programmedMetInv(0, 0) =
            lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
    triangles[index].programmedMetInv(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
    triangles[index].programmedMetInv(1, 0) = triangles[index].programmedMetInv(0, 1);
    triangles[index].programmedMetInv(1, 1) =
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


void updateDialledInverseProgrammedMetricFromAnsatz(int index,
                                                    std::vector<Triangle> &triangles,
                                                    double currDialInFactor,
                                                    const size_t &progTensorSequenceCounter,
                                                    const std::vector<std::vector<Eigen::Matrix<double, 2, 2>>> &programmedMetricSequence) {
    triangles[index].programmedMetInv =
            (1.0 - currDialInFactor) * programmedMetricSequence[progTensorSequenceCounter][index] +
            currDialInFactor * programmedMetricSequence[progTensorSequenceCounter + 1][index];

}


void updateDialledInverseProgrammedMetricFromAnsatz(Triangle &triangle,
                                                    double currDialInFactor,
                                                    const Eigen::Matrix<double, 2, 2> &programmedMetricCurrentStage,
                                                    const Eigen::Matrix<double, 2, 2> &programmedMetricNextStage) {
    triangle.programmedMetInv =
            (1.0 - currDialInFactor) * programmedMetricCurrentStage + currDialInFactor * programmedMetricNextStage;
}


void updateDialledSecondFundemantalForm(int index, std::vector<Triangle> &triangles, double currDialInFactor,
                                        double rootCurrDialInFactor, const size_t &progTensorSequenceCounter,
                                        const std::vector<std::vector<double> > &programmedTausSequence,
                                        const std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedSecondFundamentalFormSequence) {
    triangles[index].dialledProgTau =
            (1.0 - currDialInFactor) * programmedTausSequence[progTensorSequenceCounter][index] +
            currDialInFactor * programmedTausSequence[progTensorSequenceCounter + 1][index];
    triangles[index].programmedSecFF =
            (1.0 - rootCurrDialInFactor) *
            programmedSecondFundamentalFormSequence[progTensorSequenceCounter][index] +
            rootCurrDialInFactor * programmedSecondFundamentalFormSequence[progTensorSequenceCounter + 1][index];

}


void updateDialledSecondFundamentalForm(Triangle &triangle,
                                        double currDialInFactor,
                                        double rootCurrDialInFactor,
                                        const double &programmedTauCurrentStage,
                                        const double &programmedTauNextStage,
                                        const Eigen::Matrix<double, 2, 2> &programmedSecondFundamentalFormCurrentStage,
                                        const Eigen::Matrix<double, 2, 2> &programmedSecondFundamentalFormNextStage) {
    triangle.dialledProgTau =
            (1.0 - currDialInFactor) * programmedTauCurrentStage + currDialInFactor * programmedTauNextStage;
    triangle.programmedSecFF =
            (1.0 - rootCurrDialInFactor) * programmedSecondFundamentalFormCurrentStage +
            rootCurrDialInFactor * programmedSecondFundamentalFormNextStage;

}


void updateDeterminantOfDialledInverseProgrammedMetric(int index, std::vector<Triangle> &triangles) {
    triangles[index].programmedMetInvDet = triangles[index].programmedMetInv.determinant();
}


void updateDeterminantOfDialledInverseProgrammedMetric(Triangle &triangle) {
    triangle.programmedMetInvDet = triangle.programmedMetInv.determinant();
}


void updateGeometricPropertiesOfAllTriangles(std::vector<Triangle> &triangles, const std::vector<Node> &nodes) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].updateGeometricProperties(nodes);
    }
}


void
updateTriangleValues(const std::vector<Node> &nodes, std::vector<Triangle> &triangles, const SimulationStatus &status,
                     const double &dial_in_factor, const size_t &stage_counter, const SettingsNew &settings) {

    // Will use some functions of dial-in factor, pre-calculated here.
    const double dial_in_factor_root = sqrt(dial_in_factor);
    bool is_LCE_metric_used = settings.getCore().isLceModeEnabled() && !settings.getCore().isAnsatzMetricUsed();

//    updateMetrics(triangles, nodes);

#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
//        triangles[i].updateMetric(nodes);
        if (status == Dialling) {
            triangles[i].updateProgrammedQuantities(stage_counter, dial_in_factor, dial_in_factor_root, is_LCE_metric_used);
        }
//            if (is_LCE_metric_used) {
//                triangles[i].updateProgrammedMetricFromLCEInfo(stage_counter, dial_in_factor);
////                updateDialledInverseProgrammedMetricForATriangleWithoutAnsatz(i, triangles, dial_in_factor,
////                                                                              stage_counter,
////                                                                              programmedMetricInfoSequence);
//
//            } else {
//                triangles[i].updateProgrammedMetric(stage_counter, dial_in_factor);
////                updateDialledInverseProgrammedMetricFromAnsatz(triangles[i], dial_in_factor,
////                                                               programmedMetricSequence[stage_counter][i],
////                                                               programmedMetricSequence[stage_counter +
////                                                                                        1][i]);
//            }
//            updateDialledSecondFundamentalForm(triangles[i], dial_in_factor, dial_in_factor_root,
//                                               programmedTausSequence[stage_counter][i],
//                                               programmedTausSequence[stage_counter + 1][i],
//                                               programmedSecondFundamentalFormSequence[stage_counter][i],
//                                               programmedSecondFundamentalFormSequence[stage_counter +
//                                                                                       1][i]);
//
//            updateDeterminantOfDialledInverseProgrammedMetric(triangles[i]);
//        }
        triangles[i].updateGeometricProperties(nodes);
    }

//    updateGeometricPropertiesOfAllTriangles(triangles, nodes);
}