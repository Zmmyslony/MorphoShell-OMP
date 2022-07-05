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


void calc_nonVertexPatchNodes_and_MatForPatchDerivs(
        const std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings,
        CustomOutStreamClass &logStream) {

    //Some temporary variables
    int numNonBoundaryTrisThatTriedMultiplePatchChoices = 0;

#pragma omp parallel for reduction (+ : numNonBoundaryTrisThatTriedMultiplePatchChoices)
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].refCentroid = (nodes[triangles[i].vertexLabels(0)].pos + nodes[triangles[i].vertexLabels(1)].pos +
                                    nodes[triangles[i].vertexLabels(2)].pos) / 3.0;

        std::vector<double> distancesToCentroid;
        std::vector<int> possiblePatchNodeLabels;
        for (int v = 0; v < 3; ++v) {
            Eigen::VectorXi incidentTriangleLabels = nodes[triangles[i].vertexLabels(v)].incidentTriLabels;
            for (int j = 0; j < incidentTriangleLabels.size(); ++j) {
                for (int w = 0; w < 3; ++w) {
                    int thisNodeLabel = triangles[incidentTriangleLabels(j)].vertexLabels(w);
                    bool isNodeAlreadyAccountedFor = false;

                    for (int u = 0; u < 3; ++u) {
                        if (thisNodeLabel == triangles[i].vertexLabels(u)) {
                            isNodeAlreadyAccountedFor = true;
                        }
                    }


                    for (int possiblePatchNodeLabel: possiblePatchNodeLabels) {
                        if (thisNodeLabel == possiblePatchNodeLabel) {
                            isNodeAlreadyAccountedFor = true;
                        }
                    }

                    if (!isNodeAlreadyAccountedFor) {
                        double distanceToCentroid = (nodes[thisNodeLabel].pos - triangles[i].refCentroid).norm();
                        distancesToCentroid.push_back(distanceToCentroid);
                        possiblePatchNodeLabels.push_back(thisNodeLabel);
                    }
                }
            }
        }

        std::vector<size_t> idxInDistanceList(distancesToCentroid.size());

        for (size_t q = 0; q < idxInDistanceList.size(); ++q) {
            idxInDistanceList[q] = q;
        }
        std::sort(std::begin(idxInDistanceList), std::end(idxInDistanceList),
                  [&distancesToCentroid](const size_t &idx1, const size_t &idx2) -> bool {
                      return distancesToCentroid[idx1] < distancesToCentroid[idx2];
                  });

        for (int p = 0; p < 2; ++p) {
            triangles[i].nonVertexPatchNodesLabels(p) = possiblePatchNodeLabels[idxInDistanceList[p]];
        }


        for (size_t q = 2; q < possiblePatchNodeLabels.size(); ++q) {
            triangles[i].nonVertexPatchNodesLabels(2) = possiblePatchNodeLabels[idxInDistanceList[q]];

            double thisPatchSize = sqrt((
                                                (nodes[triangles[i].vertexLabels(0)].pos -
                                                 triangles[i].refCentroid).squaredNorm() +
                                                (nodes[triangles[i].vertexLabels(1)].pos -
                                                 triangles[i].refCentroid).squaredNorm() +
                                                (nodes[triangles[i].vertexLabels(2)].pos -
                                                 triangles[i].refCentroid).squaredNorm() +
                                                (nodes[triangles[i].nonVertexPatchNodesLabels(0)].pos -
                                                 triangles[i].refCentroid).squaredNorm() +
                                                (nodes[triangles[i].nonVertexPatchNodesLabels(1)].pos -
                                                 triangles[i].refCentroid).squaredNorm() +
                                                (nodes[triangles[i].nonVertexPatchNodesLabels(2)].pos -
                                                 triangles[i].refCentroid).squaredNorm()
                                        ) / 6.0);

            Eigen::Vector3d patchNodePos;
            Eigen::Matrix<double, 6, 6> tempPatchNodeDataMatrix;
            for (int n = 0; n < 6; ++n) {
                if (n < 3) {
                    patchNodePos = nodes[triangles[i].vertexLabels(n)].pos;
                } else {
                    patchNodePos = nodes[triangles[i].nonVertexPatchNodesLabels(n - 3)].pos;
                }

                tempPatchNodeDataMatrix(0, n) = 1;
                tempPatchNodeDataMatrix(1, n) = (patchNodePos(0) - triangles[i].refCentroid(0));
                tempPatchNodeDataMatrix(2, n) = (patchNodePos(1) - triangles[i].refCentroid(1));
                tempPatchNodeDataMatrix(3, n) = 0.5 * tempPatchNodeDataMatrix(1, n) * tempPatchNodeDataMatrix(1, n);
                tempPatchNodeDataMatrix(4, n) = tempPatchNodeDataMatrix(1, n) * tempPatchNodeDataMatrix(2, n);
                tempPatchNodeDataMatrix(5, n) = 0.5 * tempPatchNodeDataMatrix(2, n) * tempPatchNodeDataMatrix(2, n);
            }

            Eigen::FullPivLU<Eigen::Matrix<double, 6, 6>> tempPatchNodeDataMatrixDecomp;
            tempPatchNodeDataMatrixDecomp.compute(tempPatchNodeDataMatrix);
            bool isMatReversible = true;

            Eigen::Matrix<double, 6, 3> tempMatForPatchSecDerivs;
            double tempAbsConditionNumberDividedByThresholdValue = 0;

            if (!tempPatchNodeDataMatrixDecomp.isInvertible()) {
                isMatReversible = false;
            } else {
                Eigen::Matrix<double, 6, 6> invTempPatchNodeDataMatrix = tempPatchNodeDataMatrix.inverse();
                tempMatForPatchSecDerivs = invTempPatchNodeDataMatrix.block<6, 3>(0, 3);

                Eigen::JacobiSVD<Eigen::Matrix<double, 6, 3>> secDerivMatTempSVD;
                secDerivMatTempSVD.compute(tempMatForPatchSecDerivs);

                double tempAbsConditionNumber = secDerivMatTempSVD.singularValues()(
                        0); //This has dimensions 1/Length^2.
                tempAbsConditionNumberDividedByThresholdValue = tempAbsConditionNumber * thisPatchSize * thisPatchSize /
                                                                settings.PatchMatrixDimensionlessConditioningThreshold;
            }

            if (tempAbsConditionNumberDividedByThresholdValue >= 1.0 || !isMatReversible) {
                if (q == 2 && !triangles[i].isOnBoundary) {
                    numNonBoundaryTrisThatTriedMultiplePatchChoices += 1;
                }

                //Throw error to main if whole search has been exhausted unsuccessfully
                if (q == possiblePatchNodeLabels.size() - 1) {
                    throw std::runtime_error(
                            "At least one search for patch nodes was exhausted without success (triangle " +
                            std::to_string(triangles[i].label) +
                            "); all possible patch matrices in the search had a condition number above the acceptance "
                            "threshold. Try increasing this threshold. If that does not solve the issue, or causes "
                            "other issues, the patch node search will probably need to be extended. Please report this "
                            "issue in that case. Aborting.");
                }
                continue;
            } else {
                triangles[i].matForPatchSecDerivs = tempMatForPatchSecDerivs;
                break;
            }
        }
    }

    /*Print also the number of non-boundary triangles that had to search through
    multiple possible patch options to find one satisfying the determinant
    condition. It seems unlikely that the determinant will be a small for
    triangles in the interior of a reasonable mesh, so if this number is large,
    that suggests something suspicious. One explanation might be that
    settings.PatchMatrixDimensionlessConditioningThreshold has been set to too
    low a value*/
    logStream.open();
    logStream << "Number of non-boundary triangles that had to search through \n" <<
              "multiple possible patch options to find one \nsatisfying the condition number " <<
              "criterion was " << numNonBoundaryTrisThatTriedMultiplePatchChoices <<
              ", \nwhich should not be a large proportion of the mesh's triangles." << std::endl;
    logStream.close();
}