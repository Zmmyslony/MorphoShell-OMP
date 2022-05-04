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
#include <string> // For std::to_string
#include <omp.h>


#include "calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../SettingsStruct.hpp"
#include "../CustomOutStreamClass.hpp"


void calc_nonVertexPatchNodes_and_MatForPatchDerivs(
        const std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings,
        CustomOutStreamClass &logStream) {

    //Some temporary variables
    int numNonBoundaryTrisThatTriedMultiplePatchChoices = 0;

    omp_set_num_threads(8);
#pragma omp parallel for reduction (+ : numNonBoundaryTrisThatTriedMultiplePatchChoices)
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].refCentroid = (nodes[triangles[i].vertexLabels(0)].pos + nodes[triangles[i].vertexLabels(1)].pos +
                                    nodes[triangles[i].vertexLabels(2)].pos) / 3.0;

        std::vector<double> distancesToCentroid;
        std::vector<int> possiblePatchNodeLabels;
        for (int v = 0; v < 3; ++v) {
            for (int j = 0; j < nodes[triangles[i].vertexLabels(v)].incidentTriLabels.size(); ++j) {
                for (int w = 0; w < 3; ++w) {
                    int thisNodeLabel = triangles[nodes[triangles[i].vertexLabels(v)].incidentTriLabels(
                            j)].vertexLabels(w);
                    bool isNodeAlreadyAccountedFor = false;

                    for (int u = 0; u < 3; ++u) {
                        if (thisNodeLabel == triangles[i].vertexLabels(u)) {
                            isNodeAlreadyAccountedFor = true;
                        }
                    }


                    for (size_t u = 0; u < possiblePatchNodeLabels.size(); ++u) {
                        if (thisNodeLabel == possiblePatchNodeLabels[u]) {
                            isNodeAlreadyAccountedFor = true;
                        }
                    }

                    if (!isNodeAlreadyAccountedFor) {
                        distancesToCentroid.push_back((nodes[thisNodeLabel].pos - triangles[i].refCentroid).norm());
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
                    throw std::runtime_error("At least one search for patch nodes was exhausted without "
                                             "success (triangle " + std::to_string(triangles[i].label) +
                                             "); all possible patch matrices in the search had a  "
                                             "condition number above the acceptance threshold. Try increasing this threshold. If "
                                             "that does not solve the issue, or causes other issues, the patch node search will "
                                             "probably need to be extended. Please report this issue in that case. Aborting.");
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

//
//
//    double thisPatchSize;
//
//    std::vector<double> triangleLinearSizes(settings.NumTriangles);
//    std::vector<double> distancesToCentroid;
//    std::vector<int> possiblePatchNodeLabels;
//    int thisNodeLabel;
//    bool isNodeAlreadyAccountedFor;
//    std::vector<size_t> idxInDistanceList;
//
//    Eigen::Vector3d patchNodePos;
///*Matrix containing constants, patch node positions relative to centroid,
//their squares, and a cross term.*/
//    Eigen::Matrix<double, 6, 6> tempPatchNodeDataMatrix;
///* LU decomposition used to check invertibility of above matrix.*/
//    Eigen::FullPivLU<Eigen::Matrix<double, 6, 6> > tempPatchNodeDataMatrixDecomp;
//    Eigen::Matrix<double, 6, 6> invTempPatchNodeDataMatrix;
//    Eigen::Matrix<double, 6, 3> tempMatForPatchSecDerivs;
//
//
//    Eigen::JacobiSVD<Eigen::Matrix<double, 6, 3> > secDerivMatTempSVD;
//    double tempAbsConditionNumber;
//    double tempAbsConditionNumberDividedByThresholdValue;
//    for (int i = 0; i < settings.NumTriangles; ++i) {
//
//        /* Calculate the cartesian coordinates in the initial flat state of each
//        triangle's centroid, and store in a vector of 3-component vectors
//        (though for a flat starting state all the z-components will actually be
//        zero).*/
//        triangles[i].refCentroid = (nodes[triangles[i].vertexLabels(0)].pos + nodes[triangles[i].vertexLabels(1)].pos +
//                                    nodes[triangles[i].vertexLabels(2)].pos) / 3.0;
//
//        /* We now take a simple approach of assigning the
//        nonVertexPatchNodesLabels to nodes which are not vertices of the given
//        triangle, but are the closest nodes to the given triangle's initial
//        (flat state) centroid, (chosen from the surrounding triangles' vertices),
//        that give a reasonable answer for the matForPatchFirstDerivs. Some
//        choices of nodes lead to the inversion of a near-singular matrix en
//        route to matForPatchFirstDerivs, giving poor results, hence the
//        'reasonable' caveat. Note that if in future sheets are used which have
//        holes etc, a more complex procedure may be required/desired here.*/
//
//        // Loop over vertices of current triangles[i].
//        for (int v = 0; v < 3; ++v) {
//            // Loop over the triangles incident on each of these vertices.
//            for (int j = 0; j < nodes[triangles[i].vertexLabels(v)].incidentTriLabels.size(); ++j) {
//                // Loop over the vertices of these incident triangles.
//                for (int w = 0; w < 3; ++w) {
//                    thisNodeLabel = triangles[nodes[triangles[i].vertexLabels(v)].incidentTriLabels(j)].vertexLabels(w);
//                    isNodeAlreadyAccountedFor = false;
//
//                    /* Now check whether this node is a vertex of the
//                    triangle under consideration.*/
//                    for (int u = 0; u < 3; ++u) {
//                        if (thisNodeLabel == triangles[i].vertexLabels(u)) {
//                            isNodeAlreadyAccountedFor = true;
//                        }
//                    }
//
//                    /* Check whether current node has already been added to
//                    list of potential patch nodes.*/
//                    for (size_t u = 0; u < possiblePatchNodeLabels.size(); ++u) {
//                        if (thisNodeLabel == possiblePatchNodeLabels[u]) {
//                            isNodeAlreadyAccountedFor = true;
//                        }
//                    }
//
//                    /* If the current node is not a vertex, and is not already in
//                    the list of potentials, find its distance to the current
//                    triangle's initial centroid position, and store the
//                    distance to this node along with the others in a list of
//                    potentials. Store also this node's identity label similarly.
//                    */
//                    if (!isNodeAlreadyAccountedFor) {
//                        distancesToCentroid.push_back((nodes[thisNodeLabel].pos - triangles[i].refCentroid).norm());
//                        possiblePatchNodeLabels.push_back(thisNodeLabel);
//                    }
//                }
//            }
//        }
//
//
//        /* Now get a vector listing the *indices* (in the distancesToCentroid
//        vector) of distances to the current triangle centroid, sorted from
//        smallest to largest distance. This can then be used to extract the
//        corresponding possiblePatchNodeLabels in order of distance. */
//        idxInDistanceList.resize(distancesToCentroid.size());
//        for (size_t q = 0; q < idxInDistanceList.size(); ++q) {
//            idxInDistanceList[q] = q;
//        }
//        std::sort(std::begin(idxInDistanceList), std::end(idxInDistanceList),
//                  [&distancesToCentroid](const size_t &idx1, const size_t &idx2) -> bool {
//                      return distancesToCentroid[idx1] < distancesToCentroid[idx2];
//                  });
//        // idxInDistanceList is now sorted by distance to centroid.
//
//        // Put closest two nodes in nonVertexPatchNodesLabels.
//        for (int p = 0; p < 2; ++p) {
//            triangles[i].nonVertexPatchNodesLabels(p) = possiblePatchNodeLabels[idxInDistanceList[p]];
//        }
//
//        /* Now loop over possibilities for nonVertexPatchNodesLabels(3), finding
//        the first choice that gives reasonable results where the matrix to be
//        inverted is not close to singular. In principle the matrix could be
//        singular for all the possibilities in this loop. In that case this
//        search would need to be extended, but I think that case is extremely
//        unlikely (and perhaps impossible) for any reasonable mesh.*/
//
//        for (size_t q = 2; q < possiblePatchNodeLabels.size(); ++q) {
//
//            triangles[i].nonVertexPatchNodesLabels(2) = possiblePatchNodeLabels[idxInDistanceList[q]];
//
//
//            /* Calculate an approx `linear size' for the patch, taken to be the
//            RMS distances of the 6 patch nodes to the central triangle's
//            centroid. This measure is chosen to minimise the influence of any
//            one node, and, to be roughly unaffected by the geometric considerations
//            that lead to a badly conditioned problem. This leads to the chosen
//            threshold test having roughly the same level of `strictness' for each
//            triangles[i].*/
//            thisPatchSize = sqrt((
//                                         (nodes[triangles[i].vertexLabels(0)].pos -
//                                          triangles[i].refCentroid).squaredNorm() +
//                                         (nodes[triangles[i].vertexLabels(1)].pos -
//                                          triangles[i].refCentroid).squaredNorm() +
//                                         (nodes[triangles[i].vertexLabels(2)].pos -
//                                          triangles[i].refCentroid).squaredNorm() +
//                                         (nodes[triangles[i].nonVertexPatchNodesLabels(0)].pos -
//                                          triangles[i].refCentroid).squaredNorm() +
//                                         (nodes[triangles[i].nonVertexPatchNodesLabels(1)].pos -
//                                          triangles[i].refCentroid).squaredNorm() +
//                                         (nodes[triangles[i].nonVertexPatchNodesLabels(2)].pos -
//                                          triangles[i].refCentroid).squaredNorm()
//                                 ) / 6.0);
//
//            /* Calculate a matrix for each patch, that when multiplied
//            onto the current patch node coords, gives the coefficients specifying the
//            quadratic surface that goes through all 6 nodes in the patch. We will then
//            only retain and store the bottom three rows, which give the second
//            derivatives of position over the patch (wrt 2D parametrisation coords),
//            which is all that's needed for our 2nd F.F. approximation. If the full
//            quadratic surface was needed in future, the full matrix could easily be
//            retained and used instead.
//            If the 3D deformed state coordinates as a function of the 2D parametrisation
//            coordinates are X(x,y), Y(x,y), Z(x,y), we have the approximations:
//            X(x,y) = a_1 + b_1*(x-x_0) + c_1*(y-y_0) + d_1*(x-x_0)^2 + e_1*(y-y_0)^2 + f_1*(x-x_0)(y-y_0)
//            Y(x,y) = a_2 + b_2*(x-x_0) + ..... and equivalently for Z(x,y) etc,
//            where we approximate the surface as quadratic in the vicinity of a triangle
//            centroid at x_0, y_0. We do this for each triangle, and use the known
//            (x,y) and (X,Y,Z) coordinates of each of the patch nodes as constraints on
//            the coefficients a_1, b_1,....,f_3. Finding these coefficients then
//            simply requires inverting a 6x6 matrix, as below.*/
//
//            // Loop over patch nodes for this triangles[i].
//            for (int n = 0; n < 6; ++n) {
//
//                /* Patch nodes will always stick to this order: the three vertex
//                nodes first, then the others.*/
//                if (n < 3) {
//                    patchNodePos = nodes[triangles[i].vertexLabels(n)].pos;
//                } else {
//                    patchNodePos = nodes[triangles[i].nonVertexPatchNodesLabels(n - 3)].pos;
//                }
//
//                /* NB the following assumes that the initial state is a flat sheet,
//                the (x,y) coords on which then are the 2D parametrisation
//                coordinates, which then distort with the deforming sheet. If the
//                initial state were not flat, we would need to pick some less
//                obvious parametrisation to construct this matrix, and to define the
//                deformation gradient with respect to. The factors of 0.5 ensure
//                means that down the line, if we apply a part of the inverse of
//                this matrix, what we get are second derivatives of the patch,
//                rather than coefficients of the quadratic terms in the patch
//                expansion.*/
//                tempPatchNodeDataMatrix(0, n) = 1;
//                tempPatchNodeDataMatrix(1, n) = (patchNodePos(0) - triangles[i].refCentroid(0));
//                tempPatchNodeDataMatrix(2, n) = (patchNodePos(1) - triangles[i].refCentroid(1));
//                tempPatchNodeDataMatrix(3, n) = 0.5 * tempPatchNodeDataMatrix(1, n) * tempPatchNodeDataMatrix(1, n);
//                tempPatchNodeDataMatrix(4, n) = tempPatchNodeDataMatrix(1, n) * tempPatchNodeDataMatrix(2, n);
//                tempPatchNodeDataMatrix(5, n) = 0.5 * tempPatchNodeDataMatrix(2, n) * tempPatchNodeDataMatrix(2, n);
//            }
//
//            /* For this triangle:
//            (a_1, b_1, c_1, d_1, e_1, f_1) * tempPatchNodeDataMatrix= (X1, X2, X3, X4, X5, X6)
//            where X1, X2,..,X6 will be the X coordinates of the 6 patch nodes in any
//            3D deformed state. Similarly for the 'Y's and 'Z's with (a_2, b_2,...)
//            and (a_3, b_3,...) respectively. Thus inverting tempPatchNodeDataMatrix
//            allows us to find the 'a's, 'b's,..., 'f's from any deformed state
//            3D coords of the 6 patch nodes. We check first that
//            tempPatchNodeDataMatrix is actually invertible.
//            */
//            tempPatchNodeDataMatrixDecomp.compute(tempPatchNodeDataMatrix);
//            bool isMatReversible = true;
//            if (!tempPatchNodeDataMatrixDecomp.isInvertible()) {
//                isMatReversible = false;
//            } else {
//                invTempPatchNodeDataMatrix = tempPatchNodeDataMatrix.inverse();
//
//                /* Only one of these will end up being used - which depends on
//                which secFF estimate is being used.*/
//                tempMatForPatchSecDerivs = invTempPatchNodeDataMatrix.block<6, 3>(0, 3);
//
//
//                /* Compute absolute condition number i.e. the maximum singular
//                value of the matrix which will be used to find the patch
//                derivatives, to check that it will give reasonable results. If a
//                certain direction is not well 'sampled' by the patch for example,
//                the condition number will be high, and the results of using it
//                would be poor (and can lead to overly strict time step
//                restrictions), so an alternative choice of patch nodes should be
//                considered instead.*/
//                secDerivMatTempSVD.compute(tempMatForPatchSecDerivs);
//                tempAbsConditionNumber = secDerivMatTempSVD.singularValues()(0); //This has dimensions 1/Length^2.
//                tempAbsConditionNumberDividedByThresholdValue = tempAbsConditionNumber * thisPatchSize * thisPatchSize /
//                                                                settings.PatchMatrixDimensionlessConditioningThreshold;
//            }
//
//            /* Only proceed if this matrix will give reasonable results, i.e. has
//            a low enough condition number. The threshold can be tuned in the
//            settings file. It will be mesh dependent, with poorer meshes
//            giving higher condition numbers, and thus requiring a higher
//            threshold to actually find any patch that is considered good enough.
//            The search may be extended in future versions of the code to help
//            with this, if it becomes necessary.*/
//            if (tempAbsConditionNumberDividedByThresholdValue >= 1.0 || !isMatReversible) {
//
//                /* The first time the above condition occurs for a triangle,
//                add one to a counter if the triangle in question is *not* on the
//                boundary. This is because if lots of interior triangles are
//                going through multiple patch possibilities before satisfying the
//                determinant threshold conditions,
//                settings.PatchMatrixDimensionlessConditioningThreshold is likely
//                too high, or something may be strange about the mesh. */
//                if (q == 2 && !triangles[i].isOnBoundary) {
//                    numNonBoundaryTrisThatTriedMultiplePatchChoices += 1;
//                }
//
//                //Throw error to main if whole search has been exhausted unsuccessfully
//                if (q == possiblePatchNodeLabels.size() - 1) {
//                    throw std::runtime_error("At least one search for patch nodes was exhausted without "
//                                             "success (triangle " + std::to_string(i) +
//                                             "); all possible patch matrices in the search had a  "
//                                             "condition number above the acceptance threshold. Try increasing this threshold. If "
//                                             "that does not solve the issue, or causes other issues, the patch node search will "
//                                             "probably need to be extended. Please report this issue in that case. Aborting.");
//                }
//
//                //Move onto the next possible non-vertex patch node since this
//                //one has produced an unsuitable matrix.
//                continue;
//            } else {
//                /* Store the matrix that gives the first or second derivatives
//                of the patch, depending on which secFF approx is being used.*/
//
//                triangles[i].matForPatchSecDerivs = tempMatForPatchSecDerivs;
//
//                //Move onto next triangle as search has now been successfully
//                //completed for this one.
//                break;
//            }
//        }
//
//        //Reset containers before moving onto next triangle
//        distancesToCentroid.resize(0);
//        possiblePatchNodeLabels.resize(0);
//    }
//}