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

This file defines the member functions for the Triangle class that holds
the data for each triangular element. The class constructors are the exception;
they are left in the header file for clarity there.*/

#include <iostream>
#include <iomanip>
#include <set>
#include <unordered_set>
#include <Eigen/Dense>

#include "Triangle.hpp"


//This is a debugging tool to display the node's data
void Triangle::display() {
    std::cout << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    std::cout << "Triangle " << label << ":" << std::endl;
    std::cout << "Boundary indicator = " << isOnBoundary << std::endl;
    std::cout << "Initial (reference) area = " << initArea << std::endl;
    std::cout << "invCurrArea = " << currAreaInv << std::endl;
    std::cout << "Labels of vertices: " << vertexLabels.transpose() << std::endl;
    std::cout << "Labels of edges: " << edgeLabels.transpose() << std::endl;
    std::cout << "Initial non-boundary edge length fractions: " << initNonBoundEdgeLengthFracs.transpose()
              << std::endl;
    std::cout << "Labels of adjacent (edge-sharing) triangles: " << edgeSharingTriLabels.transpose() << std::endl;
    std::cout << "edgeAdjTriLabelSelectors: " << edgeAdjTriLabelSelectors.transpose() << std::endl;
    std::cout << "indicesIntoEdgeSharingTriLabelsOfNeighbours: "
              << indicesIntoEdgeSharingTriLabelsOfNeighbours.transpose() << std::endl;
    std::cout << "Labels of non-vertex nodes in this triangle's patch: " << nonVertexPatchNodesLabels.transpose()
              << std::endl;
    std::cout << "Current sides = " << "\n" << currSides << std::endl;
    std::cout << "faceNormal = " << "\n" << faceNormal << std::endl;
    std::cout << "Initial outward side normals = " << "\n" << initOutwardSideNormals << std::endl;
    std::cout << "invInitInPlaneSidesMat = " << "\n" << invInitSidesMat << std::endl;
    std::cout << "Dialled in inverse of programmed metric = " << "\n" << programmedMetInv << std::endl;
    std::cout << "detDialledInvProgMetric = " << "\n" << programmedMetInvDet << std::endl;
    std::cout << "Dialled in prog tau factor = " << dialledProgTau << std::endl;
    std::cout << "Dialled in programmed second fundamental form = " << "\n" << programmedSecFF << std::endl;
    std::cout << "Deformation gradient = " << "\n" << defGradient << std::endl;
    std::cout << "metric = " << "\n" << met << std::endl;
    std::cout << "Inverse of Metric = " << "\n" << metInv << std::endl;
    std::cout << "Det of inverse of Metric = " << "\n" << metInvDet << std::endl;
    std::cout << "matForPatchSecDerivs = " << "\n" << matForPatchSecDerivs << std::endl;
    std::cout << "patchSecDerivs = " << "\n" << patchSecDerivs << std::endl;
    std::cout << "Second fundamental form = " << "\n" << secFF << std::endl;
    std::cout << "Bending energy density deriv wrt secFF = " << "\n" << energyDensityDerivWRTSecFF << std::endl;
    std::cout << "bendEnergyDensityDerivWRTMetric = " << "\n" << bendEnergyDensityDerivWRTMetric << std::endl;
    std::cout << "halfPK1Stress = " << "\n" << halfPK1Stress << std::endl;
    std::cout << "-----------------------------" << std::endl;
}

void Triangle::updateHalfPK1Stress(double stretchingPrefactor) {
    halfPK1Stress = defGradient * (stretchingPrefactor * dialledProgTau *
                                   (programmedMetInv - (metInvDet / programmedMetInvDet) * metInv));
}

Eigen::Matrix<double, 3, 3> Triangle::getStretchingForces() {
    return halfPK1Stress * initOutwardSideNormals;
}

// Returns vectors that normal to triangle edges and point outward
Eigen::Matrix<double, 3, 3> Triangle::getTriangleEdgeNormals() {
    Eigen::Matrix<double, 3, 3> triangleEdgeNormals;

    triangleEdgeNormals.col(1) = faceNormal.cross(currSides.col(1));
    triangleEdgeNormals.col(2) = -faceNormal.cross(currSides.col(0));
    triangleEdgeNormals.col(0) = -triangleEdgeNormals.col(1) - triangleEdgeNormals.col(2);
    return triangleEdgeNormals;
}

Eigen::Matrix<double, 3, 1> Triangle::getBendingForce(const Eigen::Matrix<double, 3, 3> &normalDerivatives, int row) {
    Eigen::Matrix<double, 2, 2> secFFDerivPreFacMat;

    secFFDerivPreFacMat(0, 0) = matForPatchSecDerivs(row, 0);
    secFFDerivPreFacMat(0, 1) = matForPatchSecDerivs(row, 1);
    secFFDerivPreFacMat(1, 1) = matForPatchSecDerivs(row, 2);

    if (row < 3) {
        secFFDerivPreFacMat(0, 0) += normalDerivatives(0, row);
        secFFDerivPreFacMat(0, 1) += normalDerivatives(1, row);
        secFFDerivPreFacMat(1, 1) += normalDerivatives(2, row);
    }

    secFFDerivPreFacMat(1, 0) = secFFDerivPreFacMat(0, 1);

    Eigen::Vector3d bendingForce;
    bendingForce = -initArea * (energyDensityDerivWRTSecFF * secFFDerivPreFacMat).trace() * faceNormal;
    return bendingForce;
}

void Triangle::updateMetric(const std::vector<Node> &nodes) {
    currSides.col(0).noalias() = nodes[vertexLabels(1)].pos - nodes[vertexLabels(0)].pos;
    currSides.col(1).noalias() = nodes[vertexLabels(2)].pos - nodes[vertexLabels(0)].pos;
    defGradient.noalias() = currSides * invInitSidesMat;

    met.noalias() = defGradient.transpose() * defGradient;
    metInv.noalias() = met.inverse();
    metInvDet = metInv.determinant();
}

void Triangle::updateGeometricProperties(const std::vector<Node> &nodes) {
    updateMetric(nodes);
    Eigen::Matrix<double, 3, 6> matrixOfPatchNodeCoords;
    for (int n = 0; n < 6; ++n) {
        if (n < 3) {
            matrixOfPatchNodeCoords.col(n).noalias() = nodes[vertexLabels(n)].pos;
        } else {
            matrixOfPatchNodeCoords.col(n).noalias() = nodes[nonVertexPatchNodesLabels(n - 3)].pos;
        }
    }
    patchSecDerivs.noalias() = matrixOfPatchNodeCoords * matForPatchSecDerivs;

    faceNormal.noalias() = currSides.col(0).cross(currSides.col(1));
    currAreaInv = 2 / faceNormal.norm();
    faceNormal.noalias() = 0.5 * faceNormal * currAreaInv; // Normalising

    centroid = (nodes[vertexLabels(0)].pos + nodes[vertexLabels(1)].pos + nodes[vertexLabels(2)].pos) / 3;
}


void Triangle::updateSecondFundamentalForm(double bendingPreFac, double JPreFactor, double poissonRatio) {
    Eigen::Vector3d vectorOfSecFFComps = patchSecDerivs.transpose() * faceNormal;

    secFF(0, 0) = vectorOfSecFFComps(0);
    secFF(0, 1) = vectorOfSecFFComps(1);
    secFF(1, 0) = vectorOfSecFFComps(1);
    secFF(1, 1) = vectorOfSecFFComps(2);

    Eigen::Matrix<double, 2, 2> scaledSecFF;

    Eigen::Matrix<double, 2, 2> relativeSecFF = programmedMetInv * (secFF - programmedSecFF);
    double areaMultiplier = bendingPreFac * programmedMetInvDet;

    double J = areaMultiplier * JPreFactor;

    // Now calculate bending energy density for this 
    double preGentBendEnergyDensity = areaMultiplier * ((1 - poissonRatio) * (relativeSecFF * relativeSecFF).trace() +
                                                        poissonRatio * relativeSecFF.trace() * relativeSecFF.trace());
    double gentDerivFac = (1 + 2 * preGentBendEnergyDensity / J);

    /* Calculate the derivative of the bending energy density with respect to the secFF. */
    energyDensityDerivWRTSecFF = gentDerivFac * 2 * areaMultiplier *
                                 ((1 - poissonRatio) * relativeSecFF * programmedMetInv +
                                  poissonRatio * relativeSecFF.trace() * programmedMetInv);
    bendEnergyDensity = gentDerivFac * preGentBendEnergyDensity;
}

void Triangle::updateFirstFundamentalForm(double stretchingPreFac) {
    stretchEnergyDensity = stretchingPreFac * dialledProgTau *
                           ((met * programmedMetInv).trace() +
                            metInvDet / programmedMetInvDet -
                            3 / dialledProgTau);
}


void Triangle::updateAngleDeficits(std::vector<double> &angleDeficits) const {
    std::vector<Eigen::Vector3d> sides;
    sides.emplace_back(currSides.col(0));
    sides.emplace_back(currSides.col(1));
    sides.emplace_back(sides[1] - sides[0]);

    std::vector<double> sidesLength(sides.size());
    for (auto &side: sides) {
        sidesLength.emplace_back(side.norm());
    }

    angleDeficits[vertexLabels(0)] -= acos(sides[0].dot(sides[1]) / (sidesLength[0] * sidesLength[1]));
    angleDeficits[vertexLabels(1)] -= acos(sides[0].dot(sides[2]) / (sidesLength[0] * sidesLength[2]));
    angleDeficits[vertexLabels(0)] -= acos(sides[1].dot(sides[2]) / (sidesLength[1] * sidesLength[2]));
}

double Triangle::getLinearSize() const {
    double first_side_length = currSides.col(0).norm();
    double second_side_length = currSides.col(1).norm();
    double third_side_length = (currSides.col(0) - currSides.col(1)).norm();

    double longest_side = std::max({first_side_length, second_side_length, third_side_length});

    double shortest_altitude = 2 * initArea / longest_side;
    return shortest_altitude;
}

template<typename T>
T interpolate(const T &previous, const T &next, double ratio) {
    return (1 - ratio) * previous + ratio * next;
}

void Triangle::updateProgrammedMetricImplicit(double dirAngle, double lambda, double nu) {
    double cosDirAng = cos(dirAngle);
    double sinDirAng = sin(dirAngle);
    double lambdaToTheMinus2 = pow(lambda, -2);
    double lambdaToThe2Nu = pow(lambda, 2 * nu);

    programmedMetInv(0, 0) = lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
    programmedMetInv(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
    programmedMetInv(1, 0) = programmedMetInv(0, 1);
    programmedMetInv(1, 1) = lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
    programmedMetInvDet = programmedMetInv.determinant();
}

void Triangle::updateProgrammedMetricImplicit(int stage_counter, double dial_in_factor,
                                              bool is_elongation_dynamically_updated, double transfer_coefficient) {
    Eigen::Vector3d metric_current = interpolate(programmed_metric_infos[stage_counter],
                                                 programmed_metric_infos[stage_counter + 1], dial_in_factor);

    double dirAng = metric_current(0);
    double lambda = metric_current(1);
    double nu = metric_current(2);
    if (is_elongation_dynamically_updated) {
        double target_elongation = 1 - (1 - lambda) * relative_height;
        local_elongation = interpolate(local_elongation, target_elongation, transfer_coefficient);
        lambda = local_elongation;
    }
    updateProgrammedMetricImplicit(dirAng, lambda, nu);
}

void Triangle::updateProgrammedMetricExplicit(int stage_counter, double dial_in_factor) {
    programmedMetInv = interpolate(programmed_metric_inv[stage_counter], programmed_metric_inv[stage_counter + 1],
                                   dial_in_factor);
    programmedMetInvDet = programmedMetInv.determinant();
}

void Triangle::updateProgrammedSecondFundamentalForm(int stage_counter, double dial_in_factor_root) {
    programmedSecFF = interpolate(programmed_second_fundamental_form[stage_counter],
                                  programmed_second_fundamental_form[stage_counter + 1], dial_in_factor_root);
}

void Triangle::updateProgrammedTaus(int stage_counter, double dial_in_factor) {
    dialledProgTau = interpolate(programmed_taus[stage_counter], programmed_taus[stage_counter + 1], dial_in_factor);
}

void Triangle::updateProgrammedQuantities(int stage_counter, double dial_in_factor, double dial_in_factor_root,
                                          bool is_lce_metric_used, bool is_elongation_dynamically_updated,
                                          double transfer_coefficient) {
    if (is_lce_metric_used) {
        updateProgrammedMetricImplicit(stage_counter, dial_in_factor, is_elongation_dynamically_updated,
                                       transfer_coefficient);
    } else {
        updateProgrammedMetricExplicit(stage_counter, dial_in_factor);
    }
    updateProgrammedSecondFundamentalForm(stage_counter, dial_in_factor_root);
    updateProgrammedTaus(stage_counter, dial_in_factor);
}


std::vector<unsigned int> getNeighbouringNodes(const Triangle &current_triangle, const std::vector<Triangle> &triangles,
                                               const std::vector<Node> &nodes) {
    // Provides nodes, which contain edges that neighbour the nodes of the current triangle, exclusive of those of
    // considered triangle
    std::unordered_set<unsigned int> neighbouring_triangle_indices;
    std::unordered_set<unsigned int> neighbouring_node_indices;

    for (unsigned int i = 0; i < 3; i++) {
        unsigned int i_node = current_triangle.vertexLabels(i);
        const Node &node = nodes[i_node];
        for (unsigned int j_triangles: node.incidentTriLabels) {
            neighbouring_triangle_indices.insert(j_triangles);
        }
    }

    for (unsigned int j_triangles: neighbouring_triangle_indices) {
        const Triangle &triangle = triangles[j_triangles];
        for (unsigned int j: triangle.vertexLabels) {
            neighbouring_node_indices.insert(j);
        }
    }

    // We want to exclude the labels of the original nodes
    for (unsigned int i: current_triangle.vertexLabels) {
        neighbouring_node_indices.erase(i);
    }

    std::vector<unsigned int> v_node_indices;
    for (unsigned int index: neighbouring_node_indices) {
        v_node_indices.emplace_back(index);
    }
    return v_node_indices;
}

std::vector<std::pair<unsigned int, double>>
assignDistanceFromCentroid(const std::vector<unsigned int> &node_indices, const std::vector<Node> &nodes,
                           const Eigen::Vector3d &position) {
    std::vector<std::pair<unsigned int, double>> index_distance_pair;
    for (unsigned int index: node_indices) {
        double distance = (nodes[index].pos - position).norm();
        index_distance_pair.emplace_back(index, distance);
    }
    return index_distance_pair;
}

Eigen::Matrix<double, 6, 1> patchColumn(const Eigen::Vector3d &position, const Eigen::Vector3d &centroid) {
    Eigen::Matrix<double, 6, 1> patchColumn;
    patchColumn(0, 1) = 1;
    patchColumn(1, 1) = (position(0) - centroid(0));
    patchColumn(2, 1) = (position(1) - centroid(1));
    patchColumn(3, 1) = 0.5 * patchColumn(1, 1) * patchColumn(1, 1);
    patchColumn(4, 1) = patchColumn(1, 1) * patchColumn(2, 1);
    patchColumn(5, 1) = 0.5 * patchColumn(2, 1) * patchColumn(2, 1);
    return patchColumn;
}

double Triangle::getHeight() const {
    return centroid(2);
}

int Triangle::updateMatForPatchDerivs(const std::vector<Triangle> &triangles, const std::vector<Node> &nodes,
                                      double patch_threshold) {
    refCentroid = (nodes[vertexLabels(0)].pos +
                   nodes[vertexLabels(1)].pos +
                   nodes[vertexLabels(2)].pos) / 3;

    std::vector<unsigned int> possiblePatchNodeLabels = getNeighbouringNodes(*this, triangles, nodes);
    std::vector<std::pair<unsigned int, double>> indexDistancePairs = assignDistanceFromCentroid(
            possiblePatchNodeLabels, nodes, refCentroid);

    std::sort(std::begin(indexDistancePairs), std::end(indexDistancePairs),
              [](const std::pair<unsigned int, double> &p1,
                 const std::pair<unsigned int, double> &p2) -> bool {
                  return p1.second < p2.second;
              });

    int invalid_non_boundary_triangle_count = 0;
    Eigen::Matrix<double, 6, 6> patchNodeDataMatrix;
    double raw_path_size_factor = (nodes[vertexLabels(0)].pos - refCentroid).squaredNorm() +
                                  (nodes[vertexLabels(1)].pos - refCentroid).squaredNorm() +
                                  (nodes[vertexLabels(2)].pos - refCentroid).squaredNorm();

    for (int n = 0; n < 3; ++n) {
        Eigen::Vector3d candidatePatchNode;
        candidatePatchNode = nodes[vertexLabels(n)].pos;
        patchNodeDataMatrix.col(n) = patchColumn(candidatePatchNode, refCentroid);
    }

    std::set<std::vector<unsigned int>> candidate_trios;
    for (unsigned int p: possiblePatchNodeLabels) {
        for (unsigned int q: possiblePatchNodeLabels) {
            for (unsigned int r: possiblePatchNodeLabels) {
                candidate_trios.insert({p, q, r});
            }
        }
    }

    double lowestConditionNumber = patch_threshold;
    for (auto &candidateIndices: candidate_trios) {
        patchNodeDataMatrix.col(3) = patchColumn(nodes[candidateIndices[0]].pos, refCentroid);
        patchNodeDataMatrix.col(4) = patchColumn(nodes[candidateIndices[1]].pos, refCentroid);
        patchNodeDataMatrix.col(5) = patchColumn(nodes[candidateIndices[2]].pos, refCentroid);
        double patch_size = sqrt(
                (raw_path_size_factor +
                 (nodes[candidateIndices[0]].pos - refCentroid).squaredNorm() +
                 (nodes[candidateIndices[1]].pos - refCentroid).squaredNorm() +
                 (nodes[candidateIndices[2]].pos - refCentroid).squaredNorm()) / 6);

        Eigen::FullPivLU<Eigen::Matrix<double, 6, 6>> tempPatchNodeDataMatrixDecomp;
        tempPatchNodeDataMatrixDecomp.compute(patchNodeDataMatrix);
        bool isMatReversible = tempPatchNodeDataMatrixDecomp.isInvertible();

        if (isMatReversible) {
            Eigen::Matrix<double, 6, 6> invTempPatchNodeDataMatrix = patchNodeDataMatrix.inverse();
            Eigen::Matrix<double, 6, 3> candidatePatchDiv;
            candidatePatchDiv = invTempPatchNodeDataMatrix.block<6, 3>(0, 3);

            Eigen::JacobiSVD<Eigen::Matrix<double, 6, 3>> secDerivMatTempSVD;
            secDerivMatTempSVD.compute(candidatePatchDiv);

            double singular_values = secDerivMatTempSVD.singularValues()(0); // This has dimensions 1 / Length ^ 2.
            double conditionNumber = singular_values * pow(patch_size, 2);

            if (conditionNumber < lowestConditionNumber) {
                lowestConditionNumber = conditionNumber;
                nonVertexPatchNodesLabels = {candidateIndices[0], candidateIndices[1], candidateIndices[2]};
                matForPatchSecDerivs = candidatePatchDiv;
            }
        }
    }

    if (lowestConditionNumber < patch_threshold) {
        return invalid_non_boundary_triangle_count;
    }

    throw std::runtime_error(
            "T" + std::to_string(label) +
            ": Search for patch nodes was exhausted without success; "
            "all possible patch matrices in the search had a condition number above the acceptance "
            "threshold. Try increasing this threshold. If that does not solve the issue, or causes "
            "other issues, the patch node search will probably need to be extended. Please report this "
            "issue in that case. Aborting.");
}


void Triangle::setRelativeHeight(double relative_height) {
    Triangle::relative_height = relative_height;
}

void Triangle::setLocalElongation(double local_elongation) {
    Triangle::local_elongation = local_elongation;
}

Triangle::Triangle(int label, int id_0, int id_1, int id_2) : Triangle() {
    Triangle::label = label;
    vertexLabels(0) = id_0;
    vertexLabels(1) = id_1;
    vertexLabels(2) = id_2;
}
