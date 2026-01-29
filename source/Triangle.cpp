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
std::stringstream Triangle::display() {
    std::stringstream msg;
    msg << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    msg << "Triangle " << label << ":" << std::endl;
    msg << "Boundary indicator = " << isOnBoundary << std::endl;
    msg << "Initial (reference) area = " << initArea << std::endl;
    msg << "invCurrArea = " << currAreaInv << std::endl;
    msg << "Labels of vertices: " << vertexLabels.transpose() << std::endl;
    msg << "Labels of edges: " << edgeLabels.transpose() << std::endl;
    msg << "Initial non-boundary edge length fractions: " << initNonBoundEdgeLengthFracs.transpose()
        << std::endl;

    msg << "Current sides = " << "\n" << currSides << std::endl;
    msg << "Initial outward side normals = " << "\n" << initOutwardSideNormals << std::endl;
    msg << "invInitInPlaneSidesMat = " << "\n" << invInitSidesMat << std::endl;
    msg << "Dialled in inverse of programmed metric = " << "\n" << programmedMetInv << std::endl;
    msg << "detDialledInvProgMetric = " << "\n" << programmedMetInvDet << std::endl;
    msg << "Dialled in prog tau factor = " << dialledProgTau << std::endl;
    msg << "Dialled in programmed second fundamental form = " << "\n" << programmedSecFF << std::endl;
    msg << "Deformation gradient = " << "\n" << defGradient << std::endl;
    msg << "metric = " << "\n" << met << std::endl;
    msg << "Inverse of Metric = " << "\n" << metInv << std::endl;
    msg << "Det of inverse of Metric = " << "\n" << metInvDet << std::endl;
    msg << "matForPatchSecDerivs = " << "\n" << matForPatchSecDerivs << std::endl;
    msg << "Second fundamental form = " << "\n" << secFF << std::endl;
    msg << "halfPK1Stress = " << "\n" << getHalfPK1Stress(1) << std::endl;
    msg << "-----------------------------" << std::endl;
    return msg;
}

Eigen::Matrix<double, 3, 2> Triangle::getHalfPK1Stress(double stretchingPrefactor) const {
    return defGradient * (stretchingPrefactor * dialledProgTau *
        (programmedMetInv - (metInvDet / programmedMetInvDet) * metInv));
}

Eigen::Matrix<double, 3, 3> Triangle::getStretchingForces(double stretchingPrefactor) const {
    return getHalfPK1Stress(stretchingPrefactor) * initOutwardSideNormals;
}

// Returns vectors that normal to triangle edges and point outward
Eigen::Matrix<double, 3, 3> Triangle::getTriangleEdgeNormals(const Eigen::Matrix<double, 3, 2>& current_sides, const Eigen::Vector3d& face_normal) const {
    Eigen::Matrix<double, 3, 3> triangleEdgeNormals;

    triangleEdgeNormals.col(1).noalias() = face_normal.cross(current_sides.col(1));
    triangleEdgeNormals.col(2).noalias() = -face_normal.cross(current_sides.col(0));
    triangleEdgeNormals.col(0).noalias() = -triangleEdgeNormals.col(1) - triangleEdgeNormals.col(2);
    return triangleEdgeNormals;
}

Eigen::Matrix<double, 3, 1> Triangle::getBendingForceNode(const Eigen::Vector3d& normalDerivatives, int row, const Eigen::Vector3d& faceNormal, const Eigen::Matrix<double, 2, 2>
                                                          & energyDensityDerivWRTSecFF) const {
    Eigen::Matrix<double, 2, 2> secFFDerivPreFacMat;

    secFFDerivPreFacMat(0, 0) = matForPatchSecDerivs(row, 0);
    secFFDerivPreFacMat(0, 1) = matForPatchSecDerivs(row, 1);
    secFFDerivPreFacMat(1, 1) = matForPatchSecDerivs(row, 2);

    secFFDerivPreFacMat(0, 0) += normalDerivatives(0);
    secFFDerivPreFacMat(0, 1) += normalDerivatives(1);
    secFFDerivPreFacMat(1, 1) += normalDerivatives(2);

    secFFDerivPreFacMat(1, 0) = secFFDerivPreFacMat(0, 1);

    return -initArea * (energyDensityDerivWRTSecFF * secFFDerivPreFacMat).trace() * faceNormal;
}

Eigen::Matrix<double, 3, 1> Triangle::getBendingForcePatch(int row, const Eigen::Vector3d& faceNormal, const Eigen::Matrix<double, 2, 2>& energyDensityDerivWRTSecFF) const {
    Eigen::Matrix<double, 2, 2> secFFDerivPreFacMat;

    secFFDerivPreFacMat(0, 0) = matForPatchSecDerivs(row, 0);
    secFFDerivPreFacMat(0, 1) = matForPatchSecDerivs(row, 1);
    secFFDerivPreFacMat(1, 1) = matForPatchSecDerivs(row, 2);

    secFFDerivPreFacMat(1, 0) = secFFDerivPreFacMat(0, 1);

    return -initArea * (energyDensityDerivWRTSecFF * secFFDerivPreFacMat).trace() * faceNormal;
}

void Triangle::updateGeometricProperties(double bending_pre_factor, double j_pre_factor, double poisson_ratio, double stretching_prefactor) {
    const Eigen::Vector3d p0 = *corner_nodes_pos[0];
    const Eigen::Vector3d p1 = *corner_nodes_pos[1];
    const Eigen::Vector3d p2 = *corner_nodes_pos[2];

    currSides.col(0).noalias() = p1 - p0;
    currSides.col(1).noalias() = p2 - p0;
    Eigen::Vector3d faceNormal = currSides.col(0).cross(currSides.col(1));
    defGradient.noalias() = currSides * invInitSidesMat;

    met.noalias() = defGradient.transpose() * defGradient;
    metInv.noalias() = met.inverse();
    metInvDet = metInv.determinant();

    currAreaInv = 2 / faceNormal.norm();
    faceNormal *= 0.5 * currAreaInv; // Normalising

    Eigen::Matrix<double, 3, 3> patchSecDerivs = p0 * matForPatchSecDerivs.row(0)
        + p1 * matForPatchSecDerivs.row(1)
        + p2 * matForPatchSecDerivs.row(2)
        + *patch_nodes_pos[0] * matForPatchSecDerivs.row(3)
        + *patch_nodes_pos[1] * matForPatchSecDerivs.row(4)
        + *patch_nodes_pos[2] * matForPatchSecDerivs.row(5);

    centroid.noalias() = (p0 + p1 + p2) / 3;

    Eigen::Vector3d vectorOfSecFFComps = patchSecDerivs.transpose() * faceNormal;

    secFF(0, 0) = vectorOfSecFFComps(0);
    secFF(0, 1) = vectorOfSecFFComps(1);
    secFF(1, 0) = vectorOfSecFFComps(1);
    secFF(1, 1) = vectorOfSecFFComps(2);

    Eigen::Matrix<double, 2, 2> relativeSecFF = programmedMetInv * (secFF - programmedSecFF);
    double areaMultiplier = bending_pre_factor * programmedMetInvDet;
    double J = areaMultiplier * j_pre_factor;

    // Now calculate bending energy density for this
    double preGentBendEnergyDensity = areaMultiplier * ((1 - poisson_ratio) * (relativeSecFF * relativeSecFF).trace() +
        poisson_ratio * relativeSecFF.trace() * relativeSecFF.trace());
    double gentDerivFac = (1 + 2 * preGentBendEnergyDensity / J);

    /* Calculate the derivative of the bending energy density with respect to the secFF. */
    Eigen::Matrix<double, 2, 2> energyDensityDerivWRTSecFF = gentDerivFac * 2 * areaMultiplier *
    ((1 - poisson_ratio) * relativeSecFF * programmedMetInv +
        poisson_ratio * relativeSecFF.trace() * programmedMetInv);
    bendEnergyDensity = gentDerivFac * preGentBendEnergyDensity;

    Eigen::Matrix<double, 3, 3> stretchForces = getStretchingForces(stretching_prefactor);

    Eigen::Matrix<double, 3, 3> triangleEdgeNormals = getTriangleEdgeNormals(currSides, faceNormal);
    Eigen::Matrix<double, 3, 3> normalDerivPiece =
        0.5 * currAreaInv * (patchSecDerivs.transpose() * triangleEdgeNormals);

    for (int n = 0; n < 3; ++n) {
        node_triangle_force[n]->noalias() = getBendingForceNode(normalDerivPiece.col(n), n, faceNormal, energyDensityDerivWRTSecFF) + stretchForces.col(n);
        node_triangle_force[n + 3]->noalias() = getBendingForcePatch(n + 3, faceNormal, energyDensityDerivWRTSecFF);
    }
}


void Triangle::updateAngleDeficits(std::vector<double>& angleDeficits) const {
    std::vector<Eigen::Vector3d> sides;
    sides.emplace_back(currSides.col(0));
    sides.emplace_back(currSides.col(1));
    sides.emplace_back(sides[1] - sides[0]);

    std::vector<double> sidesLength(sides.size());
    for (auto& side : sides) { sidesLength.emplace_back(side.norm()); }

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

template <typename T> T interpolate(const T& previous, const T& next, double ratio) {
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

void Triangle::updateProgrammedMetricImplicit(int stage_counter, double dial_in_factor) {
    Eigen::Vector3d metric_current =
        interpolate(programmed_metric_info, next_programmed_metric_info, dial_in_factor);

    double dirAng = metric_current(0);
    double lambda = metric_current(1);
    double nu = metric_current(2);

    updateProgrammedMetricImplicit(dirAng, lambda, nu);
}

void Triangle::updateProgrammedTensorsDynamically(int stage_counter, double dial_in_factor, double transfer_coefficient,
                                                  double min_height, double max_height) {
    Eigen::Vector3d metric_current =
        interpolate(programmed_metric_info, next_programmed_metric_info, dial_in_factor);
    Eigen::Matrix<double, 2, 2> bend_programmed =
        interpolate(programmed_second_fundamental_form, next_programmed_second_fundamental_form,
                    dial_in_factor);

    double dirAng = metric_current(0);
    double lambda = metric_current(1);
    double nu = metric_current(2);

    double normalized_height = (getHeight() - min_height) / (max_height - min_height);
    if (fabs(max_height - min_height) < 1e-6) { normalized_height = 0; }

    local_magnitude = interpolate(local_magnitude, normalized_height, transfer_coefficient);
    local_elongation = 1 - (1 - lambda) * local_magnitude;

    updateProgrammedMetricImplicit(dirAng, local_elongation, nu);
    programmedSecFF = bend_programmed * (1 - local_magnitude);
    //    programmedSecFF = bend_programmed * local_magnitude;
}

void Triangle::updateProgrammedMetricExplicit(int stage_counter, double dial_in_factor) {
    programmedMetInv = interpolate(programmed_metric_inv, next_programmed_metric_inv, dial_in_factor);
    programmedMetInvDet = programmedMetInv.determinant();
}

void Triangle::updateProgrammedSecondFundamentalForm(int stage_counter, double dial_in_factor_root) {
    programmedSecFF = interpolate(programmed_second_fundamental_form,
                                  next_programmed_second_fundamental_form, dial_in_factor_root);
}

void Triangle::updateProgrammedTaus(int stage_counter, double dial_in_factor) {
    dialledProgTau = interpolate(programmed_tau, next_programmed_tau, dial_in_factor);
}

void Triangle::updateProgrammedQuantities(const int stage_counter, const double dial_in_factor,
                                          const double dial_in_factor_root,
                                          const bool is_lce_metric_used, const bool is_elongation_dynamically_updated,
                                          const double transfer_coefficient, const double min_height,
                                          const double max_height) {
    if (!is_lce_metric_used) {
        updateProgrammedMetricExplicit(stage_counter, dial_in_factor);
        updateProgrammedSecondFundamentalForm(stage_counter, dial_in_factor_root);
    } else if (!is_elongation_dynamically_updated) {
        updateProgrammedMetricImplicit(stage_counter, dial_in_factor);
        updateProgrammedSecondFundamentalForm(stage_counter, dial_in_factor_root);
    } else {
        updateProgrammedTensorsDynamically(stage_counter, dial_in_factor, transfer_coefficient, min_height,
                                           max_height);
    }
    updateProgrammedTaus(stage_counter, dial_in_factor);
}


std::vector<unsigned int> getNeighbouringNodes(const Triangle& current_triangle, const std::vector<Triangle>& triangles,
                                               const std::vector<Node>& nodes) {
    // Provides nodes, which contain edges that neighbour the nodes of the current triangle, exclusive of those of
    // considered triangle
    std::unordered_set<unsigned int> neighbouring_triangle_indices;
    std::unordered_set<unsigned int> neighbouring_node_indices;

    for (unsigned int i = 0; i < 3; i++) {
        unsigned int i_node = current_triangle.vertexLabels(i);
        const Node& node = nodes[i_node];
        for (unsigned int j_triangles : node.incidentTriLabels) { neighbouring_triangle_indices.insert(j_triangles); }
    }

    for (unsigned int j_triangles : neighbouring_triangle_indices) {
        const Triangle& triangle = triangles[j_triangles];
        for (unsigned int j : triangle.vertexLabels) { neighbouring_node_indices.insert(j); }
    }

    // We want to exclude the labels of the original nodes
    for (unsigned int i : current_triangle.vertexLabels) { neighbouring_node_indices.erase(i); }

    std::vector<unsigned int> v_node_indices;
    for (unsigned int index : neighbouring_node_indices) { v_node_indices.emplace_back(index); }
    return v_node_indices;
}

std::vector<std::pair<unsigned int, double>> assignDistanceFromCentroid(const std::vector<unsigned int>& node_indices,
                                                                        const std::vector<Node>& nodes,
                                                                        const Eigen::Vector3d& position) {
    std::vector<std::pair<unsigned int, double>> index_distance_pair;
    for (unsigned int index : node_indices) {
        double distance = (nodes[index].pos - position).norm();
        index_distance_pair.emplace_back(index, distance);
    }
    return index_distance_pair;
}

Eigen::Matrix<double, 6, 1> patchColumn(const Eigen::Vector3d& position, const Eigen::Vector3d& centroid) {
    Eigen::Matrix<double, 6, 1> patchColumn;
    patchColumn(0) = 1;
    patchColumn(1) = (position(0) - centroid(0));
    patchColumn(2) = (position(1) - centroid(1));
    patchColumn(3) = 0.5 * patchColumn(1) * patchColumn(1);
    patchColumn(4) = patchColumn(1) * patchColumn(2);
    patchColumn(5) = 0.5 * patchColumn(2) * patchColumn(2);
    return patchColumn;
}

double Triangle::getHeight() const { return centroid(2); }

double Triangle::updateMatForPatchDerivs(const std::vector<Triangle>& triangles, const std::vector<Node>& nodes) {
    Eigen::Vector3d refCentroid = (nodes[vertexLabels(0)].pos +
        nodes[vertexLabels(1)].pos +
        nodes[vertexLabels(2)].pos) / 3;

    std::vector<unsigned int> possiblePatchNodeLabels = getNeighbouringNodes(*this, triangles, nodes);
    std::vector<std::pair<unsigned int, double>> indexDistancePairs = assignDistanceFromCentroid(
        possiblePatchNodeLabels, nodes, refCentroid);

    std::sort(std::begin(indexDistancePairs), std::end(indexDistancePairs),
              [](const std::pair<unsigned int, double>& p1,
                 const std::pair<unsigned int, double>& p2) -> bool {
                  return p1.second < p2.second;
              });

    Eigen::Matrix<double, 6, 6> patchNodeDataMatrix;
    double inner_patch_size = (nodes[vertexLabels(0)].pos - refCentroid).squaredNorm() +
        (nodes[vertexLabels(1)].pos - refCentroid).squaredNorm() +
        (nodes[vertexLabels(2)].pos - refCentroid).squaredNorm();

    for (int n = 0; n < 3; ++n) {
        Eigen::Vector3d candidatePatchNode;
        candidatePatchNode = nodes[vertexLabels(n)].pos;
        patchNodeDataMatrix.col(n) = patchColumn(candidatePatchNode, refCentroid);
    }

    std::set<std::vector<int>> candidate_trios;
    for (unsigned int p : possiblePatchNodeLabels) {
        for (unsigned int q : possiblePatchNodeLabels) {
            for (unsigned int r : possiblePatchNodeLabels) { candidate_trios.insert({(int)p, (int)q, (int)r}); }
        }
    }

    auto patch_condition_number = DBL_MAX;
    for (auto& candidateIndices : candidate_trios) {
        patchNodeDataMatrix.col(3) = patchColumn(nodes[candidateIndices[0]].pos, refCentroid);
        patchNodeDataMatrix.col(4) = patchColumn(nodes[candidateIndices[1]].pos, refCentroid);
        patchNodeDataMatrix.col(5) = patchColumn(nodes[candidateIndices[2]].pos, refCentroid);

        Eigen::FullPivLU<Eigen::Matrix<double, 6, 6>> patchNodeDecomposition;
        patchNodeDecomposition.compute(patchNodeDataMatrix);
        if (!patchNodeDecomposition.isInvertible()) { continue; }

        double outer_patch_size = (nodes[candidateIndices[0]].pos - refCentroid).squaredNorm() +
            (nodes[candidateIndices[1]].pos - refCentroid).squaredNorm() +
            (nodes[candidateIndices[2]].pos - refCentroid).squaredNorm();

        double patch_size = sqrt((inner_patch_size + outer_patch_size) / 6);
        Eigen::Matrix<double, 6, 6> invTempPatchNodeDataMatrix = patchNodeDataMatrix.inverse();
        Eigen::Matrix<double, 6, 3> candidatePatchDiv;
        candidatePatchDiv = invTempPatchNodeDataMatrix.block<6, 3>(0, 3);

        Eigen::JacobiSVD<Eigen::Matrix<double, 6, 3>> secDerivMatTempSVD;
        secDerivMatTempSVD.compute(candidatePatchDiv);

        double singular_values = secDerivMatTempSVD.singularValues()(0); // This has dimensions 1 / Length ^ 2.
        double current_condition_number = singular_values * pow(patch_size, 2);

        if (current_condition_number < patch_condition_number) {
            patch_condition_number = current_condition_number;
            nonVertexPatchNodesLabels[0] = candidateIndices[0];
            nonVertexPatchNodesLabels[1] = candidateIndices[1];
            nonVertexPatchNodesLabels[2] = candidateIndices[2];

            patch_nodes_pos[0] = &nodes[candidateIndices[0]].pos;
            patch_nodes_pos[1] = &nodes[candidateIndices[1]].pos;
            patch_nodes_pos[2] = &nodes[candidateIndices[2]].pos;

            matForPatchSecDerivs = candidatePatchDiv;
        }
    }
    return patch_condition_number;
}


void Triangle::setLocalElongation(double local_elongation) { Triangle::local_elongation = local_elongation; }

Triangle::Triangle(int label, int id_0, int id_1, int id_2, const std::vector<Node>& nodes) : Triangle() {
    Triangle::label = label;
    vertexLabels(0) = id_0;
    vertexLabels(1) = id_1;
    vertexLabels(2) = id_2;
    corner_nodes_pos[0] = &nodes[id_0].pos;
    corner_nodes_pos[1] = &nodes[id_1].pos;
    corner_nodes_pos[2] = &nodes[id_2].pos;
    reference_node_positions[0] = nodes[id_0].pos;
    reference_node_positions[1] = nodes[id_1].pos;
    reference_node_positions[2] = nodes[id_2].pos;
}

void Triangle::updateMagneticForce(const Eigen::Vector3d& magnetic_field) {
    for (int i = 0; i < 3; i++) {
        Eigen::Vector3d magnetic_force = {0, 0, 0};
        magnetic_force[0] = magnetic_field[0] * reference_magnetisation_density[1] * (
                reference_node_positions[(i + 1) % 3](0) - reference_node_positions[(i + 2) % 3](0)) -
            magnetic_field[0] * reference_magnetisation_density[0] * (
                reference_node_positions[(i + 1) % 3](1) - reference_node_positions[(i + 2) % 3](1)) +
            magnetic_field[1] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(2) - (*corner_nodes_pos[(i + 2) % 3])(2)) -
            magnetic_field[2] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(1) - (*corner_nodes_pos[(i + 2) % 3])(1));

        magnetic_force[1] = magnetic_field[1] * reference_magnetisation_density[1] * (
                reference_node_positions[(i + 1) % 3](0) - reference_node_positions[(i + 2) % 3](0)) -
            magnetic_field[1] * reference_magnetisation_density[0] * (
                reference_node_positions[(i + 1) % 3](1) - reference_node_positions[(i + 2) % 3](1)) +
            magnetic_field[2] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(0) - (*corner_nodes_pos[(i + 2) % 3])(0)) -
            magnetic_field[0] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(2) - (*corner_nodes_pos[(i + 2) % 3])(2));

        magnetic_force[2] = magnetic_field[2] * reference_magnetisation_density[1] * (
                reference_node_positions[(i + 1) % 3](0) - reference_node_positions[(i + 2) % 3](0)) -
            magnetic_field[2] * reference_magnetisation_density[0] * (
                reference_node_positions[(i + 1) % 3](1) - reference_node_positions[(i + 2) % 3](1)) +
            magnetic_field[0] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(1) - (*corner_nodes_pos[(i + 2) % 3])(1)) -
            magnetic_field[1] * reference_magnetisation_density[2] * (
                (*corner_nodes_pos[(i + 1) % 3])(0) - (*corner_nodes_pos[(i + 2) % 3])(0));

        node_triangle_force[i]->noalias() = -magnetic_force;
    }
    for (int i = 3; i < 6; i++) { node_triangle_force[i]->noalias() = Eigen::Vector3d({0, 0, 0}); }
}

void Triangle::setReferenceMagnetisationDensity(const Eigen::Vector3d& magnetisation) {
    reference_magnetisation_density = magnetisation;
}

void Triangle::setNodeForceAddress(unsigned int index, Eigen::Vector3d* address) {
    node_triangle_force[index] = address;
}
