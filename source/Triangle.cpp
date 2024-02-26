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
#include <Eigen/Dense>

#include "Triangle.hpp"


//This is a debugging tool to display the node's data
void Triangle::display() {
    triLogStream.open();
    triLogStream << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    triLogStream << "Triangle " << label << ":" << std::endl;
    triLogStream << "Boundary indicator = " << isOnBoundary << std::endl;
    triLogStream << "Initial (reference) area = " << initArea << std::endl;
    triLogStream << "invCurrArea = " << currAreaInv << std::endl;
    triLogStream << "Labels of vertices: " << vertexLabels.transpose() << std::endl;
    triLogStream << "Labels of edges: " << edgeLabels.transpose() << std::endl;
    triLogStream << "Initial non-boundary edge length fractions: " << initNonBoundEdgeLengthFracs.transpose()
                 << std::endl;
    triLogStream << "Labels of adjacent (edge-sharing) triangles: " << edgeSharingTriLabels.transpose() << std::endl;
    triLogStream << "edgeAdjTriLabelSelectors: " << edgeAdjTriLabelSelectors.transpose() << std::endl;
    triLogStream << "indicesIntoEdgeSharingTriLabelsOfNeighbours: "
                 << indicesIntoEdgeSharingTriLabelsOfNeighbours.transpose() << std::endl;
    triLogStream << "Labels of non-vertex nodes in this triangle's patch: " << nonVertexPatchNodesLabels.transpose()
                 << std::endl;
    triLogStream << "Current sides = " << "\n" << currSides << std::endl;
    triLogStream << "faceNormal = " << "\n" << faceNormal << std::endl;
    triLogStream << "Initial outward side normals = " << "\n" << initOutwardSideNormals << std::endl;
    triLogStream << "invInitInPlaneSidesMat = " << "\n" << invInitSidesMat << std::endl;
    triLogStream << "Dialled in inverse of programmed metric = " << "\n" << programmedMetInv << std::endl;
    triLogStream << "detDialledInvProgMetric = " << "\n" << programmedMetInvDet << std::endl;
    triLogStream << "Dialled in prog tau factor = " << dialledProgTau << std::endl;
    triLogStream << "Dialled in programmed second fundamental form = " << "\n" << programmedSecFF << std::endl;
    triLogStream << "Deformation gradient = " << "\n" << defGradient << std::endl;
    triLogStream << "metric = " << "\n" << met << std::endl;
    triLogStream << "Inverse of Metric = " << "\n" << metInv << std::endl;
    triLogStream << "Det of inverse of Metric = " << "\n" << metInvDet << std::endl;
    triLogStream << "matForPatchSecDerivs = " << "\n" << matForPatchSecDerivs << std::endl;
    triLogStream << "patchSecDerivs = " << "\n" << patchSecDerivs << std::endl;
    triLogStream << "Second fundamental form = " << "\n" << secFF << std::endl;
    triLogStream << "Bending energy density deriv wrt secFF = " << "\n" << energyDensityDerivWRTSecFF << std::endl;
    triLogStream << "bendEnergyDensityDerivWRTMetric = " << "\n" << bendEnergyDensityDerivWRTMetric << std::endl;
    triLogStream << "halfPK1Stress = " << "\n" << halfPK1Stress << std::endl;
    triLogStream << "-----------------------------" << std::endl;
    triLogStream.close();
}

Eigen::Matrix<double, 3, 2> Triangle::updateHalfPK1Stress(double stretchingPrefactor) {
    halfPK1Stress = defGradient * (stretchingPrefactor * dialledProgTau *
                                   (programmedMetInv - (metInvDet / programmedMetInvDet) * metInv));
    return halfPK1Stress;
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
    currSides.col(0) = nodes[vertexLabels(1)].pos - nodes[vertexLabels(0)].pos;
    currSides.col(1) = nodes[vertexLabels(2)].pos - nodes[vertexLabels(0)].pos;
    defGradient = currSides * invInitSidesMat;

    met = defGradient.transpose() * defGradient;
    metInvDet = 1 / met.determinant();
    metInv = met.inverse();
}

void Triangle::updateGeometricProperties(const std::vector<Node> &nodes) {
    updateMetric(nodes);
    Eigen::Matrix<double, 3, 6> matrixOfPatchNodeCoords;
    for (int n = 0; n < 6; ++n) {
        if (n < 3) {
            matrixOfPatchNodeCoords.col(n) = nodes[vertexLabels(n)].pos;
        } else {
            matrixOfPatchNodeCoords.col(n) = nodes[nonVertexPatchNodesLabels(n - 3)].pos;
        }
    }
    patchSecDerivs = matrixOfPatchNodeCoords * matForPatchSecDerivs;

    faceNormal = currSides.col(0).cross(currSides.col(1));
    currAreaInv = 2 / faceNormal.norm();
    faceNormal = 0.5 * faceNormal * currAreaInv; // Normalising
}


void Triangle::updateSecondFundamentalForm(double bendingPreFac, double JPreFactor, double poissonRatio) {
    Eigen::Vector3d vectorOfSecFFComps = patchSecDerivs.transpose() * faceNormal;

    secFF(0, 0) = vectorOfSecFFComps(0);
    secFF(0, 1) = vectorOfSecFFComps(1);
    secFF(1, 0) = vectorOfSecFFComps(1);
    secFF(1, 1) = vectorOfSecFFComps(2);

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
                            3.0 / dialledProgTau);
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

void Triangle::updateProgrammedMetricFromLCEInfo(int stage_counter, double dial_in_factor) {
    Eigen::Vector3d metric_previous = programmed_metric_infos[stage_counter];
    Eigen::Vector3d metric_next = programmed_metric_infos[stage_counter + 1];
    double dirAng = (1.0 - dial_in_factor) * metric_previous(0) + dial_in_factor * metric_next(0);
    double cosDirAng = cos(dirAng);
    double sinDirAng = sin(dirAng);
    double lambda = (1.0 - dial_in_factor) * metric_previous(1) + dial_in_factor * metric_next(1);
    double nu = (1.0 - dial_in_factor) * metric_previous(2) + dial_in_factor * metric_next(2);
    double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
    double lambdaToThe2Nu = pow(lambda, 2.0 * nu);

    programmedMetInv(0, 0) =
            lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
    programmedMetInv(0, 1) = (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
    programmedMetInv(1, 0) = programmedMetInv(0, 1);
    programmedMetInv(1, 1) =
            lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
}

void Triangle::updateProgrammedMetric(int stage_counter, double dial_in_factor) {
    Eigen::Matrix<double, 2, 2> metric_previous = programmed_metric_inv[stage_counter];
    Eigen::Matrix<double, 2, 2> metric_next = programmed_metric_inv[stage_counter + 1];
    programmedMetInv = (1.0 - dial_in_factor) * metric_previous + dial_in_factor * metric_next;
    programmedMetInvDet = programmedMetInv.determinant();
}

void Triangle::updateProgrammedSecondFundamentalForm(int stage_counter, double dial_in_factor_root) {
    Eigen::Matrix<double, 2, 2> fundamental_form_previous = programmed_second_fundamental_form[stage_counter];
    Eigen::Matrix<double, 2, 2> fundamental_form_next = programmed_second_fundamental_form[stage_counter + 1];

    programmedSecFF =
            (1.0 - dial_in_factor_root) * fundamental_form_previous +
            dial_in_factor_root * fundamental_form_next;
}

void Triangle::updateProgrammedTaus(int stage_counter, double dial_in_factor) {
    double tau_previous = programmed_taus[stage_counter];
    double tau_next = programmed_taus[stage_counter + 1];
    dialledProgTau = (1.0 - dial_in_factor) * tau_previous + dial_in_factor * tau_next;
}

void Triangle::updateProgrammedQuantities(int stage_counter, double dial_in_factor, double dial_in_factor_root, bool is_LCE_metric_used) {
    if (is_LCE_metric_used) {
        updateProgrammedMetricFromLCEInfo(stage_counter, dial_in_factor);
    } else {
        updateProgrammedMetric(stage_counter, dial_in_factor);
    }
    updateProgrammedSecondFundamentalForm(stage_counter, dial_in_factor_root);
    updateProgrammedTaus(stage_counter, dial_in_factor);
}
