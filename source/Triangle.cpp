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
    triLogStream << "invCurrArea = " << invCurrArea << std::endl;
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
    triLogStream << "Dialled in inverse of programmed metric = " << "\n" << dialledInvProgMetric << std::endl;
    triLogStream << "detDialledInvProgMetric = " << "\n" << detDialledInvProgMetric << std::endl;
    triLogStream << "Dialled in prog tau factor = " << dialledProgTau << std::endl;
    triLogStream << "Dialled in programmed second fundamental form = " << "\n" << dialledProgSecFF << std::endl;
    triLogStream << "Deformation gradient = " << "\n" << defGradient << std::endl;
    triLogStream << "metric = " << "\n" << metric << std::endl;
    triLogStream << "Inverse of Metric = " << "\n" << invMetric << std::endl;
    triLogStream << "Det of inverse of Metric = " << "\n" << detInvMetric << std::endl;
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
                                   (dialledInvProgMetric - (detInvMetric / detDialledInvProgMetric) * invMetric));
    return halfPK1Stress;
}

Eigen::Matrix<double, 3, 3> Triangle::getStretchingForces() {
    return halfPK1Stress * initOutwardSideNormals;
}

Eigen::Matrix<double, 3, 3> Triangle::getOutwardTriangleNormals() {
    Eigen::Matrix<double, 3, 3> outwardTriangleNormals;

    outwardTriangleNormals.col(1) = faceNormal.cross(currSides.col(1));
    outwardTriangleNormals.col(2) = -faceNormal.cross(currSides.col(0));
    outwardTriangleNormals.col(0) = -outwardTriangleNormals.col(1) - outwardTriangleNormals.col(2);
    return outwardTriangleNormals;
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

    // Calculate corresponding metric, and its det and inverse.
    metric = defGradient.transpose() * defGradient;
    detInvMetric = 1 / metric.determinant();
//    Eigen::Matrix<double, 2, 2> metricAdjMatrix;
//    metricAdjMatrix << metric(1, 1), -metric(0, 1),
//            -metric(0, 1), metric(0, 0);
//    invMetric = detInvMetric * metricAdjMatrix;
    invMetric = metric.inverse();
}

void Triangle::updateGeometricProperties(const std::vector<Node> &nodes) {
    Eigen::Matrix<double, 3, 6> matrixOfPatchNodeCoords;
    // Loop over patch nodes for this triangle and get their positions.
    for (int n = 0; n < 6; ++n) {
        if (n < 3) {
            matrixOfPatchNodeCoords.col(n) = nodes[vertexLabels(n)].pos;
        } else {
            matrixOfPatchNodeCoords.col(n) = nodes[nonVertexPatchNodesLabels(n - 3)].pos;
        }
    }

    patchSecDerivs = matrixOfPatchNodeCoords * matForPatchSecDerivs;

    faceNormal = currSides.col(0).cross(currSides.col(1));
    invCurrArea = 2.0 / faceNormal.norm();
    faceNormal = 0.5 * faceNormal * invCurrArea; // Normalising
}


void Triangle::calculateSecondFundamentalForm(double bendingPreFac, double JPreFactor) {
    Eigen::Vector3d vectorOfSecFFComps = patchSecDerivs.transpose() * faceNormal;

    secFF(0, 0) = vectorOfSecFFComps(0);
    secFF(0, 1) = vectorOfSecFFComps(1);
    secFF(1, 0) = vectorOfSecFFComps(1);
    secFF(1, 1) = vectorOfSecFFComps(2);

    Eigen::Matrix<double, 2, 2> tempMat1 =
            dialledInvProgMetric * (secFF - dialledProgSecFF);
    Eigen::Matrix<double, 2, 2> tempMat2 = tempMat1 * dialledInvProgMetric;
    double tr_tempMat1 = tempMat1.trace();
    double tempScalar = bendingPreFac * detDialledInvProgMetric;

    double J = tempScalar * JPreFactor;

    // Now calculate bending energy density for this 
    double preGentBendEnergyDensity = tempScalar * ((tempMat1 * tempMat1).trace() + tr_tempMat1 * tr_tempMat1);
    double gentDerivFac = (1 + 2 * preGentBendEnergyDensity / J);

    /* Calculate the derivative of the bending energy density with respect
    to the secFF. */
    energyDensityDerivWRTSecFF =
            gentDerivFac * 2 * tempScalar * (tempMat2 + tr_tempMat1 * dialledInvProgMetric);
}


void Triangle::updateAngleDeficits(std::vector<double> &angleDeficits) const {
    std::vector<Eigen::Vector3d> sides;
    sides.emplace_back(currSides.col(0));
    sides.emplace_back(currSides.col(1));
    sides.emplace_back(sides[1] - sides[0]);

    std::vector<double> sidesLength(sides.size());
    for (auto &side : sides) {
        sidesLength.emplace_back(side.norm());
    }

    angleDeficits[vertexLabels(0)] -= acos((sides[0].dot(sides[1]) / (sidesLength[0] * sidesLength[1]));
    angleDeficits[vertexLabels(1)] -= acos((sides[0].dot(sides[2]) / (sidesLength[0] * sidesLength[2]));
    angleDeficits[vertexLabels(0)] -= acos((sides[1].dot(sides[2]) / (sidesLength[1] * sidesLength[2]));
}