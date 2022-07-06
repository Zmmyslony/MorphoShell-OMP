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

This is the header file for the class that will contain data for each
triangular element, such as vertices, area etc.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef _TRIANGLE_CLASS_TAG_
#define _TRIANGLE_CLASS_TAG_

#include <vector>
#include <Eigen/Dense>
#include "Node.hpp"
#include "SettingsStruct.hpp"

#include "CustomOutStreamClass.hpp"

class Triangle {
public:
    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
    CustomOutStreamClass triLogStream;

    // Label so this triangle 'knows' which it is.
    int label;

    /* Boolean representing whether the triangle has a vertex on the boundary of
    the sample (true) or not (false).*/
    bool isOnBoundary;

    // Initial (reference) area and 1/(current area) for this triangle.
    double initArea;
    double invCurrArea;

    /* Indices (labels) of the nodes at the vertices of the triangle. These are
    in no particular order*/
    Eigen::Vector3i vertexLabels;

    /* Labels (and indices in the triangles' container vector) of the 3 edges of
    this triangle.*/
    Eigen::Vector3i edgeLabels;

    /* Vector storing, for each non-boundary edge, its initial
    length divided by sum of initial non-boundary edge lengths for the triangle.
    These fractions are used for weighting the edge normals used in the secFF
    calculation.*/
    Eigen::VectorXd initNonBoundEdgeLengthFracs;

    /* Labels (and indices in the triangles' container vector) of the other
    triangles that share an edge with this triangle.*/
    Eigen::VectorXi edgeSharingTriLabels;

    /* A vector with one component for each non-boundary edge of this triangle,
    with each entry being either +1.0 or -1.0. A value of +1.0 means that
    this triangle corresponds to adjTriLabels(0) for the corresponding edge.
    -1.0 implies adjTriLabels(1) similarly. */
    Eigen::VectorXi edgeAdjTriLabelSelectors;

    /* Vector where the i'th element is the index corresponding to this triangle
    in the edgeSharingTriLabels vector of the i'th edge-sharing triangle
    (neighbour) of THIS triangle.*/
    Eigen::VectorXi indicesIntoEdgeSharingTriLabelsOfNeighbours;

    /* Indices (labels) of the 3 nodes that are not vertices of the triangle, but
    are part of the estimation of the 2nd F.F. for each triangle. For
    non-boundary triangles these are the non-shared-edge nodes from each of the
    triangles sharing an edge with this triangle. For boundary triangles one or
    two of these are missing, and different nodes are chosen based on proximity
    to this triangle's centroid. There is again no particular order.*/
    Eigen::Vector3i nonVertexPatchNodesLabels;

    /* Position of triangle's centroid in initial (reference) x-y plane.*/
    Eigen::Vector3d refCentroid;

    /* Matrix (of doubles) where each *column* is a vector describing a side of
    the triangle. We only need two sides per triangle for the algorithm
    hence the 3x2 matrices.*/
    Eigen::Matrix<double, 3, 2> currSides;

    // Current unit normal to face.
    Eigen::Vector3d faceNormal;

    /* Inverse of 2x2 matrix that has (two) initial sides of the triangle as
    columns. Those sides correspond to the two current sides stored in currSides.*/
    Eigen::Matrix<double, 2, 2> invInitSidesMat;

    /* The reference (initial) state in-triangle-plane *outward* normals of
    each triangle's sides (with lengths = corresponding side lengths).
    initOutwardTriNormals.col(v) is opposite vertexLabels(v).*/
    Eigen::Matrix<double, 2, 3> initOutwardSideNormals;

    /* Matrix representing the inverse of the *energetically* favoured metric for
    this triangle, induced by a programmed nematic director field, for example.
    The matrix stored here is the 'dialled in' value, which is used to prevent
    anything too explosive happening. */
    Eigen::Matrix<double, 2, 2> dialledInvProgMetric;

    // Determinant of the above dialledInvProgMetric matrix.
    double detDialledInvProgMetric;

    /* Dialled in programmed scalar 'tau' factor.*/
    double dialledProgTau;

    /* Matrix representing the components of the *energetically* favoured
    (programmed) Second Fundamental Form in the x-y cartesian coordinate system
    of the initial flat state. */
    Eigen::Matrix<double, 2, 2> dialledProgSecFF;

    /* Matrix representing the total deformation gradient of this triangle. This
    maps the in-plane reference (initial) state (for which no z components are
    stored) to a triangle in 3D space, so it is 3x2.*/
    Eigen::Matrix<double, 3, 2> defGradient;

    /* The (1st order) approximation for the metric for this triangle, which
    equals defGradient.transpose9) * defGradient.*/
    Eigen::Matrix<double, 2, 2> metric;

    /* Inverse of metric.*/
    Eigen::Matrix<double, 2, 2> invMetric;

    /* Determinant of the inverse of the metric.*/
    double detInvMetric;

    /* Matrix representing (1st Piola-Kirchoff stress tensor)/2 for this triangle.*/
    Eigen::Matrix<double, 3, 2> halfPK1Stress;

    /*Matrix that is pre-calculated and then used repeatedly in finding the
    components of the second fundamental form estimated for this triangle. */
    Eigen::Matrix<double, 6, 3> matForPatchSecDerivs;

    /* Matrix of the second position derivatives, used in calculating the
    secFF estimate. */
    Eigen::Matrix<double, 3, 3> patchSecDerivs;

    /* Estimated Second Fundamental Form matrix (secFF) of the deformed surface,
    defined (as with the deformation gradient) with respect to the 'material'
    coordinate system, i.e. the coordinate chart that used to be the (x,y)
    cartesians of the flat initial state, and then deformed with the sheet. This
    is a 2x2 symmetric matrix. */
    Eigen::Matrix<double, 2, 2> secFF;

    /* Derivative of the bending energy density with respect to the secFF.*/
    Eigen::Matrix<double, 2, 2> energyDensityDerivWRTSecFF;

    /* The derivative of the bending energy density with respect to the metric.*/
    Eigen::Matrix<double, 2, 2> bendEnergyDensityDerivWRTMetric;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for logStream).
    I ensure that default data values are recognisable values,
    for debugging.
    edgeSharingTriLabels, nonSharedVerticesOfEdgeSharingTris,
    edgeAdjTriLabelSelectors, indicesIntoEdgeSharingTriLabelsOfNeighbours,
    and usefulTermsForSecFFDeriv
    are left with zero size at initialisation. */
    Triangle() {
        label = -12345;
        isOnBoundary = false;
        vertexLabels.fill(123456789);
        nonVertexPatchNodesLabels.fill(987654321);
        edgeLabels.fill(-4321);
        initArea = 98765;
        invCurrArea = 56789;
        currSides.fill(5678);
        faceNormal.fill(8765);
        initOutwardSideNormals.fill(8765432);
        invInitSidesMat.fill(-5678);
        dialledInvProgMetric.fill(-6789);
        detDialledInvProgMetric = -9.87654321;
        dialledProgTau = 123456;
        dialledProgSecFF.fill(-5678);
        defGradient.fill(123456);
        metric.fill(1234567.89);
        detInvMetric = -1234.567;
        invMetric.fill(-1234567.89);
        halfPK1Stress.fill(-1234567);
        patchSecDerivs.fill(12345.6789);
        secFF.fill(-7654);
        energyDensityDerivWRTSecFF.fill(456789);
        bendEnergyDensityDerivWRTMetric.fill(-456789);
        matForPatchSecDerivs.fill(54321);
    }

    // Declare other member functions.

    // Debugging function to display all member data.
    void display();

    Eigen::Matrix<double, 3, 2> updateHalfPK1Stress(double stretchingPrefactor);

    Eigen::Matrix<double, 3, 3> getStretchingForces();

    Eigen::Matrix<double, 3, 3> getOutwardTriangleNormals();

    Eigen::Matrix<double, 3, 1> getBendingForce(const Eigen::Matrix<double, 3, 3> &normalDerivatives, int row);

    void updateMetric(const std::vector<Node> &nodes);

    void updateGeometricProperties(const std::vector<Node> &nodes);

    void calculateSecondFundamentalForm(double bendingPreFac, double JPreFactor);

    void updateAngleDeficits(std::vector<double> &angleDeficits) const;
};

#endif
