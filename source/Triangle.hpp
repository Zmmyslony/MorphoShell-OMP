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
// #include "Settings.hpp"


class Triangle {
    double local_elongation = 1;
    const Node *corner_nodes[3];
    const Node *patch_nodes[3];
    Eigen::Vector3d node_elastic_force[6];

public:
    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/

    /// Label so this triangle 'knows' which it is.
    int label;

    /* Boolean representing whether the triangle has a vertex on the boundary of
    the sample (true) or not (false).*/
    bool isOnBoundary;

    /// Initial (reference) area and 1/(current area) for this triangle.
    double initArea;
    double currAreaInv;

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
    Eigen::Matrix<unsigned int, 3, 1> nonVertexPatchNodesLabels;

    /* Position of triangle's centroid in initial (reference) x-y plane.*/
    Eigen::Vector3d refCentroid;
    Eigen::Vector3d centroid;

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
    Eigen::Matrix<double, 2, 2> programmedMetInv;

    // Determinant of the above dialledInvProgMetric matrix.
    double programmedMetInvDet;

    /* Dialled in programmed scalar 'tau' factor.*/
    double dialledProgTau;

    std::vector<Eigen::Vector3d> programmed_metric_infos = {{1, 1, 1}};
    std::vector<Eigen::Matrix<double, 2, 2>> programmed_metric_inv {Eigen::Matrix<double, 2, 2>::Identity()};
    std::vector<Eigen::Matrix<double, 2, 2>> programmed_second_fundamental_form = {Eigen::Matrix<double, 2, 2>::Zero()};
    std::vector<double> programmed_taus {1};

    /* Matrix representing the components of the *energetically* favoured
    (programmed) Second Fundamental Form in the x-y cartesian coordinate system
    of the initial flat state. */
    Eigen::Matrix<double, 2, 2> programmedSecFF;

    /* Matrix representing the total deformation gradient of this triangle. This
    maps the in-plane reference (initial) state (for which no z components are
    stored) to a triangle in 3D space, so it is 3x2.*/
    Eigen::Matrix<double, 3, 2> defGradient;

    /* The (1st order) approximation for the metric for this triangle, which
    equals defGradient.transpose) * defGradient.*/
    Eigen::Matrix<double, 2, 2> met;

    /* Inverse of metric.*/
    Eigen::Matrix<double, 2, 2> metInv;

    /* Determinant of the inverse of the metric.*/
    double metInvDet;

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

private:
    void updateProgrammedMetricExplicit(int stage_counter, double dial_in_factor);

    void updateProgrammedMetricImplicit(int stage_counter, double dial_in_factor);

    void updateProgrammedSecondFundamentalForm(int stage_counter, double dial_in_factor_root);

    void updateProgrammedTaus(int stage_counter, double dial_in_factor);

    void updateProgrammedMetricImplicit(double dirAngle, double lambda, double nu);

    void updateProgrammedMetricDynamically(int stage_counter, double dial_in_factor, double transfer_coefficient,
                                                     double min_height, double max_height);

public:
    double bendEnergyDensity;
    double stretchEnergyDensity;

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
        label = INT_MAX;
        isOnBoundary = false;
        vertexLabels.fill(INT_MAX);
        nonVertexPatchNodesLabels.fill(INT_MAX);
        edgeLabels.fill(INT_MAX);
        initArea = DBL_MAX;
        currAreaInv = DBL_MAX;
        currSides.fill(DBL_MAX);
        faceNormal.fill(DBL_MAX);
        initOutwardSideNormals.fill(DBL_MAX);
        invInitSidesMat.fill(DBL_MAX);
        programmedMetInv.fill(DBL_MAX);
        programmedMetInvDet = DBL_MAX;
        dialledProgTau = DBL_MAX;
        programmedSecFF.fill(DBL_MAX);
        defGradient.fill(DBL_MAX);
        met.fill(DBL_MAX);
        metInvDet = DBL_MAX;
        metInv.fill(DBL_MAX);
        halfPK1Stress.fill(DBL_MAX);
        patchSecDerivs.fill(DBL_MAX);
        secFF.fill(DBL_MAX);
        energyDensityDerivWRTSecFF.fill(DBL_MAX);
        bendEnergyDensityDerivWRTMetric.fill(DBL_MAX);
        matForPatchSecDerivs.fill(DBL_MAX);
        bendEnergyDensity = DBL_MAX;
        stretchEnergyDensity = DBL_MAX;
    }

    Triangle(int label, int id_0, int id_1, int id_2, const std::vector<Node> &nodes);

    // Declare other member functions.

    // Debugging function to display all member data.
    void display();

    void updateHalfPK1Stress(double stretchingPrefactor);

    Eigen::Matrix<double, 3, 3> getStretchingForces();

    Eigen::Matrix<double, 3, 3> getTriangleEdgeNormals();

    Eigen::Matrix<double, 3, 1> getBendingForce(const Eigen::Matrix<double, 3, 3> &normalDerivatives, int row);

    void updateMetric();

    void updateGeometricProperties();

    void updateSecondFundamentalForm(double bendingPreFac, double JPreFactor, double poissonRatio);

    void updateAngleDeficits(std::vector<double> &angleDeficits) const;

    void updateFirstFundamentalForm(double stretchingPreFac);

    double updateMatForPatchDerivs(const std::vector<Triangle> &triangles, const std::vector<Node> &nodes);

    /**
     * Returns the linear size of the triangle calculated as the shortest triangle altitude.
     * @return
     */
    double getLinearSize() const;

    double getHeight() const;

    void setLocalElongation(double local_elongation);

    void updateProgrammedQuantities(int stage_counter, double dial_in_factor, double dial_in_factor_root,
                                    bool is_lce_metric_used, bool is_elongation_dynamically_updated,
                                    double transfer_coefficient, double min_height, double max_height);

    void updateElasticForce(double bendingPreFac, double JPreFactor, double stretchingPreFac, double poisson_ratio);

    Eigen::Vector3d getNodeForce(unsigned int index) const;
};

#endif
