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

#ifndef _TRIANGLE_CLASS_TAG_
#define _TRIANGLE_CLASS_TAG_

#include <vector>
#include <Eigen/Dense>
#include "Node.hpp"
// #include "Settings.hpp"


class Triangle {
public:
    const Eigen::Vector3d *corner_nodes_pos[3];
    Eigen::Vector3d centroid;
    const Eigen::Vector3d *patch_nodes_pos[3];
    double programmedMetInvDet;     // Determinant of the programmed metric inverse .

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

    Eigen::Vector3d programmed_metric_info = {1, 1, 1};
    Eigen::Vector3d next_programmed_metric_info = {1, 1, 1};

    Eigen::Matrix2d programmed_metric_inv = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix2d next_programmed_metric_inv = Eigen::Matrix<double, 2, 2>::Identity();

    Eigen::Matrix2d programmed_second_fundamental_form = Eigen::Matrix<double, 2, 2>::Identity();
    Eigen::Matrix2d next_programmed_second_fundamental_form = Eigen::Matrix<double, 2, 2>::Identity();

    double programmed_tau = 1;
    double next_programmed_tau = 1;

    Eigen::Matrix<double, 2, 2> programmedSecFF;

    /*Matrix that is pre-calculated and then used repeatedly in finding the
    components of the second fundamental form estimated for this triangle. */
    Eigen::Matrix<double, 6, 3> matForPatchSecDerivs;

    /* Determinant of the inverse of the metric.*/
    double metInvDet;
    double dialledProgTau; // Dialled in programmed scalar 'tau' factor.

    /* Estimated Second Fundamental Form matrix (secFF) of the deformed surface,
    defined (as with the deformation gradient) with respect to the 'material'
    coordinate system, i.e. the coordinate chart that used to be the (x,y)
    cartesians of the flat initial state, and then deformed with the sheet. This
    is a 2x2 symmetric matrix. */

    // Forces exerted by this triangle on nodes associated with it
    Eigen::Vector3d *node_triangle_force[6];
    // Magnetisation density expressed in the reference system coordinates.
    Eigen::Vector3d reference_magnetisation_density = {0, 0, 0};
    Eigen::Vector3d reference_node_positions[3];

    double local_elongation = 1;
    double local_magnitude = 1;
    double bendEnergyDensity;
    double stretchEnergyDensity;
    double initArea;
    double currAreaInv;

    Eigen::Vector3i vertexLabels;
    Eigen::Vector3i edgeLabels;
    uint32_t nonVertexPatchNodesLabels[3];
    int label;

    /* Labels (and indices in the triangles' container vector) of the other
    triangles that share an edge with this triangle.*/
    Eigen::VectorXi edgeSharingTriLabels;

    bool isOnBoundary;

private:

    void updateProgrammedMetricExplicit(int stage_counter, double dial_in_factor);

    void updateProgrammedMetricImplicit(int stage_counter, double dial_in_factor);

    void updateProgrammedSecondFundamentalForm(int stage_counter, double dial_in_factor_root);

    void updateProgrammedTaus(int stage_counter, double dial_in_factor);

    void updateProgrammedMetricImplicit(double dirAngle, double lambda, double nu);

    void updateProgrammedTensorsDynamically(int stage_counter, double dial_in_factor, double transfer_coefficient,
                                            double min_height, double max_height);

    Eigen::Matrix<double, 3, 1> getBendingForcePatch(int row, const Eigen::Vector3d& faceNormal, const Eigen::Matrix<double, 2, 2>& energyDensityDerivWRTSecFF) const;

    Eigen::Matrix<double, 3, 1> getBendingForceNode(const Eigen::Vector3d& normalDerivatives, int row, const Eigen::Vector3d& faceNormal, const Eigen::Matrix<double, 2, 2>
                                                    & energyDensityDerivWRTSecFF) const;

    Eigen::Matrix<double, 3, 3> getStretchingForces(double stretchingPrefactor, const Eigen::Matrix<double, 2, 2>& metInv, const Eigen::Matrix<double, 3, 2>&
                                                    deformationGradient) const;

    Eigen::Matrix<double, 3, 3> getTriangleEdgeNormals(const Eigen::Matrix<double, 3, 2>& current_sides, const Eigen::Vector3d& face_normal) const;

public:
    Eigen::Matrix<double, 3, 2> getCurrentSides() const;
    Eigen::Matrix<double, 3, 2> getHalfPK1Stress(double stretchingPrefactor, const Eigen::Matrix<double, 2, 2>& metInv, const Eigen::Matrix<double, 3, 2>&
                                                 deformationGradient) const;

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
        label = -1;
        isOnBoundary = false;
        vertexLabels.fill(-1);
        edgeLabels.fill(INT_MAX);
        initArea = DBL_MAX;
        currAreaInv = DBL_MAX;
        initOutwardSideNormals.fill(DBL_MAX);
        invInitSidesMat.fill(DBL_MAX);
        programmedMetInv.fill(DBL_MAX);
        programmedMetInvDet = DBL_MAX;
        dialledProgTau = DBL_MAX;
        programmedSecFF.fill(DBL_MAX);
        metInvDet = DBL_MAX;
        matForPatchSecDerivs.fill(DBL_MAX);
        bendEnergyDensity = DBL_MAX;
        stretchEnergyDensity = DBL_MAX;
    }

    Triangle(int label, int id_0, int id_1, int id_2, const std::vector<Node> &nodes);
    Eigen::Matrix<double, 3, 2> getDeformationGradient() const;

    std::stringstream display();
    Eigen::Matrix2d getSecondFundamentalForm() const ;
    Eigen::Matrix2d getMetric() const;
    void updateGeometricProperties(double bending_pre_factor, double j_pre_factor, double poisson_ratio, double stretching_prefactor);

    void updateAngleDeficits(std::vector<double> &angleDeficits) const;

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

    void setNodeForceAddress(unsigned int index, Eigen::Vector3d *address);

    void updateMagneticForce(const Eigen::Vector3d &magnetic_field);

    void setReferenceMagnetisationDensity(const Eigen::Vector3d &magnetisation);
};

#endif
