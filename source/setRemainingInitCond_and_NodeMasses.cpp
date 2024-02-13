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

Function to finish setting up initial data for the flat LCE sheet: setting
node velocities to zero, and storing the initial in-plane sides' components for
the triangles (using the x-y plane basis), which will then not change. The
initial areas are also calculated and stored. Also, calculate node masses by
having each triangle contribute 1/3 of its initial mass to each of its vertcies.
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <Eigen/Dense>

#include "setRemainingInitCond_and_NodeMasses.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "Edge.hpp"
#include "Settings.hpp"

void setRemainingInitCond_and_NodeMasses(
        std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        std::vector<Edge> &edges,
        std::vector<std::vector<Eigen::Vector3d> > &programmedMetricInfos,
        std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedInvertedMetrics,
        std::vector<std::vector<double> > &programmedTaus,
        std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &programmedSecondFundamentalForms,
        const Settings &settings) {

    // Temp matrix used to hold initial triangle sides.
    Eigen::Matrix<double, 2, 2> initSidesMat;
    // Temp LU decomp of the above, used to check invertibility.
    Eigen::FullPivLU<Eigen::Matrix<double, 2, 2> > tempinitSidesMatDecomp;


    for (int i = 0; i < settings.num_nodes; ++i) {
        //Set all initial node velocities to zero
        nodes[i].vel.fill(0.0);
        //Set all nodes masses to zero before calculating them next
        nodes[i].mass = 0.0;
    }

    //Resize the vector to hold the first ('trivial') programmed tensors, that
    //is populated in the next loop.
    programmedMetricInfos[0].resize(settings.num_triangles);
    programmedInvertedMetrics[0].resize(settings.num_triangles);
    programmedTaus[0].resize(settings.num_triangles);
    programmedSecondFundamentalForms[0].resize(settings.num_triangles);

    for (int i = 0; i < settings.num_triangles; ++i) {
        /*set triangle sides' initial in-plane x-y basis components.*/
//        initSidesMat(0, 0) = triangles[i].currSides(0, 0);
//        initSidesMat(1, 0) = triangles[i].currSides(1, 0);
//        initSidesMat(0, 1) = triangles[i].currSides(0, 1);
//        initSidesMat(1, 1) = triangles[i].currSides(1, 1);
        initSidesMat = triangles[i].currSides.block<2, 2>(0, 0);

        //Store initial (reference) area
        triangles[i].initArea = 0.5 * initSidesMat.determinant();

        /*The reference (initial) state in-triangle-plane *outward* normals of
        each triangle's sides (with lengths = corresponding side lengths).
        initOutwardTriNormals.col(v) is opposite vertexLabels(v).*/
        triangles[i].initOutwardSideNormals.col(1) << -initSidesMat(1, 1), initSidesMat(0, 1);
        triangles[i].initOutwardSideNormals.col(2) << initSidesMat(1, 0), -initSidesMat(0, 0);
        triangles[i].initOutwardSideNormals.col(0) =
                -triangles[i].initOutwardSideNormals.col(1) - triangles[i].initOutwardSideNormals.col(2);

        // Now store inverse of initSidesMat permanently.
        tempinitSidesMatDecomp.compute(initSidesMat);
        if (!tempinitSidesMatDecomp.isInvertible()) {
            throw std::runtime_error(
                    "At least one triangle had a non-invertible initial sides matrix. \n"
                    "This should not occur in a reasonable mesh. Aborting.");
        } else {
            triangles[i].invInitSidesMat = initSidesMat.inverse();
        }

        /* Calculate vector storing, for each non-boundary edge, its initial
        length divided by sum of initial non-boundary edge lengths for the triangle.
        These could in future be used to weight the least squares fit, for instance.*/
        double initTotNonBoundaryLength = 0;
        triangles[i].initNonBoundEdgeLengthFracs.resize(triangles[i].edgeSharingTriLabels.size());
        for (int e = 0; e < 3; ++e) {
            if (!edges[triangles[i].edgeLabels(e)].isOnBoundary) {
                triangles[i].initNonBoundEdgeLengthFracs(e) = (
                        nodes[edges[triangles[i].edgeLabels(e)].nodeLabels(0)].pos -
                        nodes[edges[triangles[i].edgeLabels(e)].nodeLabels(1)].pos).norm();
                initTotNonBoundaryLength += triangles[i].initNonBoundEdgeLengthFracs(e);
            }
        }
        // Convert these initial edge lengths to fractions of the sum.
        triangles[i].initNonBoundEdgeLengthFracs = triangles[i].initNonBoundEdgeLengthFracs / initTotNonBoundaryLength;

        /* Add 1/3 of the mass of each triangle to each of its vertices.*/
        for (int v = 0; v < 3; ++v) {
            nodes[triangles[i].vertexLabels(v)].mass +=
                    settings.init_density * triangles[i].initArea * settings.sheet_thickness / 3.0;
        }

        /* Set first set of programmed tensors to be the trivial ones for the
        flat plane. This may be overridden later if
        settings.isDialingFromAnsatzEnabled == true. See main().*/
        if (settings.is_lce_mode_enabled) {
            programmedMetricInfos[0][i] << programmedMetricInfos[1][i](0), 1.0, programmedMetricInfos[1][i](
                    2);
            programmedTaus[0][i] = 1.0;
            programmedSecondFundamentalForms[0][i] = Eigen::Matrix<double, 2, 2>::Zero();
        } else {
            programmedInvertedMetrics[0][i] = Eigen::Matrix<double, 2, 2>::Identity();
            programmedTaus[0][i] = 1.0;
            programmedSecondFundamentalForms[0][i] = Eigen::Matrix<double, 2, 2>::Zero();
        }
    }
}




















