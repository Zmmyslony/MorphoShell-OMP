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

Function to include damping and other non-deformation contributions to forces,
so the total force on each node is accounted for. Then BCs (e.g. clamping) are
accounted for.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>
#include <vector>
#include <cmath>

#include "calcNonDeformationForces_and_ImposeBCS.hpp"
#include "../Node.hpp"

std::pair<double, double> calcNonDeformationForces_and_ImposeBCS(std::vector<Node> &nodes, const double &time,
                                                                 const Settings &settings) {

    double totUpperSlideForce = 0;
    double totLowerSlideForce = 0;

#pragma omp parallel for reduction (+: totUpperSlideForce, totLowerSlideForce)
    for (int i = 0; i < nodes.size(); i++) {
        if (!settings.isGradientDescentDynamicsEnabled) {
            nodes[i].add_damping(settings);
        }
        nodes[i].add_gravity(settings);

        if (time < settings.ProdForceTime) {
            nodes[i].add_prod_force(settings);
        }

        // Simple load force.
        if (time < settings.LoadForceTime) {
            nodes[i].add_load_force(settings, time, totUpperSlideForce, totLowerSlideForce);
        }

        // FOR CONE SQUASHING/BUCKLING BETWEEN TWO SLIDES.
        // Force from "glass slides".
        if (!settings.GlassCones) {
            nodes[i].add_slide_force(settings, settings.initSlideZCoord_lower, true, totLowerSlideForce);
            nodes[i].add_slide_force(settings, settings.initSlideZCoord_upper, false, totUpperSlideForce);
        }
            // GLASS CONES
        else {
            nodes[i].add_cone_force(settings, settings.currSlideZCoord_upper, true, totLowerSlideForce);
            nodes[i].add_cone_force(settings, settings.initSlideZCoord_upper, false, totUpperSlideForce);

            /*
            // SEIDE'S HORIZONTAL CLAMP ON ALL BOUNDARY NODES.
            if( nodes[i].isOnBoundary == true ){
                nodes[i].force(0) = 0.0;
                nodes[i].force(1) = 0.0;
            }
            */

            // For all nodes within a certain distance,
            // of the top or bottom in terms of vertical distance, we here kill any force component
            // not tangent to the simple expected base state (Seide) of a perfect cone with only
            // a single in-plane stress component.
            if (nodes[i].isOnBoundary) {
                double polar_angle = atan2(nodes[i].pos(1), nodes[i].pos(0));
                Eigen::Vector3d theoryNormalVec;
                theoryNormalVec << cos(settings.ConeAngle) * cos(polar_angle),
                        cos(settings.ConeAngle) * sin(polar_angle), sin(settings.ConeAngle);
                nodes[i].force -= (nodes[i].force.dot(theoryNormalVec)) * theoryNormalVec;
            }
        }

        // BOUNDARY CONDITIONS
        nodes[i].apply_boundary_conditions();
    }
    return {totUpperSlideForce, totLowerSlideForce};
}


