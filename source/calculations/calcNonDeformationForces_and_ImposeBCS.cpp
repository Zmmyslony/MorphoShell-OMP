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

//std::pair<double, double> calcNonDeformationForces_and_ImposeBCS(std::vector<Node> &nodes, const double &time,
//                                                                 const Settings &settings) {
//
//    double totUpperSlideForce = 0;
//    double totLowerSlideForce = 0;
//
//#pragma omp parallel for reduction (+: totUpperSlideForce, totLowerSlideForce)
//    for (int i = 0; i < nodes.size(); i++) {
//        if (!settings.is_gradient_descent_dynamics_enabled) {
//            nodes[i].add_damping(settings);
//        }
//        nodes[i].add_gravity(settings);
//
//        if (time < settings.prod_force_time) {
//            nodes[i].add_prod_force(settings);
//        }
//
//        // Simple load force.
//        if (time < settings.load_force_time) {
//            nodes[i].add_load_force(settings, time, totUpperSlideForce, totLowerSlideForce);
//        }
//
//        // If slide stiffness is positive, apply forces from the slides.
//        if (settings.slide_stiffness_prefactor > 0) {
//            if (!settings.glass_cones) {
//                nodes[i].add_slide_force(settings, settings.init_slide_z_coord_lower, true, totLowerSlideForce);
//                nodes[i].add_slide_force(settings, settings.init_slide_z_coord_upper, false, totUpperSlideForce);
//            }
//                // GLASS CONES
//            else {
//                nodes[i].add_cone_force(settings, settings.curr_slide_z_coord_upper, true, totLowerSlideForce);
//                nodes[i].add_cone_force(settings, settings.init_slide_z_coord_upper, false, totUpperSlideForce);
//
//                if (nodes[i].isOnBoundary) {
//                    double polar_angle = atan2(nodes[i].pos(1), nodes[i].pos(0));
//                    Eigen::Vector3d theoryNormalVec;
//                    theoryNormalVec << cos(settings.cone_angle) * cos(polar_angle),
//                            cos(settings.cone_angle) * sin(polar_angle), sin(settings.cone_angle);
//                    nodes[i].force -= (nodes[i].force.dot(theoryNormalVec)) * theoryNormalVec;
//                }
//            }
//        }
//
//        // BOUNDARY CONDITIONS
//        nodes[i].apply_boundary_conditions();
//    }
//    return {totUpperSlideForce, totLowerSlideForce};
//}


