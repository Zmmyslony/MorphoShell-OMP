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

Function to advance the dynamics of the system by evolving node
positions and velocities according to the total forces on them. The scheme used
is so-called 'semi-implicit Euler'. This is first order, but it is stable and
energy conserving, which is good enough, since we are not interested in
resolving the details of the dynamics, only the final equilibrium state. If
chosen in settings, Gradient Descent dynamics are used instead.

The forces are also checked, to catch code crashed in which the forces usually
'blow up'.*/

#include <iostream>

#include "advanceDynamics.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "Settings.hpp"
#include "CustomOutStreamClass.hpp"


bool isForceThresholdExceeded(const Node &node, const Settings &settings) {
    return node.force.norm() >= 1e5 * settings.char_force_scale;
}


void logForceThresholdExceeded(Node &node, std::vector<Triangle> &triangles, Settings &settings,
                               CustomOutStreamClass &logStream) {
    logStream.open();
    logStream << " ----------------------------------------" << std::endl;
    logStream << " ------------CRASH REPORT----------------" << std::endl;
    logStream << " ----------------------------------------" << std::endl;
    logStream << "Offending node and its incident triangles: " << std::endl;
    logStream.close();
    node.display();
    for (int t = 0; t < node.incidentTriLabels.size(); ++t) {
        triangles[node.incidentTriLabels(t)].display();
    }
    throw std::runtime_error("suspiciously_high_force");
}


void advanceDynamics(std::vector<Node> &nodes, std::vector<Triangle> &triangles, Settings &settings,
                     CustomOutStreamClass &logStream) {

    for (int i = 0; i < settings.num_nodes; ++i) {
        if (isForceThresholdExceeded(nodes[i], settings)) {
            logForceThresholdExceeded(nodes[i], triangles, settings, logStream);
        }
        /*  Check the force is well-behaved. If not, throw error, and
            display offending nodes and its incident triangles.
            Error is caught in main(). */


        /* Set velocities based on either Gradient Descent (over dampened) dynamics
        or Newtonian dynamics. Then advance positions accordingly. */
        if (settings.is_gradient_descent_dynamics_enabled) {
            // Gradient Descent dynamics.
            nodes[i].vel = settings.init_density * nodes[i].force / (settings.num_damp_factor * nodes[i].mass);
        } else {
            // Newtonian Dynamics.
            /* Advance velocity and *then* position (Semi-Implicit Euler, also
            called Symplectic Euler).*/
            nodes[i].vel += (settings.time_step / nodes[i].mass) * nodes[i].force;
        }

        // Advance position.
        nodes[i].pos += settings.time_step * nodes[i].vel;
    }

    if (settings.is_controlled_force_enabled) {
        /*
        double gravAccel = settings.ApproxMinInitElemSize / (settings.TimeStep * settings.TimeStep);// Set g to be related to a characteristic acceleration in simulation.
        double upperSlideMass = settings.upperSlideWeight / gravAccel;
        // Advance upper slide too, displacement and velocity taken downwards.
        settings.upperSlideVel += (settings.TimeStep/upperSlideMass) * (settings.upperTotSlideForce + settings.upperSlideWeight);
        settings.upperSlideDisplacement += settings.TimeStep * settings.upperSlideVel;
        */


        // Advance upper slide too, displacement and velocity taken downwards.
        settings.upper_slide_vel = (settings.upper_tot_slide_force + settings.upper_slide_weight) / settings.slide_damping_param;

        // If doing constant weight experiment, need spacer to avoid squashing completely flat.
        if (settings.const_slide_weight_fac > 0 && settings.upper_slide_vel > 0 &&
                (settings.curr_slide_z_coord_upper - settings.init_slide_z_coord_lower) <
                settings.spacer_height * settings.sample_char_length) {
            settings.upper_slide_vel = 0.0;
        }

        settings.upper_slide_displacement += settings.time_step * settings.upper_slide_vel;
    }
}
