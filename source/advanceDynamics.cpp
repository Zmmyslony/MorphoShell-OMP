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


bool isForceThresholdExceeded(const Node &node, const SettingsNew &settings) {
    return node.force.norm() >= 1e5 * settings.getForceScale();
}


void logForceThresholdExceeded(Node &node, std::vector<Triangle> &triangles, const SettingsNew &settings) {
    std::stringstream msg;
    msg << " ----------------------------------------" << std::endl;
    msg << " ------------CRASH REPORT----------------" << std::endl;
    msg << " ----------------------------------------" << std::endl;
    msg << "First offending node and its incident triangles: " << std::endl;

    msg << node.display().str();
    for (int t = 0; t < node.incidentTriLabels.size(); ++t) {
        msg << triangles[node.incidentTriLabels(t)].display().str();
    }
    std::cout << msg.str();
    throw std::runtime_error(msg.str() + "Suspiciously high force at node " + std::to_string(node.label) + std::string(" (") +
                             std::to_string(node.force.norm() / 1e5 * settings.getForceScale()) + " of limit).");
}


void advanceDynamics(std::vector<Node> &nodes, std::vector<Triangle> &triangles, SettingsNew &settings) {
#pragma omp parallel for
    for (int i = 0; i < nodes.size(); ++i) {
        if (isForceThresholdExceeded(nodes[i], settings)) {
            logForceThresholdExceeded(nodes[i], triangles, settings);
        }
        /*  Check the force is well-behaved. If not, throw error, and
            display offending nodes and its incident triangles.
            Error is caught in main(). */


        /* Set velocities based on either Gradient Descent (over dampened) dynamics
        or Newtonian dynamics. Then advance positions accordingly. */
        if (settings.getCore().isGradientDescentDynamics()) {
            // Gradient Descent dynamics.
            nodes[i].vel =
                    settings.getCore().getDensity() * nodes[i].force / (settings.getDampingFactor() * nodes[i].mass);
        } else {
            // Newtonian Dynamics.
            /* Advance velocity and *then* position (Semi-Implicit Euler, also
            called Symplectic Euler).*/
            nodes[i].vel += (settings.getTimeStepSize() / nodes[i].mass) * nodes[i].force;
        }

        // Advance position.
        nodes[i].pos += settings.getTimeStepSize() * nodes[i].vel;
    }
}
