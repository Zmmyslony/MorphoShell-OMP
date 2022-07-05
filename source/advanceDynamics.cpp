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
#include "SettingsStruct.hpp"
#include "CustomOutStreamClass.hpp"


bool isForceThresholdExceeded(const Node &node, const SettingsStruct &settings) {
    return node.force.norm() >= 100000.0 * settings.charForceScale;
}


void logForceThresholdExceeded(Node &node, std::vector<Triangle> &triangles, SettingsStruct &settings,
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


void advanceDynamics(std::vector<Node> &nodes, std::vector<Triangle> &triangles, SettingsStruct &settings,
                     CustomOutStreamClass &logStream) {

    for (int i = 0; i < settings.NumNodes; ++i) {
        if (isForceThresholdExceeded(nodes[i], settings)) {
            logForceThresholdExceeded(nodes[i], triangles, settings, logStream);
        }
        /*  Check the force is well-behaved. If not, throw error, and
            display offending nodes and its incident triangles.
            Error is caught in main(). */


        /* Set velocities based on either Gradient Descent (over dampened) dynamics
        or Newtonian dynamics. Then advance positions accordingly. */
        if (settings.isGradientDescentDynamicsEnabled) {
            // Gradient Descent dynamics.
            nodes[i].vel = settings.InitDensity * nodes[i].force / (settings.NumDampFactor * nodes[i].mass);
        } else {
            // Newtonian Dynamics.
            /* Advance velocity and *then* position (Semi-Implicit Euler, also
            called Symplectic Euler).*/
            nodes[i].vel += (settings.TimeStep / nodes[i].mass) * nodes[i].force;
        }

        // Advance position.
        nodes[i].pos += settings.TimeStep * nodes[i].vel;
    }

    if (settings.isControlledForceEnabled) {
        /*
        double gravAccel = settings.ApproxMinInitElemSize / (settings.TimeStep * settings.TimeStep);// Set g to be related to a characteristic acceleration in simulation.
        double upperSlideMass = settings.upperSlideWeight / gravAccel;
        // Advance upper slide too, displacement and velocity taken downwards.
        settings.upperSlideVel += (settings.TimeStep/upperSlideMass) * (settings.upperTotSlideForce + settings.upperSlideWeight);
        settings.upperSlideDisplacement += settings.TimeStep * settings.upperSlideVel;
        */


        // Advance upper slide too, displacement and velocity taken downwards.
        settings.upperSlideVel = (settings.upperTotSlideForce + settings.upperSlideWeight) / settings.slideDampingParam;

        // If doing constant weight experiment, need spacer to avoid squashing completely flat.
        if (settings.constSlideWeightFac > 0 && settings.upperSlideVel > 0 &&
            (settings.currSlideZCoord_upper - settings.initSlideZCoord_lower) <
            settings.SpacerHeight * settings.SampleCharLength) {
            settings.upperSlideVel = 0.0;
        }

        settings.upperSlideDisplacement += settings.TimeStep * settings.upperSlideVel;
    }
}
