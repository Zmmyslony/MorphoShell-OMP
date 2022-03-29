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

Function to determine whether the system is in equilibrium or not (to within
chosen thresholds), based on:
1) The ratio of the largest node force to a typical simulation force scale given
by shear modulus * smallest initial element size * SheetThickness.
2) The largest value of the ratio of node velocity to a typical local velocity
scale - the (stretching) sound speed.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include "../Eigen/Dense"
#include <vector>
#include <cmath>

#include "equilCheck.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "StatusEnum.hpp"
#include "SettingsStruct.hpp"
#include "CustomOutStreamClass.hpp"

StatusEnum equilCheck(
        const std::vector<Node> &nodes,
        const std::vector<Triangle> &triangles,
        const SettingsStruct &settings,
        CustomOutStreamClass &logStream) {

    // Maximum (over nodes) non-damping force.
    double maxNonDampForce = 0.0;
    double tempNonDampForce;
    int maxNonDampForceNode = -1;

    // Max( node speed divided by local stretching sound (wave) speed ).
    double max_NodeSpeedOverLocalStretchingWaveSpeed = 0.0;
    int maxRelativeSpeedNode = -1;

    // Calculate node non-damping forces and speed ratios.
    for (int i = 0; i < settings.NumNodes; ++i) {

        if (!settings.isGradientDescentDynamicsEnabled) {
            tempNonDampForce = (nodes[i].force +
                                (settings.NumDampFactor * nodes[i].mass * nodes[i].vel / settings.InitDensity)).norm();
        } else {
            tempNonDampForce = nodes[i].force.norm();
        }

        if (maxNonDampForce < tempNonDampForce) {
            maxNonDampForce = tempNonDampForce;
            maxNonDampForceNode = i;
        }

        /* Calculate the local node speed : stretching wave speed ratio, based
        on the smallest (most stringent) progTau (and hence stiffness) of any of
        the triangles incident on this node. Note that for Gradient Descent
        dynamics, this quantity is essentially meaningless, because only the
        step size has meaning. You can see this in the fact that the 'timestep'
        and the velocity both have a damping factor in them, but it cancels out
        when they are combined to find the distance to move a node by. So only
        their combination is non-arbitrary in size.*/

        double smallestIncidentProgTau = triangles[nodes[i].incidentTriLabels(0)].dialledProgTau;

        for (int t = 0; t < nodes[i].incidentTriLabels.size(); ++t) {
            if (smallestIncidentProgTau > triangles[nodes[i].incidentTriLabels(t)].dialledProgTau) {
                smallestIncidentProgTau = triangles[nodes[i].incidentTriLabels(t)].dialledProgTau;
            }
        }

        double nodeSpeedOverLocalStretchingWaveSpeed =
                nodes[i].vel.norm() / sqrt(smallestIncidentProgTau * settings.ShearModulus / settings.InitDensity);

        if (max_NodeSpeedOverLocalStretchingWaveSpeed < nodeSpeedOverLocalStretchingWaveSpeed) {
            max_NodeSpeedOverLocalStretchingWaveSpeed = nodeSpeedOverLocalStretchingWaveSpeed;
            maxRelativeSpeedNode = i;
        }

    }

    /* Now check whether non-damping force and speed ratios are below chosen
    `equilibrium' thresholds.*/
    logStream.open();
    logStream << "Ratio of max non-damping force to characteristic force = "
              << maxNonDampForce / settings.charForceScale << " (node " << maxNonDampForceNode << ")" << "\n";
    logStream << "Max ratio of node speed to local stretching wave speed = "
              << max_NodeSpeedOverLocalStretchingWaveSpeed << " (node " << maxRelativeSpeedNode << ")" << "\n";

    if ((maxNonDampForce / settings.charForceScale) < settings.Equil_Force_To_CharForce_Ratio_Threshold
        && max_NodeSpeedOverLocalStretchingWaveSpeed < settings.Equil_Speed_To_SoundSpeed_Ratio_Threshold) {

        logStream << "Equilibrium reached" << std::endl;
        logStream.close();
        return equil_reached;
    } else {
        logStream << "Equilibrium not reached" << std::endl;
        logStream.close();
        return waiting_for_equil;
    }
}
