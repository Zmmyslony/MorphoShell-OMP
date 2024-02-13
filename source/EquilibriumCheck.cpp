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

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

#include "EquilibriumCheck.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "SimulationStatus.hpp"
#include "Settings.hpp"
#include "CustomOutStreamClass.hpp"

SimulationStatus equilibriumCheck(
        const std::vector<Node> &nodes,
        const std::vector<Triangle> &triangles,
        const Settings &settings,
        CustomOutStreamClass &logStream) {

    std::vector<double> nodeNonDampingForce(nodes.size());
    std::vector<double> incidentProgTau(nodes.size());
    std::vector<double> velocity(nodes.size());

#pragma omp parallel for
    for (int i = 0; i < settings.NumNodes; ++i) {
        if (!settings.isGradientDescentDynamicsEnabled) {
            nodeNonDampingForce[i] = (nodes[i].force +
                                (settings.NumDampFactor * nodes[i].mass * nodes[i].vel / settings.InitDensity)).norm();
        } else {
            nodeNonDampingForce[i] = nodes[i].force.norm();
        }

        /* Calculate the local node speed : stretching wave speed ratio, based
        on the smallest (most stringent) progTau (and hence stiffness) of
        the triangles incident on this node. Note that for Gradient Descent
        dynamics, this quantity is essentially meaningless, because only the
        step size has meaning. You can see this in the fact that the 'timestep'
        and the velocity both have a damping factor in them, but it cancels out
        when they are combined to find the distance to move a node by. So only
        their combination is non-arbitrary in size.*/

        incidentProgTau[i] = triangles[nodes[i].incidentTriLabels(0)].dialledProgTau;

        for (int t = 1; t < nodes[i].incidentTriLabels.size(); ++t) {
            double currentProgTau = triangles[nodes[i].incidentTriLabels(t)].dialledProgTau;
            if (incidentProgTau[i] > currentProgTau) {
                incidentProgTau[i] = currentProgTau;
            }
        }

        velocity[i] = nodes[i].vel.norm();
    }

    int maxNonDampForceNode = std::max_element(nodeNonDampingForce.begin(), nodeNonDampingForce.end()) - nodeNonDampingForce.begin();
    double maxNonDampForce = *std::max_element(nodeNonDampingForce.begin(), nodeNonDampingForce.end());

    double smallestIncidentProgTau = *std::max_element(incidentProgTau.begin(), incidentProgTau.end());
    double stretchingWaveSpeed = sqrt(smallestIncidentProgTau * settings.ShearModulus / settings.InitDensity);

    int maxRelativeSpeedNode = std::max_element(velocity.begin(), velocity.end()) - velocity.begin();
    double maxRelativeSpeed = *std::max_element(velocity.begin(), velocity.end()) / stretchingWaveSpeed;

    /* Now check whether non-damping force and speed ratios are below chosen
    `equilibrium' thresholds.*/
    logStream.open();
    logStream << "\tRatio of max non-damping force to characteristic force = "
              << maxNonDampForce / settings.charForceScale << " (node " << maxNonDampForceNode << ")." << std::endl;
    logStream << "\tMax ratio of node speed to local stretching wave speed = "
              << maxRelativeSpeed << " (node " << maxRelativeSpeedNode << ")." << std::endl;

    if ((maxNonDampForce / settings.charForceScale) < settings.Equil_Force_To_CharForce_Ratio_Threshold
        && maxRelativeSpeed < settings.Equil_Speed_To_SoundSpeed_Ratio_Threshold) {

        logStream << "Equilibrium reached." << std::endl;
        logStream.close();
        return equilibriumReached;
    } else {
        logStream << "Equilibrium not reached." << std::endl;
        logStream.close();
        return waitingForEquilibrium;
    }
}
