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
#include "configuration/core_config.h"

SimulationStatus
equilibriumCheck(const std::vector<Node> &nodes, const SettingsNew &settings, const std::vector<Triangle> &triangles,
                 CustomOutStreamClass &logStream) {

    std::vector<double> nodeNonDampingForce(nodes.size());
    std::vector<double> incidentProgTau(nodes.size());
    std::vector<double> velocity(nodes.size());

#pragma omp parallel for
    for (int i = 0; i < nodes.size(); ++i) {
        if (!settings.getCore().isGradientDescentDynamics()) {
            nodeNonDampingForce[i] = (nodes[i].force +
                                      (settings.getDampingFactor() * nodes[i].mass * nodes[i].vel /
                                       settings.getCore().getDensity())).norm();
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

    int max_non_damp_force_node =
            std::max_element(nodeNonDampingForce.begin(), nodeNonDampingForce.end()) - nodeNonDampingForce.begin();
    double max_non_damp_force = *std::max_element(nodeNonDampingForce.begin(), nodeNonDampingForce.end());

    double smallest_incident_prog_tau = *std::max_element(incidentProgTau.begin(), incidentProgTau.end());
    double stretching_wave_speed = sqrt(smallest_incident_prog_tau * settings.getCore().getShearModulus() /
                                        settings.getCore().getDensity());

    int max_relative_speed_node = std::max_element(velocity.begin(), velocity.end()) - velocity.begin();
    double max_relative_speed = *std::max_element(velocity.begin(), velocity.end()) / stretching_wave_speed;

    /* Now check whether non-damping force and speed ratios are below chosen
    `equilibrium' thresholds.*/

    double max_relative_force = max_non_damp_force / settings.getForceScale();
    logStream.open();

    SimulationStatus status;
    if (max_relative_force < settings.getForceScale()
        && max_relative_speed < settings.getCore().getEquilibriumSpeedScale()) {

        logStream << "\tEquilibrium reached." << std::endl;
        status = EquilibriumReached;
    } else {
        logStream << "\tEquilibrium not reached." << std::endl;
        status = WaitingForEquilibrium;
    }
    logStream << "\tRatio of max non-damping force to characteristic force = "
              << max_relative_force / settings.getForceScale()
              << " (node " << max_non_damp_force_node << ")." << std::endl
              << "\tRatio of max node speed to local stretching wave speed = "
              << max_relative_speed / settings.getCore().getEquilibriumSpeedScale()
              << " (node " << max_relative_speed_node << ")." << std::endl;
    logStream.close();
    return status;
}
