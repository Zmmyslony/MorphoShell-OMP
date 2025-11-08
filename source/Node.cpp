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

This file defines the member functions for the Node class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there.*/

#include <iostream>
#include <iomanip>

// #include "Settings.hpp"
#include "Node.hpp"
#include "configuration/gravity_config.h"

#define GRAVITY_MAGNITUDE (9.80665 * 1e3)

//This is a debugging tool to display the node's data
std::stringstream Node::display()
{
    std::stringstream msg;
    msg << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    msg << "Node " << label << ":" << std::endl;
    msg << "Labels of incident triangles: " << "\n" << incidentTriLabels << std::endl;
    msg << "neighbourNodeLabels = " << "\n" << neighbourNodeLabels << std::endl;
    msg << "Position = " << "\n" << pos << std::endl;
    msg << "Velocity = " << "\n" << vel << std::endl;
    msg << "Force = " << "\n" << force << std::endl;
    msg << "Mass = " << mass << std::endl;
    msg << "Boundary indicator: " << isOnBoundary << std::endl;
    msg << "Clamp indicator: " << (is_x_clamped || is_y_clamped || is_z_clamped) << std::endl;
    msg << "Load indicator: " << isLoadForceEnabled << std::endl;
    msg << "-----------------------------" << std::endl;
    return msg;
}

void Node::add_gravity(const GravityConfig& config)
{
    if (config.isGravityEnabled())
    {
        force(0) += config.getXGravityComponent() * mass * GRAVITY_MAGNITUDE;
        force(1) += config.getYGravityComponent() * mass * GRAVITY_MAGNITUDE;
        force(2) += config.getZGravityComponent() * mass * GRAVITY_MAGNITUDE;
    }
}


double Node::add_damping(const SettingsNew& settings_new)
{
    if (settings_new.getCore().isGradientDescentDynamics()) { return 0; }

    force += -settings_new.getDampingFactor() * mass * vel /
        settings_new.getCore().getDensity();
    return settings_new.getDampingFactor() * mass * pow(vel.norm(), 2) / settings_new.getCore().getDensity();
}


//void Node::add_prod_force(const Settings &settings) {
//    force(2) += -settings.prod_strength * settings.shear_modulus * settings.sheet_thickness *
//                sqrt(pos(0) * pos(0) + pos(1) * pos(1));
//}
//
//
//void Node::add_load_force(const Settings &settings, double time, double &upper_slide_force, double &lower_slide_force) {
//    if (isLoadForceEnabled) {
//        double pullForce = settings.load_strength * settings.char_force_scale * time / settings.bending_long_time;
//        if (pos(0) < 0) {
//            force(0) += -pullForce;
//            upper_slide_force += pullForce;
//        } else {
//            force(0) += pullForce;
//            lower_slide_force += pullForce;
//        }
//    }
//}
//
//void Node::add_slide_force(const Settings &settings, double height, bool is_bottom_slide, double &total_slide_force) {
//    bool is_interacting = is_bottom_slide && pos(2) < height || !is_bottom_slide && pos(2) > height;
//
//    if (is_interacting) {
//        double slice_vert_force =
//                settings.slide_stiffness_prefactor * settings.shear_modulus * settings.sheet_thickness *
//                (height - pos(2));
//        force(2) += slice_vert_force;
//        total_slide_force += slice_vert_force;
//
//        // Friction
//        double in_plane_force = sqrt(force(0) * force(0) + force(1) * force(1));
//        double friction_force = settings.slide_friction_coefficient * fabs(slice_vert_force);
//        if (in_plane_force < friction_force) {
//            force(0) = 0;
//            force(1) = 0;
//            vel(0) = 0;
//            vel(1) = 0;
//        } else {
//            force(0) -= friction_force * force(0) / in_plane_force;
//            force(1) -= friction_force * force(1) / in_plane_force;
//        }
//    }
//}
//
//void Node::add_cone_force(const Settings &settings, double tip_height, bool is_bottom_cone, double &total_cone_force) {
//    double r = sqrt(pos(0) * pos(0) + pos(1) * pos(1));
//    double polar_angle = atan2(pos(1), pos(0));
//
//    double distance_from_cone = (pos(2) - tip_height - r * tan(settings.cone_angle)) * cos(settings.cone_angle);
//    bool is_interacting = is_bottom_cone && distance_from_cone < 0 || !is_bottom_cone && distance_from_cone > 0;
//
//    if (is_interacting) {
//        double slide_force = settings.slide_stiffness_prefactor * settings.shear_modulus * settings.sheet_thickness *
//                             distance_from_cone;
//        if (is_bottom_cone) {
//            slide_force *= -1;
//        }
//        force(0) += slide_force * sin(settings.cone_angle) * cos(polar_angle);
//        force(1) += slide_force * sin(settings.cone_angle) * sin(polar_angle);
//        force(2) += -slide_force * cos(settings.cone_angle);
//        total_cone_force += -slide_force * cos(settings.cone_angle);
//    }
//}

void Node::apply_boundary_conditions()
{
    if (is_x_clamped)
    {
        force(0) = 0;
        vel(0) = 0;
    }
    if (is_y_clamped)
    {
        force(1) = 0;
        vel(1) = 0;
    }
    if (is_z_clamped)
    {
        force(2) = 0;
        vel(2) = 0;
    }
}

void Node::clamp(const CoreConfig& config)
{
    is_x_clamped = config.isXFixedBc();
    is_y_clamped = config.isYFixedBc();
    is_z_clamped = config.isZFixedBc();
}

void Node::addNodeForceAddress(Eigen::Vector3d* force_address)
{
    force_pointers.push_back(force_address);
}

void Node::updateForce()
{
    force = {0, 0, 0};
    for (auto& pt : force_pointers)
    {
        force.noalias() += *pt;
        *pt = {0, 0, 0};
    }
}

void Node::addForce()
{
    for (auto& pt : force_pointers)
    {
        force.noalias() += *pt;
        *pt = {0, 0, 0};
    }
}
