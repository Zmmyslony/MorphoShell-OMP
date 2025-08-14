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

This is the header file for the class that will be placed at each node,
containing the node's position, velocity etc.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef _NODE_CLASS_TAG_
#define _NODE_CLASS_TAG_

#include <Eigen/Dense>
#include <cfloat>

// #include "Settings.hpp"
//#include "CustomOutStreamClass.hpp"
#include "configuration/gravity_config.h"
#include "settings_new.h"

class Node {
    /// Pointers to a field in a triangles which correspond to this node.
    std::vector<Eigen::Vector3d *> force_pointers;
public:
    void addNodeForceAddress(Eigen::Vector3d *force_address);
    void updateForce();
    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
//    CustomOutStreamClass nodeLogStream;

    // Label of this node, (also its index in the nodes' container vector).
    int label = INT_MAX;

    /* Labels of the triangles with this node as a vertex ('incident triangles').
    Ordering is arbitrary.*/
    Eigen::VectorXi incidentTriLabels;

    /* Labels of this node's neighbours, i.e. those nodes connected to this node
    by triangle edges. Ordering is arbitrary.*/
    Eigen::VectorXi neighbourNodeLabels;

    /* Boolean representing whether the node is on a boundary of the sample
    (true) or not (false).*/
    bool isOnBoundary = false;

    /* Boolean representing whether the node's position is to be held clamped
    (true) or not (false).*/
    bool is_x_clamped = false;
    bool is_y_clamped = false;
    bool is_z_clamped = false;

    // Bool indicating whether to impose a Seide displacement to this node or not.
    bool isSeideDisplacementEnabled = false;

    /* Bool representing whether the load force should be applied to this node
    (true) or not (false). For more general loading (multiple, spatially varying,
    non-trivial directions etc.) this will need generalising a lot.
    It would  require a major think about how to even mathematically represent
    general loading, and then how to put it in data format, code etc.*/
    bool isLoadForceEnabled = false;

    // Position vector (x, y and z coordinates).
    Eigen::Vector3d pos;
    Eigen::Vector3d prev_pos;

    // Velocity vector.
    Eigen::Vector3d vel;

    // Force vector.
    Eigen::Vector3d force;

    // Mass assigned to node.
    double mass = DBL_MAX;

    // Area calculated from the mass
    double area = DBL_MAX;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for logStream).
    I ensure that default data values are recognisable values,
    for debugging.
    incidentTriLabels and neighbourNodeLabels are left with zero size at
    initialisation. */
    Node() {
        pos.fill(DBL_MAX);
        vel.fill(DBL_MAX);
        force.fill(DBL_MAX);
        prev_pos.fill(DBL_MAX);
    }

    explicit Node(int n_label, const double positions[3]) {
        label = n_label;
        pos(0) = positions[0];
        pos(1) = positions[1];
        pos(2) = positions[2];
        prev_pos = pos;
        vel.fill(DBL_MAX);
        force.fill(DBL_MAX);
    }

    // Declare other member functions.

    // Debugging function to display all member data.
    std::stringstream display();

    /**
     * Adds gravitational force using the multiplier from the settings file. Use multiplier = 1 for normal pulling in
     * z-direction.
     * @param config
     */
    void add_gravity(const GravityConfig &config);

    /**
     * Adds damping force that is proportional to the velocity and numerical damping factor and returns the power loss:
     * F = - a * v * m / rho
     * @param settings_new
     * @return power loss
     */
    double add_damping(const SettingsNew &settings_new);

    /** Perturbing 'prod' force, to prompt the sheet to buckle in the upward
    *  direction, and ensure evolution actually begins. The particular shape
    *  has no special justification, except that many target shapes have a
    * bulge or protrusion vaguely in the middle, e.g. a nose. The shear
    *  modulus factor ensures the force has the correct dimensions and an order
    *  of magnitude roughly that of typical forces in the simulation. Note, an
    *  ansatz is probably a better way to choose a buckling direction than a
    *  prod.
    */
//    void add_prod_force(const Settings &settings);

    /**
    * Adds load force to the nodes that have isLoadForceEnabled = True.
    */
//    void add_load_force(const Settings &settings, double time, double &upper_slide_force, double &lower_slide_force);

    void apply_boundary_conditions();

    /**
     * Adds force coming from interaction with the "glass" slides. The slides can either approach from the top or the
     * bottom.
     */
//    void add_slide_force(const Settings &settings, double height, bool is_bottom_slide, double &total_slide_force);

    /**
     * Adds force coming from the interaction with the "glass" cone. The cones can either approach from the top or the
     * bottom and are centred as x = 0, y = 0.
     * @param settings
     * @param tip_height
     * @param is_bottom_cone
     * @param total_cone_force
     */
//    void add_cone_force(const Settings &settings, double tip_height, bool is_bottom_cone, double &total_cone_force);

    void clamp(const CoreConfig &config);
};

#endif
