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

#include "CustomOutStreamClass.hpp"

class Node {
public:

    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
    CustomOutStreamClass nodeLogStream;

    // Label of this node, (also its index in the nodes' container vector).
    int label;

    /* Labels of the triangles with this node as a vertex ('incident triangles').
    Ordering is arbitrary.*/
    Eigen::VectorXi incidentTriLabels;

    /* Labels of this node's neighbours, i.e. those nodes connected to this node
    by triangle edges. Ordering is arbitrary.*/
    Eigen::VectorXi neighbourNodeLabels;

    /* Boolean representing whether the node is on a boundary of the sample
    (true) or not (false).*/
    bool isOnBoundary;

    /* Boolean representing whether the node's position is to be held clamped
    (true) or not (false).*/
    bool isClamped;

    // Bool indicating whether to impose a Seide displacement to this node or not.
    bool isSeideDisplacementEnabled;

    /* Bool representing whether the load force should be applied to this node
    (true) or not (false). For more general loading (multiple, spatially varying,
    non-trivial directions etc.) this will need generalising a lot.
    It would  require a major think about how to even mathematically represent
    general loading, and then how to put it in data format, code etc.*/
    bool isLoadForceEnabled;

    // Position vector (x, y and z coordinates).
    Eigen::Vector3d pos;

    // Velocity vector.
    Eigen::Vector3d vel;

    // Force vector.
    Eigen::Vector3d force;

    // Mass assigned to node.
    double mass;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for logStream).
    I ensure that default data values are recognisable values,
    for debugging.
    incidentTriLabels and neighbourNodeLabels are left with zero size at
    initialisation. */
    Node() {
        label = -1234;
        isSeideDisplacementEnabled = false;
        isOnBoundary = false;
        isClamped = false;
        isLoadForceEnabled = false;
        pos.fill(-123456);
        vel.fill(654321);
        force.fill(56789);
        mass = -56789;
    }

    // Declare other member functions.

    // Debugging function to display all member data.
    void display();

};

#endif
