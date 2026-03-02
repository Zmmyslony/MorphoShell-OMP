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

This is the header file for the class that will correspond to each edge,
containing the labels of the nodes on the edge and at the other corners of the
two triangles sharing the edge, and the labels of those triangles.*/

#ifndef _EDGE_CLASS_TAG_
#define _EDGE_CLASS_TAG_

#include <Eigen/Dense>

//#include "CustomOutStreamClass.hpp"

class Edge {
public:

    /* Custom output stream allowing the debugging display function to print to
    a particular file in addition to std::cout.*/
//    CustomOutStreamClass edgeLogStream;

    // Label so this edge 'knows' which it is
    int label;

    /* Labels (and indexes in the nodes' container vector) of the nodes that the
    edge is defined to start and end at.*/
    std::pair<unsigned int, unsigned int> nodeLabels;

    /* Labels (and indices in the triangles' container vector) of the (either 1 or
    2) triangles that this edge is an edge of. We term these triangles 'adjacent'
    to the edge.*/
    std::vector<unsigned int> adjTriLabels;

    /*Constructor, taking a single argument which is an output file name
    that gets the debugging display function to print to a particular file, as
    well as to std::out. This should usually be the log file (as for logStream).
    I ensure that default data values are recognisable values,
    for debugging. */
    Edge() {
        label = INT_MAX;
        nodeLabels = {UINT_MAX, UINT_MAX};
    }

    Edge(unsigned int edge_label, unsigned int first_node_label, unsigned int second_node_label, unsigned int triangle_label);

    bool isBoundary();

    // Debugging function to display all member data.
    void display();

    bool operator==(const Edge &rhs) const;
};


struct Hash {public:
    std::size_t operator()(const Edge& edge) const {
        return UINT32_MAX * static_cast<std::size_t>(edge.nodeLabels.first) + static_cast<std::size_t>(edge.nodeLabels.second) ;

        return std::hash<unsigned int>()(edge.nodeLabels.first) ^ std::hash<unsigned int>()(edge.nodeLabels.second);
        // std::size_t h1 = std::hash<unsigned int>{}(edge.nodeLabels.first);
        // std::size_t h2 = std::hash<unsigned int>{}(edge.nodeLabels.first);
        // return h1 ^ (h2 << 1); // or use boost::hash_combine (see Discussion) https://en.cppreference.com/w/Talk:cpp/utility/hash
    }
};


#endif
