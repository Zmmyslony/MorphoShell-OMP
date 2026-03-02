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

Function to calculate the node neighbours to each node (those connected to it
by edges).
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#include <cstddef>
#include <vector>

#include "calcNodeNeighbours.hpp"
#include "../Node.hpp"
#include "../Edge.hpp"

void configureNodeAdjacency(std::vector<Node> &nodes, const std::vector<Edge> &edges) {
    /* Temporarily, we will store the node neighbour labels in a std::vector of
    std::vectors. At the end, the same information will then be stored as a
    dynamic-sized Eigen::VectorXd at each node. We use std::vector here because
    the push_back() function is helpful. There may well be a more efficient
    approach with less vector resizing, but this only occurs once and is
    therefore unlikely to be a bottleneck worth worrying about. We use
    Eigen::VectorXd eventually just because it i) provides an easy way to turn
    bounds checking on and off, and ii) already has overloaded functions set up
    for easy printing out to std::cout etc. */

    for (auto &edge: edges) {
        nodes[edge.nodeLabels.first].neighbourNodeLabels.emplace_back(edge.nodeLabels.second);
        nodes[edge.nodeLabels.second].neighbourNodeLabels.emplace_back(edge.nodeLabels.first);
    }
}
