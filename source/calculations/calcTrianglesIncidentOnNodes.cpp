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

Function to determine the labels of the triangles incident on each node (i.e.
having it as a vertex), and store these as node member data.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <cstddef>
#include "../../Eigen/Dense"
#include <vector>

#include "calcTrianglesIncidentOnNodes.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../SettingsStruct.hpp"

void calcTrianglesIncidentOnNodes(std::vector<Node> &nodes, const std::vector<Triangle> &triangles,
                                  const SettingsStruct &settings) {

    /* Temporarily, we will store the incident triangle labels for nodes in a
    std::vector of std::vectors. At the end, the same information will then be
    stored as a dynamic-sized Eigen::VectorXd at each node. We use std::vector
    here because the push_back() function is helpful. There may well be a
    more efficient approach with less vector resizing, but this only occurs once
    and is therefore unlikely to be a bottleneck worth worrying about. We use
    Eigen::VectorXd eventually just because it i) provides an easy way to turn
    bounds checking on and off, and ii) already has overloaded functions set up
    for easy printing out to std::cout etc. */

    std::vector<std::vector<int> > tempNodesIncidentTriLabels(settings.NumNodes);
    for (int i = 0; i < settings.NumTriangles; ++i) {
        //Add this triangle's label to each of its vertices in turn
        for (int v = 0; v < 3; ++v) {
            tempNodesIncidentTriLabels[triangles[i].vertexLabels(v)].push_back(i);
        }
    }


    /* Now we know how many triangles are incident on each node, we put the
    labels in the corresponding node member data Eigen::VectorXd.*/
//    omp_set_num_threads(8);
//#pragma omp parallel for
    for (int i = 0; i < settings.NumNodes; ++i) {

        nodes[i].incidentTriLabels.resize(tempNodesIncidentTriLabels[i].size());

        for (size_t j = 0; j < tempNodesIncidentTriLabels[i].size(); ++j) {
            nodes[i].incidentTriLabels(j) = tempNodesIncidentTriLabels[i][j];
        }
    }
}
