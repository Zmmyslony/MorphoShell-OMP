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

Function to calculate, for each triangle, a list of the triangles it shares
edges with, and store this list as member data, as well as a list of edge labels
for the triangle. For a given triangle, the n elements of the edge-sharing
triangles list correspond to the first n elements of the edge labels list, with
the remaining edge labels corresponding to boundary edges.

The edges data structure is also set up, where each edge stores its two end
nodes, its adjacent triangles etc.
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)

#include <vector>
#include <unordered_set>
#include <algorithm> // For std::find
#include <limits> // For INT_MAX
#include <stdexcept>

#include "calcTriangleAdjacencies_And_Edges.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../Edge.hpp"



int calcTriangleAdjacencies_And_Edges(const std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                                      std::vector<Edge> &edges) {
    std::unordered_set<Edge, Hash> edges_set;
    for (int i = 0; i < triangles.size(); i++) {
        for (int j = 0; j < 3; j++) {
            Edge edge(edges_set.size(), triangles[i].vertexLabels[j], triangles[i].vertexLabels[(j + 1) % 3], i);
            auto is_inserted = edges_set.insert(edge);
            if (!is_inserted.second) {
                auto find = *edges_set.find(edge);
                edge.adjTriLabels.push_back(find.adjTriLabels[0]);
                edges_set.erase(find);
                edges_set.insert(edge);
            }
        }
    }
    edges.reserve(edges_set.size());
    int index = 0;
    std::vector<unsigned int> assigned_labels(triangles.size(), 0);

    for (const Edge &edge_const : edges_set) {
        Edge edge = edge_const;
        assert ((edge.adjTriLabels.size() <= 2, "Edge is shared by more than two triangles"));
        edge.label = index;
        edges.emplace_back(edge);
        int i = edge.adjTriLabels[0];
        triangles[i].edgeLabels[assigned_labels[i]] = index;
        assigned_labels[i]++;
        if (edge.isBoundary()) {
            triangles[i].isOnBoundary = true;
        } else {
            int j = edge.adjTriLabels[1];
            triangles[j].edgeLabels[assigned_labels[j]] = index;
            assigned_labels[j]++;

            triangles[i].edgeSharingTriLabels.push_back(j);
            triangles[j].edgeSharingTriLabels.push_back(i);
        }

        index++;
    }

#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        assert((assigned_labels[i] == 3, "Not all triangles have had all their edges assigned to them"));
    }

    return edges.size();
}
