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
#include <algorithm> // For std::find
#include <limits> // For INT_MAX
#include <stdexcept>

#include "calcTriangleAdjacencies_And_Edges.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../Edge.hpp"

int calcTriangleAdjacencies_And_Edges(const std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                                      std::vector<Edge> &edges) {

    /* Initially we store the adjacent triangle labels for each triangle in a
    std::vector of std::vectors. At the end, the same information will then be
    stored as a dynamic-sized Eigen::VectorXd at each triangle. We use std::vector
    here because the emplace_back() function is helpful here. There may well be a
    more efficient approach with less vector resizing, but this only occurs once
    and is therefore unlikely to be a bottleneck worth worrying about. We use
    Eigen::VectorXd eventually just because it i) provides an easy way to turn
    bounds checking on and off, and ii) already has overloaded functions set up
    for easy printing out to std::cout etc. */
    std::vector<std::vector<int> > tempAdjTriLabels(triangles.size());

    // We do the same thing for the triangles' edge labels' member data
    std::vector<std::vector<int> > tempTriEdgeLabels(triangles.size());

    /* We don't yet know how many edges there are, but it will be a lot of most of
    the time, so we use reserve suitable space using an upper bound on the
    number of edges. */
    edges.reserve(3 * triangles.size());
    int e = 0; //Index for edges std::vector container

    // Begin main loop over triangles.
    for (int i = 0; i < triangles.size(); ++i) {

        //Loop over vertices of this triangle
        for (int v = 0; v < 3; ++v) {

            int vertLabel = triangles[i].vertexLabels(v);

            /*Loop over triangles touching this vertex. For each one that isn't
            triangle i, see if it shares a second vertex with triangle i. If it
            does, they must share an edge, so add that triangle to i's list of
            adjacent (:= edge sharing) triangles.*/
            for (int t = 0; t < nodes[vertLabel].incidentTriLabels.size(); ++t) {

                int incidentTriLabel = nodes[vertLabel].incidentTriLabels(t);

                if (incidentTriLabel != i) {

                    for (int p = 0; p < 3; ++p) {
                        for (int q = 0; q < 3; ++q) {
                            if (p != v && triangles[incidentTriLabel].vertexLabels(q) == triangles[i].vertexLabels(p)) {

                                /*At this point we have found an edge-sharing
                                triangle. We only add it to triangle i's list if
                                it is not already in that list. */
                                if (std::find(tempAdjTriLabels[i].begin(), tempAdjTriLabels[i].end(),
                                              incidentTriLabel) == tempAdjTriLabels[i].end()) {

                                    tempAdjTriLabels[i].emplace_back(incidentTriLabel);

                                    /*We may have also found a new non-boundary
                                    edge, so also fill the next unfilled entry
                                    in 'edges', as long as this adjacency was
                                    new for the other triangle too, meaning no
                                    entry as been made in 'edges' for this edge
                                    yet. We also add the new edge's label to the
                                    edge label lists of both triangles.*/
                                    if (std::find(tempAdjTriLabels[incidentTriLabel].begin(),
                                                  tempAdjTriLabels[incidentTriLabel].end(), i) ==
                                        tempAdjTriLabels[incidentTriLabel].end()) {

                                        tempTriEdgeLabels[incidentTriLabel].emplace_back(e);
                                        tempTriEdgeLabels[i].emplace_back(e);

                                        edges.emplace_back();
                                        edges[e].label = e;
                                        edges[e].nodeLabels(0) = vertLabel;
                                        edges[e].nodeLabels(1) = triangles[i].vertexLabels(p);
                                        edges[e].adjTriLabels.resize(2);
                                        edges[e].adjTriLabels << i, incidentTriLabel;
                                        edges[e].isOnBoundary = false;
                                        ++e;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*We have filled entries in 'edges' for all the non-boundary edges; now we
    do the entries for all the boundary edges. Note we could maybe be less
    wasteful here, cutting the j loop off once all edges accounted for, but
    this whole operation only happens once and therefore life is too short to
    worry about such small optimisations. Also, readability is more important
    here anyway.*/
    for (int i = 0; i < triangles.size(); ++i) {

        // i.e. if(this triangle has some boundary edges)...
        if (tempTriEdgeLabels[i].size() < 3) {

            triangles[i].isOnBoundary = true;

            std::size_t numTriEdgesAlreadyAccountedFor = tempTriEdgeLabels[i].size();

            /*See whether each pair of vertices corresponds to one of the
            already-existing edges. If not, we create that edge. Each j value
            picks out a vertex pair using mod arithmetic. */
            for (int j = 0; j < 3; ++j) {

                /*These will correspond to each possible pair of 'i's vertices
                in turn. */
                int node1 = triangles[i].vertexLabels(j);
                int node2 = triangles[i].vertexLabels((j + 1) % 3);

                bool addNewEdge = true;

                for (std::size_t k = 0; k < numTriEdgesAlreadyAccountedFor; ++k) {
                    /*i.e. 'If this entry in edges does not correspond to this
                    pair of nodes (i.e. to the current j value)': */
                    if ((edges[tempTriEdgeLabels[i][k]].nodeLabels(0) == node1 &&
                         edges[tempTriEdgeLabels[i][k]].nodeLabels(1) == node2)
                        || (edges[tempTriEdgeLabels[i][k]].nodeLabels(1) == node1 &&
                            edges[tempTriEdgeLabels[i][k]].nodeLabels(0) == node2)) {

                        /*If this already-accounted-for edge (k) does
                        correspond to the node pair given by j, then don't
                        create a new edge entry for that node pair. */
                        addNewEdge = false;
                    }
                }

                /* If this node pair didn't match any pre-existing edges for
                this triangle, then it's a new boundary edge, and the relevant
                data is now set up for it. Note that this ensures that the
                boundary edge labels are always at the back of edgeLabels, with
                the non-boundary edges at the front.*/
                if (addNewEdge == true) {

                    tempTriEdgeLabels[i].emplace_back(e);

                    edges.emplace_back();
                    edges[e].label = e;
                    edges[e].nodeLabels(0) = node1;
                    edges[e].nodeLabels(1) = node2;
                    edges[e].adjTriLabels.resize(1);
                    edges[e].adjTriLabels << i;
                    edges[e].isOnBoundary = true;
                    ++e;
                }
            }
        }
    }

    /* Now we now put the temporary triangle adjacency and edge labels into
    permanent member data storage. We also sort the edge labels for each
    triangle so that:
    For a given triangle, the n elements of the edge-sharing
    triangles list correspond to the first n elements of the edge labels list,
    with the remaining edge labels corresponding to boundary edges.

    We also calculate the indicesIntoEdgeSharingTriLabelsOfNeighbours member
    data for triangles here - see Triangle.hpp for explanation. */
    for (int i = 0; i < triangles.size(); ++i) {

        triangles[i].edgeSharingTriLabels.resize(tempAdjTriLabels[i].size());
        triangles[i].indicesIntoEdgeSharingTriLabelsOfNeighbours.resize(tempAdjTriLabels[i].size());

        for (std::size_t j = 0; j < tempAdjTriLabels[i].size(); ++j) {
            triangles[i].edgeSharingTriLabels(j) = tempAdjTriLabels[i][j];
        }

        for (int j = 0; j < 3; ++j) {
            triangles[i].edgeLabels(j) = tempTriEdgeLabels[i][j];
        }

        /* Now sort triangles[i].edgeLabels and calculate
        indicesIntoEdgeSharingTriLabelsOfNeighbours.*/
        std::vector<int> newEdgeLabels;
        for (int t = 0; t < triangles[i].edgeSharingTriLabels.size(); ++t) {

            int neighbourTri = triangles[i].edgeSharingTriLabels(t);

            /* Find which edge corresponds to this neighbouring triangle, and
            put its label in the newEdgeLabels container.*/
            for (int ed = 0; ed < 3; ++ed) {

                for (int n = 0; n < edges[triangles[i].edgeLabels(ed)].adjTriLabels.size(); ++n) {

                    if (edges[triangles[i].edgeLabels(ed)].adjTriLabels(n) == neighbourTri) {

                        newEdgeLabels.emplace_back(triangles[i].edgeLabels(ed));
                    }
                }
            }

            /* Now calculate indicesIntoEdgeSharingTriLabelsOfNeighbours for
            triangle i.*/
            for (std::size_t u = 0; u < tempAdjTriLabels[neighbourTri].size(); ++u) {

                if (tempAdjTriLabels[neighbourTri][u] == i) {

                    // Catch error if overflow occurs converting size_t to int.
                    if (u > INT_MAX) {
                        throw std::overflow_error(
                                "std::size_t loop variable u overflowed upon conversion to int (i.e. u was larger than INT_MAX)");
                    }

                    triangles[i].indicesIntoEdgeSharingTriLabelsOfNeighbours(t) = static_cast<int>(u);
                }
            }
        }

        /* Now find the other edges, which are the boundary edges, and add these
        to the back of the newEdgeLabels container.*/
        for (int ed = 0; ed < 3; ++ed) {

            if (edges[triangles[i].edgeLabels(ed)].isOnBoundary) {

                newEdgeLabels.emplace_back(triangles[i].edgeLabels(ed));
            }
        }
        /* Now update triangles[i].edgeLabels to the new, sorted version.*/
        for (int ed = 0; ed < 3; ++ed) {
            triangles[i].edgeLabels(ed) = newEdgeLabels[ed];
        }
    }

    /* Now shrink capacity of std::vector edges to match its actual number of
    elements (edges) (i.e. just allow discard of extra unused storage). */
    edges.shrink_to_fit();


    /* Calculate edgeAdjTriLabelSelectors for each triangle. See
    Triangle.hpp for more info.*/
    for (int i = 0; i < triangles.size(); ++i) {
        triangles[i].edgeAdjTriLabelSelectors.resize(triangles[i].edgeSharingTriLabels.size());

        for (int ed = 0; ed < triangles[i].edgeAdjTriLabelSelectors.size(); ++ed) {

            if (edges[triangles[i].edgeLabels(ed)].adjTriLabels(0) == i) {
                triangles[i].edgeAdjTriLabelSelectors(ed) = 1.0;
            } else {
                triangles[i].edgeAdjTriLabelSelectors(ed) = -1.0;
            }
        }
    }

    /* Finally, do some checking of some simple things as a bug check.
    Note the third condition is not mathematically prohibited but would
    indicate a very concerning mesh where some triangles share no edges with
    other triangles. More checks would be possible here; some are better than
    none.*/
    for (int i = 0; i < triangles.size(); ++i) {
        try {
            if (edges.size() != e) { throw std::runtime_error("incorrect edge size"); }
            if (triangles[i].edgeLabels.size() != 3) {
                throw std::runtime_error(
                        "has " + std::to_string(triangles[i].edgeLabels.size()) + " edges instead of 3.");
            }
            if (triangles[i].edgeSharingTriLabels.size() < 1 || triangles[i].edgeSharingTriLabels.size() > 3) {
                std::string neighbouring_triangle_labels;
                for (auto &label : triangles[i].edgeSharingTriLabels) {
                    neighbouring_triangle_labels += std::to_string(label) + " ";
                }
                throw std::runtime_error(
                        "shares edges with " + std::to_string(triangles[i].edgeSharingTriLabels.size()) +
                        " other triangles: " + neighbouring_triangle_labels + ".");
            }
            if ((triangles[i].isOnBoundary && triangles[i].edgeSharingTriLabels.size() >= 3)
                || (!triangles[i].isOnBoundary &&
                    (edges[triangles[i].edgeLabels(0)].isOnBoundary ||
                     edges[triangles[i].edgeLabels(1)].isOnBoundary ||
                     edges[triangles[i].edgeLabels(2)].isOnBoundary))
                    ) {
                throw std::runtime_error("boundary inconsistency.");
            }
        }
        catch (std::runtime_error &err) {
            throw std::runtime_error("Error: T" + std::to_string(i) + " " + err.what());
        }
    }
    // If no problems have been spotted, store the total number of edges
    return e;
}
