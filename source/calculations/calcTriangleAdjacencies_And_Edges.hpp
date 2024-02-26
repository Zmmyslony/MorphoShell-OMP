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

Header file for calcTriangleAdjacencies_And_Edges.cpp function.
*/

#include <vector>

#include "../Node.hpp"
#include "../Triangle.hpp"
#include "../Edge.hpp"
#include "../Settings.hpp"

/* Calculate triangle edge-sharing adjacencies (triangle member data), and
set up an 'edges' data structure to store further edge information. */
int calcTriangleAdjacencies_And_Edges(const std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                                      std::vector<Edge> &edges);
