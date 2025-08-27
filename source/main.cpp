/* Author: Daniel Duffy, University of Cambridge, dld34@cam.ac.uk

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
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.`
/////////////////////////////////////////////////////

Libconfig++ is distributed under the Lesser GPL license (2.1 or later), and copyright is held by Mark A Lindner (at least in large part).
Version 1.7.2 was used for this code, which can be accessed via: http://hyperrealm.github.io/libconfig/

Eigen is distributed under the Mozilla Public License v. 2.0 , and copyright is held by Gael Guennebaud and Benoit Jacob (at least in large part).
Version 3.3.7 was used for this code, which can be accessed via: http://eigen.tuxfamily.org/

//////////////////////////////////////////////////////////////////////////

Main file for Shellmorph, the main code developed for my PhD.

This code simulates a thin shape-morphing 2D sheet with a programmed metric and
second fundamental form.
The main object in the simulations is an unstructured triangulated mesh, with
data stored in three std:vectors: 'nodes', 'triangles' and 'edges', each element
of which holds a class, in turn holding the data for that node/triangle/edge.
So for example nodes[0] is a Node instance holding the data for the first
node, which is labelled zero, matching the indexing of 'nodes'. The initial node
and triangle data is read in from a VTK Legacy file, which is also the output
format. Paraview is recommended for visualising these files.

The code executable takes either two or three command line arguments. The first
is a settings file, and the second is an input data file containing the planar
2D initial (reference) state mesh and the programmed quantities. The optional
third argument is an `ansatz' file (also .vtk), which can be an output file from
a previous simulation. If this is provided, the nodes will be moved to the
positions specified in the ansatz file before dynamics begins from that state.
Thus you can terminate a simulation and restart it from where it left off with
a new settings file (though the velocities start from zero upon restarting).

An elastic energy that penalises stretching and bending deviations from the
programmed metric and secondfundamental form is minimised (using analytical
derivatives) via either linearly-damped Newtonian dynamics or gradient descent.
The former is on the whole recommended, while the latter can be useful for
smoothing out meshes that have significant angles between triangles i.e. small
wavelength oscillations. Equilibrium is deemed achieved when two dimensionless
numbers fall below thresholds specified in the settings file. These numbers are
maxima (over the mesh) of non-dimensionalised node speed and elastic force.
*/

// Turn Eigen bounds checking off for speed (after running with checks naturally).
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#define _USE_MATH_DEFINES

#include "simulation.h"

int main(int argc, char *argv[]) {
    Simulation simulation(argc, argv);
    return simulation.run_simulation();
}
