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

Function to slightly perturb the node coordinates away from the positions read
in as initial data. This can be useful as a gentle way of allowing the
simulation to proceed if all nodes in the initial data lay in a plane for
example. Then, the random disturbance applied here will introduce small z
components in the node positions, allowing the simulation to 'break out' of the
plane. The size of the disturbance is chosen to be small relative to the
approximate smallest element size, to ensure the right kind of scale.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <vector>
#include <random>
#include <ctime>

#include "perturbInitialPositionsWithRandomNoise.hpp"
#include "../Node.hpp"

void perturbInitialPositionsWithRandomNoise(std::vector<Node> &nodes, double element_size) {

    /*Set random number generator . A simple and common one is chosen here:
    there is little point worrying about obtaining extremely 'good' random
    numbers for such a simple task, in which the quality of randomness is
    very unlikely to be important. The seed is chosen to be the current calendar
    time.*/
    std::minstd_rand aSimpleEngine(static_cast<long unsigned int>(std::time(nullptr)));

    /*Set distribution to be symmetric about zero, and extend to a small
    distance relative to the mesh spacing. The 0.1 factor is hard coded but
    could be changed if desired.*/
    std::uniform_real_distribution<double> distr(-element_size * 0.001,
                                                 element_size * 0.001);

    for (auto & node : nodes) {
        for (int c = 0; c < 3; ++c) {
            node.pos(c) += distr(aSimpleEngine);
        }
    }
}
