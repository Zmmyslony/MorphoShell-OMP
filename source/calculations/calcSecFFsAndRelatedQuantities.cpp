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

Function to loop over triangular elements, and assign an approximate second
fundamental form (secFF) to each triangles[i]. This is stored as triangle member
data, as is the derivative of the bending energy with respect to the secFF, which
will be used later to calculate bending forces.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>
#include <vector>

#include "calcSecFFsAndRelatedQuantities.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"


void calcSecFFsAndRelatedQuantities(
        std::vector<Triangle> &triangles,
        const Settings &settings) {

    // Temp objects used to calculate bending forces once 2nd F.F. is found


    double bendingPreFac =
            pow(settings.SheetThickness, 3) * settings.ShearModulus / 12.0;

    double JPreFactor = settings.GentFactor / (settings.SheetThickness * settings.SheetThickness);

    // secFF estimate
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].calculateSecondFundamentalForm(bendingPreFac, JPreFactor);
    }
}