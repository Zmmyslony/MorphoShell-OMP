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
#include <omp.h>

#include "calcSecFFsAndRelatedQuantities.hpp"
#include "../Node.hpp"
#include "../Triangle.hpp"


void calcSecFFsAndRelatedQuantities(
        std::vector<Triangle> &triangles,
        const SettingsStruct &settings) {

    // Temp objects used to calculate bending forces once 2nd F.F. is found

//    Eigen::Matrix<double, 6, 1> lsqDataVec;

    double bendingPreFac =
            settings.SheetThickness * settings.SheetThickness * settings.SheetThickness * settings.ShearModulus / 12.0;

    // secFF estimate
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        Eigen::Vector3d vectorOfSecFFComps = triangles[i].patchSecDerivs.transpose() * triangles[i].faceNormal;

        triangles[i].secFF(0, 0) = vectorOfSecFFComps(0);
        triangles[i].secFF(0, 1) = vectorOfSecFFComps(1);
        triangles[i].secFF(1, 0) = vectorOfSecFFComps(1);
        triangles[i].secFF(1, 1) = vectorOfSecFFComps(2);

        Eigen::Matrix<double, 2, 2> tempMat1 =
                triangles[i].dialledInvProgMetric * (triangles[i].secFF - triangles[i].dialledProgSecFF);
        Eigen::Matrix<double, 2, 2> tempMat2 = tempMat1 * triangles[i].dialledInvProgMetric;
        double tr_tempMat1 = tempMat1.trace();
        double tempScalar = bendingPreFac * triangles[i].detDialledInvProgMetric;

        double J = settings.GentFactor * tempScalar / (settings.SheetThickness * settings.SheetThickness);

        // Now calculate bending energy density for this triangles[i].
        double preGentBendEnergyDensity = tempScalar * ((tempMat1 * tempMat1).trace() + tr_tempMat1 * tr_tempMat1);
        double gentDerivFac = (1.0 + 2.0 * preGentBendEnergyDensity / J);

        /* Calculate the derivative of the bending energy density with respect
        to the secFF. */
        triangles[i].energyDensityDerivWRTSecFF =
                gentDerivFac * 2.0 * tempScalar * (tempMat2 + tr_tempMat1 * triangles[i].dialledInvProgMetric);
    }
}