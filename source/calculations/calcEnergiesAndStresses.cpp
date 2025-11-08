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

Function to calculate stretching and bending potential energy density for each
triangle, and then the corresponding energy, filling NumTriangles-long vectors
with these. Also does the same for the node kinetic energies. We also calculate
my strain tensor and a corresponding scalar strain measure in this function,
which is in a sense a kind of non-dimensionalised stretch energy! And also some
stress stuff.
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>
#include <vector>
#include <cmath>

#include "calcEnergiesAndStresses.hpp"
#include "../Triangle.hpp"
#include "../Node.hpp"

void updateFirstFundamentalForms(
        std::vector<Triangle> &triangles,
        const CoreConfig &core_config) {

    double stretchingPreFac = 0.5 * core_config.getThickness() * core_config.getShearModulus();
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].updateFirstFundamentalForm(stretchingPreFac);
    }
}

void updateSecondFundamentalForms(
        std::vector<Triangle> &triangles,
        const CoreConfig &core_config) {

    double bendingPreFac = 0.5 * pow(core_config.getThickness(), 3) * core_config.getShearModulus() / 12;
    double JPreFactor = core_config.getGentFactor() / pow(core_config.getThickness(), 2);

#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].updateSecondFundamentalForm(bendingPreFac, JPreFactor, core_config.getPoissonRatio());
    }
}


void calcEnergiesAndStresses(const std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                             std::vector<double> &stretchEnergies, std::vector<double> &bendEnergies,
                             std::vector<double> &kineticEnergies, std::vector<double> &strainMeasures,
                             std::vector<Eigen::Vector2d> &cauchyStressEigenvals,
                             std::vector<Eigen::Matrix<double, 3, 2> > &cauchyStressEigenvecs,
                             const CoreConfig &core_config) {
    double stretchingPreFac = 0.5 * core_config.getThickness() * core_config.getShearModulus();
    double bendingPreFac = 0.5 * pow(core_config.getThickness(), 3) * core_config.getShearModulus() / 12;
    double JPreFactor = core_config.getGentFactor() / pow(core_config.getThickness(), 2);

    // Loop over triangles and calculate potential energies and energy densities.
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); ++i) {
        Triangle *triangle = &triangles[i];
        triangle->updateFirstFundamentalForm(stretchingPreFac);
        triangle->updateSecondFundamentalForm(bendingPreFac, JPreFactor, core_config.getPoissonRatio());

        stretchEnergies[i] = triangle->initArea * triangle->stretchEnergyDensity;
        bendEnergies[i] = triangle->initArea * triangle->bendEnergyDensity;

        Eigen::Matrix<double, 2, 3> defGradPseudoInv =
                (triangle->defGradient.transpose() * triangle->defGradient).inverse() *
                (triangle->defGradient.transpose());

        Eigen::Matrix<double, 3, 3> myStrainTensor =
                0.5 * defGradPseudoInv.transpose() * (triangle->met - triangle->programmedMetInv.inverse()) *
                defGradPseudoInv;

        strainMeasures[i] = sqrt((myStrainTensor.transpose() * myStrainTensor).trace() / 2.0);

        Eigen::Matrix<double, 3, 3> cauchyStress =
                sqrt(triangle->metInvDet) * (2 * triangle->getHalfPK1Stress(stretchingPreFac)) *
                triangle->defGradient.transpose();

        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> eigenSolver(cauchyStress);

        cauchyStressEigenvals[i](0) = eigenSolver.eigenvalues()(1);
        cauchyStressEigenvecs[i].col(0) = eigenSolver.eigenvectors().col(1);
        cauchyStressEigenvals[i](1) = eigenSolver.eigenvalues()(2);
        cauchyStressEigenvecs[i].col(1) = eigenSolver.eigenvectors().col(2);
    }

#pragma omp parallel for
    for (int n = 0; n < nodes.size(); ++n) {
        kineticEnergies[n] = 0.5 * nodes[n].mass * nodes[n].vel.dot(nodes[n].vel);
    }
}


