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

Function to include damping and other non-deformation contributions to forces,
so the total force on each node is accounted for. Then BCs (e.g. clamping) are
accounted for.*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include "../../Eigen/Dense"
#include <vector>
#include <cmath>
#include <omp.h>

#include "calcNonDeformationForces_and_ImposeBCS.hpp"
#include "../Node.hpp"
//#include "../SettingsStruct.hpp"

std::pair<double, double> calcNonDeformationForces_and_ImposeBCS(std::vector<Node> &nodes, const double &time,
                                                                 const SettingsStruct &settings) {

    double totUpperSlideForce = 0;
    double totLowerSlideForce = 0;

    omp_set_num_threads(8);
#pragma omp parallel for reduction (+: totUpperSlideForce, totLowerSlideForce)
    for (int i = 0; i < nodes.size(); i++) {
        if (!settings.isGradientDescentDynamicsEnabled) {
            nodes[i].force += -settings.NumDampFactor * nodes[i].mass * nodes[i].vel / settings.InitDensity;
        }

        /* Perturbing 'prod' force, to prompt the sheet to buckle in the upward
        direction, and ensure evolution actally begins. The particular shape
        has no special justification, except that many target shapes have a
        bulge or protrusion vaguely in the middle, e.g. a nose. The shear
        modulus factor ensures the force has the correct dimensions and an order
        of magnitude roughly that of typical forces in the simulation. Note, an
        ansatz is probably a better way to choose a buckling direction than a
        prod.*/
        if (time < settings.ProdForceTime) {
            nodes[i].force(2) += -settings.ProdStrength * settings.ShearModulus * settings.SheetThickness *
                                 sqrt(nodes[i].pos(0) * nodes[i].pos(0) + nodes[i].pos(1) * nodes[i].pos(1));
        }

        // Simple load force.
        if (time < settings.LoadForceTime) {
            if (nodes[i].isLoadForceEnabled) {

                //nodes[i].force(2) += -settings.LoadStrength * settings.ShearModulus * settings.ApproxMinInitElemSize * settings.SheetThickness;

                double pullForce = settings.LoadStrength * settings.charForceScale * time / settings.BendingLongTime;
                if (nodes[i].pos(0) < 0) {
                    nodes[i].force(0) += -pullForce;
                    totUpperSlideForce += pullForce;
                } else {
                    nodes[i].force(0) += pullForce;
                    totLowerSlideForce += pullForce;
                }
            }
        }


        // FOR CONE SQUASHING/BUCKLING BETWEEN TWO SLIDES.
        // Force from "glass slides".
        double slideVertForce;
        if (settings.GlassCones == false) {
            // Use the if() statement with the && for single ridge experiment if you want to allow tip to go through the lower slide.
            //if( nodes[i].pos(2) < settings.initSlideZCoord_lower && nodes[i].pos(0)*nodes[i].pos(0) + nodes[i].pos(1)*nodes[i].pos(1) > (0.95*1.8)*(0.95*1.8) ){
            if (nodes[i].pos(2) < settings.initSlideZCoord_lower) {
                slideVertForce = settings.slideStiffnessPrefactor * settings.ShearModulus * settings.SheetThickness *
                                 (settings.initSlideZCoord_lower - nodes[i].pos(2));
                nodes[i].force(2) += slideVertForce;
                totLowerSlideForce += slideVertForce;

                // Friction
                double nonFrictionInPlaneForceSize = sqrt(
                        nodes[i].force(0) * nodes[i].force(0) + nodes[i].force(1) * nodes[i].force(1));
                if (nonFrictionInPlaneForceSize < settings.slideFrictionCoeff * fabs(slideVertForce)) {
                    nodes[i].force(0) = 0;
                    nodes[i].force(1) = 0;
                    nodes[i].vel(0) = 0;
                    nodes[i].vel(1) = 0;
                } else {
                    nodes[i].force(0) += -settings.slideFrictionCoeff * fabs(slideVertForce) * nodes[i].force(0) /
                                         nonFrictionInPlaneForceSize;
                    nodes[i].force(1) += -settings.slideFrictionCoeff * fabs(slideVertForce) * nodes[i].force(1) /
                                         nonFrictionInPlaneForceSize;

                }

            }
            // Use the if() statement with the && for single ridge experiment if you want to apply point(ish) load to tip.
            if (nodes[i].pos(2) > settings.currSlideZCoord_upper) {
                //if( nodes[i].pos(2) > settings.currSlideZCoord_upper && nodes[i].pos(0)*nodes[i].pos(0) + nodes[i].pos(1)*nodes[i].pos(1) < 2.0*settings.SheetThickness*2.0*settings.SheetThickness ){
                slideVertForce = settings.slideStiffnessPrefactor * settings.ShearModulus * settings.SheetThickness *
                                 (settings.currSlideZCoord_upper - nodes[i].pos(2));
                nodes[i].force(2) += slideVertForce;
                totUpperSlideForce += slideVertForce;

                // Friction
                double nonFrictionInPlaneForceSize = sqrt(
                        nodes[i].force(0) * nodes[i].force(0) + nodes[i].force(1) * nodes[i].force(1));
                if (nonFrictionInPlaneForceSize < settings.slideFrictionCoeff * fabs(slideVertForce)) {
                    nodes[i].force(0) = 0;
                    nodes[i].force(1) = 0;
                    nodes[i].vel(0) = 0;
                    nodes[i].vel(1) = 0;
                } else {
                    nodes[i].force(0) += -settings.slideFrictionCoeff * fabs(slideVertForce) * nodes[i].force(0) /
                                         nonFrictionInPlaneForceSize;
                    nodes[i].force(1) += -settings.slideFrictionCoeff * fabs(slideVertForce) * nodes[i].force(1) /
                                         nonFrictionInPlaneForceSize;

                }
            }
        }
            // GLASS CONES
        else {
            // Cylindrical polar coords.
            double r = sqrt(nodes[i].pos(0) * nodes[i].pos(0) + nodes[i].pos(1) * nodes[i].pos(1));
            double polarAng = atan2(nodes[i].pos(1), nodes[i].pos(0));

            // UPPER GLASS CONE
            // Distance this node is from the glass surface. If positive the node is
            // inside the glass.
            double distFromGlass = (nodes[i].pos(2) - r * tan(settings.ConeAngle) - settings.currSlideZCoord_upper) *
                                   cos(settings.ConeAngle);

            if (distFromGlass > 0) {
                double slideForceSize =
                        settings.slideStiffnessPrefactor * settings.ShearModulus * settings.SheetThickness *
                        distFromGlass;

                nodes[i].force(0) += slideForceSize * sin(settings.ConeAngle) * cos(polarAng);
                nodes[i].force(1) += slideForceSize * sin(settings.ConeAngle) * sin(polarAng);
                nodes[i].force(2) += -slideForceSize * cos(settings.ConeAngle);
                totUpperSlideForce += -slideForceSize * cos(settings.ConeAngle);
            }

            // LOWER GLASS CONE
            // Distance this node is from the glass surface. If negative the node is
            // inside the glass.
            distFromGlass = (nodes[i].pos(2) - r * tan(settings.ConeAngle) - settings.initSlideZCoord_lower) *
                            cos(settings.ConeAngle);

            if (distFromGlass < 0) {
                double slideForceSize =
                        -settings.slideStiffnessPrefactor * settings.ShearModulus * settings.SheetThickness *
                        distFromGlass;

                nodes[i].force(0) += -slideForceSize * sin(settings.ConeAngle) * cos(polarAng);
                nodes[i].force(1) += -slideForceSize * sin(settings.ConeAngle) * sin(polarAng);
                nodes[i].force(2) += slideForceSize * cos(settings.ConeAngle);
                totLowerSlideForce += slideForceSize * cos(settings.ConeAngle);
            }


            /*
            // SEIDE'S HORIZONTAL CLAMP ON ALL BOUNDARY NODES.
            if( nodes[i].isOnBoundary == true ){
                nodes[i].force(0) = 0.0;
                nodes[i].force(1) = 0.0;
            }
            */

            // For all nodes within a certain distance,
            // of the top or bottom in terms of vertical distance, we here kill any force component
            // not tangent to the simple expected base state (Seide) of a perfect cone with only
            // a single in-plane stress component.
            if (nodes[i].isOnBoundary) {
                Eigen::Vector3d theoryNormalVec;
                theoryNormalVec << cos(settings.ConeAngle) * cos(polarAng), cos(settings.ConeAngle) *
                                                                            sin(polarAng), sin(settings.ConeAngle);
                nodes[i].force = nodes[i].force - (nodes[i].force.dot(theoryNormalVec)) * theoryNormalVec;
            }

        }

        // BOUNDARY CONDITIONS
        if (nodes[i].isClamped) {
            // Full clamp.
            //nodes[i].force << 0.0,0.0,0.0;
            // Only clamp in z direction.
            //nodes[i].force(2) = 0.0;
            // Only clamp in z and x directions.
            nodes[i].force(0) = 0.0;
            nodes[i].force(1) = 0.0;
            nodes[i].force(2) = 0.0;
        }
    }
    return {totUpperSlideForce, totLowerSlideForce};
}


