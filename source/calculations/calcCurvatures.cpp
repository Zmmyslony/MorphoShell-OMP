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

Function to use the metric and second fundamental form for each triangle
to calculate corresponding Gauss and mean curvatures, filling
NumTriangles-long vectors with these. If specified in the settings file, an
angle deficit for each node is also calculated, which is a measure of
integrated Gauss and geodesic curvature (i.e. the LHS of Gauss-Bonnet).*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <Eigen/Dense>
#include <vector>
#include <cmath> // For sqrt and acos

#include "calcCurvatures.hpp"
#include "../Triangle.hpp"

void calcCurvatures(
        const std::vector<Node> &nodes,
        const std::vector<Triangle> &triangles,
        std::vector<double> &gaussCurvatures,
        std::vector<double> &meanCurvatures,
        std::vector<double> &angleDeficits,
        std::vector<double> &interiorNodeAngleDeficits,
        std::vector<double> &boundaryNodeAngleDeficits,
        const SettingsStruct &settings) {


    // Objects used to calculate angle deficits.
    const double PI = 3.14159265358979323846;

#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        /* Old, more naive and direct approach using usual shape operator. */
        /*
        shapeOperator = (triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse() * triangles[i].secFF;
        gaussCurvatures[i] = shapeOperator.determinant();
        meanCurvatures[i] = 0.5 * shapeOperator.trace();
        */


        /* Approach using 'symmetric shape operator', which should be more
        accurate and stable - see Wang, Clark, Jiao 2009 and Jiao, Zha 2008.
        It works by essentially factoring out and thereby ignoring the pure
        rotation part of the deformation, which has no effect on curvatures.
        This is accomplished via QR factorisation of the deformation gradient
        following the 2009 paper, though the 2008 one uses SVD instead.

        The sym. shape op., being symmetric has orthogonal eigendirections,
        so they clearly don't tell you anything about curvature directions - I
        don't think they have any useful meaning. But the eigenvalues are the
        principal curvatures, just like the usual shape operator.*/
        Eigen::HouseholderQR<Eigen::Matrix<double, 3, 2> > qrDecomp(3, 2);
        qrDecomp.compute(triangles[i].defGradient);

        /* Extracting the 'R' part relies on the special and under-documented
        way Eigen stores the QR-decomposed matrix, which I believe is the LAPACK
        way. Inconveniently, it returns R as 3x2 upper triangular, not 2x2.*/
        Eigen::Matrix<double, 3, 2> stretchPart_3x2 = qrDecomp.matrixQR().triangularView<Eigen::Upper>();
        Eigen::Matrix<double, 2, 2> invStretchPart = stretchPart_3x2.block<2, 2>(0, 0).inverse(); // If it's not invertible you'll notice elsewhere!

        Eigen::Matrix<double, 2, 2> symShapeOp = invStretchPart.transpose() * triangles[i].secFF * invStretchPart;

        /* In case you want principle curvatures. For principal directions too,
        you could follow eqn(8) in the 2009 paper, or maybe easier: construct
        Curvature tensor using the pseudo-inverse of F, which is inv(F^T F)*F^T.
        Then find the eigendirections and discard the one with the ~0 eigenvalue.
        */
        /*
        double sqrtDiscrim = sqrt( (symShapeOp(0,0) - symShapeOp(1,1)) * (symShapeOp(0,0) - symShapeOp(1,1)) + 4.0 * symShapeOp(0,1) * symShapeOp(0,1) );
        double princCurv1 = 0.5 * (symShapeOp(0,0) + symShapeOp(1,1) + sqrtDiscrim);
        double princCurv2 = 0.5 * (symShapeOp(0,0) + symShapeOp(1,1) - sqrtDiscrim);
        gaussCurvatures[i] = princCurv1 * princCurv2;
        meanCurvatures[i] = 0.5 * (princCurv1 + princCurv2);
        */

        gaussCurvatures[i] = symShapeOp.determinant();
        meanCurvatures[i] = 0.5 * symShapeOp.trace();
    }


    // Calculate angle deficits at nodes, if specified in settings file.
    if (settings.isAngleDeficitsPrinted) {
#pragma omp parallel for
        for (int n = 0; n < settings.NumNodes; ++n) {
            if (!nodes[n].isOnBoundary) {
                angleDeficits.at(n) = 2.0 * PI;
            } else {
                angleDeficits.at(n) = PI;
            }
        }

#pragma omp parallel for
        for (int i = 0; i < triangles.size(); i++) {
            Eigen::Vector3d side2 = triangles[i].currSides.col(1) - triangles[i].currSides.col(0);

            double sideLength0 = triangles[i].currSides.col(0).norm();
            double sideLength1 = triangles[i].currSides.col(1).norm();
            double sideLength2 = side2.norm();

            angleDeficits[triangles[i].vertexLabels(0)] -= acos(
                    (triangles[i].currSides.col(0).dot(triangles[i].currSides.col(1))) / (sideLength0 * sideLength1));

            angleDeficits[triangles[i].vertexLabels(1)] -= acos(
                    -(triangles[i].currSides.col(0).dot(side2)) / (sideLength0 * sideLength2));

            angleDeficits[triangles[i].vertexLabels(2)] -= acos(
                    (triangles[i].currSides.col(1).dot(side2)) / (sideLength1 * sideLength2));
        }

        /* Also fill vectors holding interior and exterior angle deficits
        separately.*/
        int idxInto_interiorNodeAngleDeficits = 0;
        int idxInto_boundaryNodeAngleDeficits = 0;
        for (int n = 0; n < settings.NumNodes; ++n) {

            if (!nodes[n].isOnBoundary) {
                interiorNodeAngleDeficits.at(idxInto_interiorNodeAngleDeficits) = angleDeficits[n];
                idxInto_interiorNodeAngleDeficits += 1;
            } else {
                boundaryNodeAngleDeficits.at(idxInto_boundaryNodeAngleDeficits) = angleDeficits[n];
                idxInto_boundaryNodeAngleDeficits += 1;
            }
        }

        // OLD APPROACH, NOT LIKELY TO BE EVER USED AGAIN:
        /* A good definition of angular deficit for a boundary node is not clear.
        Just using Pi instead of 2*Pi is not great because that does not give
        zero on a non-straight boundary. Instead, we just assign boundary nodes
        an angle deficit that is the average of the deficits of neighbouring
        non-boundary nodes. If the only neighbours are boundary nodes (which
        should be rare), we don't assign an angle deficit.*/
        /*
        for(int n = 0; n < settings.NumNodes; ++n){

            if( nodes[n].isOnBoundary == true ){

                // 'Reset' angle deficit to zero.
                angleDeficits[n] = 0.0;

                int numNonBoundaryNeighbours = 0;

                for(int nn = 0; nn < nodes[n].neighbourNodeLabels.size(); ++nn){
                    if( nodes[nodes[n].neighbourNodeLabels(nn)].isOnBoundary == false ){
                        angleDeficits[n] += angleDeficits[nodes[n].neighbourNodeLabels(nn)];
                        numNonBoundaryNeighbours += 1;
                    }
                }
                angleDeficits[n] /= numNonBoundaryNeighbours;

            }
        }
        */ //remoooooooooooovvvvvvvveeeeeeeeeeeeee??????????
    }
}
