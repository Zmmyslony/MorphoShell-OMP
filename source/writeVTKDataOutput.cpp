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

Function to write output data to a text file. To put integers in filenames
we use std::to_string from <string>. This is a feature from c++11. If that is
not acceptable for some reason, one could use a std::stringstream object from
<sstream> followed by a .str(), which would give a std::string
containing the filename, which could then be passed as outFile(filename.c_str()).
The data format is Legacy VTK PolyData, to match the read-in input files.
See the explanation in readVTKData.cpp for the format etc.
The 'title' in the data file will contain the current stepcount, time and
lambda_Dynamic.
The directory the files are put in is created in main().
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <fstream>
#include <iomanip> //for setting output precision etc
#include <vector>
#include <stdexcept>

#include "writeVTKDataOutput.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "Settings.hpp"
#include "functions/kahanSum.hpp"

void writeVTKDataOutput(const std::vector<Node> &nodes, const std::vector<Triangle> &triangles, const int &stepcount,
                        const double &time, const double &currDialInFactor, const size_t &progTensorSequenceCounter,
                        const std::vector<double> &gaussCurvatures, const std::vector<double> &meanCurvatures,
                        const std::vector<double> &angleDeficits, const std::vector<double> &interiorNodeAngleDeficits,
                        const std::vector<double> &boundaryNodeAngleDeficits,
                        const std::vector<double> &stretchEnergyDensities,
                        const std::vector<double> &bendEnergyDensities, const std::vector<double> &stretchEnergies,
                        const std::vector<double> &bendEnergies, const std::vector<double> &kineticEnergies,
                        const std::vector<double> &strainMeasures,
                        const std::vector<Eigen::Vector2d> &cauchyStressEigenvals,
                        const std::vector<Eigen::Matrix<double, 3, 2> > &cauchyStressEigenvecs,
                        const SettingsNew &settings, const std::string &outputDirName) {

    std::ofstream outFile(outputDirName + "/stepcount_" + std::to_string(stepcount) + "_output.vtk");
//    std::stringstream outFile;
    if (!outFile) {
        throw std::runtime_error("Error: Problem creating or opening output file.");
    }


    outFile << std::scientific << std::setprecision(14)
            << "# vtk DataFile Version 4.2" << "\n"
            << "dial_in_factor = " << currDialInFactor << " dialling programmed tensors _"
            << progTensorSequenceCounter + 1 << ", time = " << time << ", stepcount = " << stepcount;

    /* If specified in settings file, also print total stretching and bending
    energies here.*/
    if (settings.getCore().isEnergyPrinted()) {
        double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.getStretchEnergyScale();
        double nonDimBendEnergy = kahanSum(bendEnergies) / settings.getStretchEnergyScale();
        double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.getStretchEnergyScale();
        outFile << ", non-dimensionalised stretch, bend, kinetic, and total energies: "
                << nonDimStretchEnergy << ", " << nonDimBendEnergy << ", " << nonDimKineticEnergy
                << ", " << nonDimKineticEnergy << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy;
    }

    /* If angle deficits were calculated, print the sum of the non-boundary
    and boundary node deficits separately, corresponding to the integrated
    Gauss curvature and integrated geodesic curvature terms in Gauss-Bonnet
    respectively.*/
    if (settings.getCore().isAngleDeficitPrinted()) {
        outFile << ", total interior angle deficit = " << kahanSum(interiorNodeAngleDeficits)
                << ", total boundary angle deficit = " << kahanSum(boundaryNodeAngleDeficits);
    }


    // Continue with rest of preamble and then output the relevant data to file.
    outFile << "\n" << "ASCII" << "\n"
            << "DATASET POLYDATA" << "\n"
            << "POINTS " << nodes.size() << " double" << "\n";

    for (int i = 0; i < nodes.size(); ++i) {
        outFile << nodes[i].pos(0) << " " << nodes[i].pos(1) << " " << nodes[i].pos(2) << "\n";
    }

    outFile << "POLYGONS " << triangles.size() << " " << 4 * triangles.size() << "\n";

    for (int i = 0; i < triangles.size(); ++i) {
        outFile << "3 " << triangles[i].vertexLabels(0) << " " << triangles[i].vertexLabels(1) << " "
                << triangles[i].vertexLabels(2) << "\n";
    }

    outFile << "CELL_DATA " << triangles.size() << "\n";

    outFile << "SCALARS gaussCurv double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < triangles.size(); ++i) {
        outFile << gaussCurvatures[i] << "\n";
    }

    outFile << "SCALARS meanCurv double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < triangles.size(); ++i) {
        outFile << meanCurvatures[i] << "\n";
    }

    /* Now print non-dimensionalised stretch and bend energies, if specified in
    settings file. Print also our strain measure, which is not completely
    dissimilar to a non-dimensional stretch energy. Also, Cauchy stress info.*/
    if (settings.getCore().isEnergyPrinted()) {

        outFile << "SCALARS nonDimStretchEnergyDensity double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < triangles.size(); ++i) {
            outFile << stretchEnergyDensities[i] / settings.getStretchEnergyDensityScale()  << "\n";
        }

        outFile << "SCALARS nonDimBendEnergyDensity double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < triangles.size(); ++i) {
            outFile << bendEnergyDensities[i] / settings.getStretchEnergyDensityScale() << "\n";
        }

        outFile << "SCALARS strainMeasure double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < triangles.size(); ++i) {
            outFile << strainMeasures[i] << "\n";
        }


        outFile << "FIELD cauchyStressInfo 4" << "\n";

        outFile << "dimlessCauchyStressEigenval1 1 " << triangles.size() << "double" << "\n";
        for (int i = 0; i < triangles.size(); ++i) {
            outFile << cauchyStressEigenvals[i](0) / (settings.getCore().getShearModulus() * settings.getCore().getThickness()) << "\n";
        }
        outFile << "cauchyStressEigenvec1 3 " << triangles.size() << "double" << "\n";
        for (int i = 0; i < triangles.size(); ++i) {
            outFile << cauchyStressEigenvecs[i](0, 0) << " " << cauchyStressEigenvecs[i](1, 0) << " "
                    << cauchyStressEigenvecs[i](2, 0) << "\n";
        }
        outFile << "dimlessCauchyStressEigenval2 1 " << triangles.size() << "double" << "\n";
        for (int i = 0; i < triangles.size(); ++i) {
            outFile << cauchyStressEigenvals[i](1) / (settings.getCore().getShearModulus() * settings.getCore().getThickness()) << "\n";
        }
        outFile << "cauchyStressEigenvec2 3 " << triangles.size() << "double" << "\n";
        for (int i = 0; i < triangles.size(); ++i) {
            outFile << cauchyStressEigenvecs[i](0, 1) << " " << cauchyStressEigenvecs[i](1, 1) << " "
                    << cauchyStressEigenvecs[i](2, 1) << "\n";
        }
    }

    // Now print current triangle areas if specified in settings file.
    if (settings.getCore().isTriangleAreaPrinted()) {

        outFile << "SCALARS triArea double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < triangles.size(); ++i) {
            outFile << 1.0 / triangles[i].currAreaInv << "\n";
        }
    }

    // Now print radii triangle centroids are at in reference state.
    outFile << "SCALARS triRefRadialCoord double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < triangles.size(); ++i) {
        outFile << triangles[i].refCentroid.norm() << "\n";
    }

    // Now print angle deficits at nodes, if specified in settings file.
    if (settings.getCore().isAngleDeficitPrinted()) {

        outFile << "POINT_DATA " << nodes.size() << "\n";

        outFile << "SCALARS angleDeficit double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < nodes.size(); ++i) {
            outFile << angleDeficits[i] << "\n";
        }
    }

    outFile.close();
}
