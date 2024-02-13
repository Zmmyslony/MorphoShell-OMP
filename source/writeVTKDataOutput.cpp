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

void writeVTKDataOutput(
        const std::vector<Node> &nodes,
        const std::vector<Triangle> &triangles,
        const int &stepcount,
        const double &time,
        const double &currDialInFactor,
        const size_t &progTensorSequenceCounter,
        const std::vector<double> &gaussCurvatures,
        const std::vector<double> &meanCurvatures,
        const std::vector<double> &angleDeficits,
        const std::vector<double> &interiorNodeAngleDeficits,
        const std::vector<double> &boundaryNodeAngleDeficits,
        const std::vector<double> &stretchEnergyDensities,
        const std::vector<double> &bendEnergyDensities,
        const std::vector<double> &stretchEnergies,
        const std::vector<double> &bendEnergies,
        const std::vector<double> &kineticEnergies,
        const std::vector<double> &strainMeasures,
        const std::vector<Eigen::Vector2d> &cauchyStressEigenvals,
        const std::vector<Eigen::Matrix<double, 3, 2> > &cauchyStressEigenvecs,
        const Settings &settings,
        const std::string &outputDirName) {

    std::ofstream outFile(outputDirName + "/stepcount_" + std::to_string(stepcount) + "_output.vtk");
//    std::stringstream outFile;
    if (!outFile) {
        throw std::runtime_error("Error: Problem creating or opening output file.");
    }


    outFile << std::scientific << std::setprecision(14)
            << "# vtk DataFile Version 4.2" << "\n"
            << "currDialInFactor = " << currDialInFactor << " dialling programmed tensors _"
            << progTensorSequenceCounter + 1 << ", time = " << time << ", stepcount = " << stepcount;

    /* If specified in settings file, also print total stretching and bending
    energies here.*/
    if (settings.is_energy_densities_printed) {
        double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.char_stretch_energy_scale;
        double nonDimBendEnergy = kahanSum(bendEnergies) / settings.char_stretch_energy_scale;
        double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.char_stretch_energy_scale;
        outFile << ", non-dimensionalised stretch, bend, kinetic, and total energies: "
                << nonDimStretchEnergy << ", " << nonDimBendEnergy << ", " << nonDimKineticEnergy
                << ", " << nonDimKineticEnergy << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy;
    }

    /* If angle deficits were calculated, print the sum of the non-boundary
    and boundary node deficits separately, corresponding to the integrated
    Gauss curvature and integrated geodesic curvature terms in Gauss-Bonnet
    respectively.*/
    if (settings.is_angle_deficits_printed) {
        outFile << ", total interior angle deficit = " << kahanSum(interiorNodeAngleDeficits)
                << ", total boundary angle deficit = " << kahanSum(boundaryNodeAngleDeficits);
    }


    // Continue with rest of preamble and then output the relevant data to file.
    outFile << "\n" << "ASCII" << "\n"
            << "DATASET POLYDATA" << "\n"
            << "POINTS " << settings.num_nodes << " double" << "\n";

    for (int i = 0; i < settings.num_nodes; ++i) {
        outFile << nodes[i].pos(0) << " " << nodes[i].pos(1) << " " << nodes[i].pos(2) << "\n";
    }

    outFile << "POLYGONS " << settings.num_triangles << " " << 4 * settings.num_triangles << "\n";

    for (int i = 0; i < settings.num_triangles; ++i) {
        outFile << "3 " << triangles[i].vertexLabels(0) << " " << triangles[i].vertexLabels(1) << " "
                << triangles[i].vertexLabels(2) << "\n";
    }

    outFile << "CELL_DATA " << settings.num_triangles << "\n";

    outFile << "SCALARS gaussCurv double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < settings.num_triangles; ++i) {
        outFile << gaussCurvatures[i] << "\n";
    }

    outFile << "SCALARS meanCurv double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < settings.num_triangles; ++i) {
        outFile << meanCurvatures[i] << "\n";
    }

    /* Now print non-dimensionalised stretch and bend energies, if specified in
    settings file. Print also our strain measure, which is not completely
    dissimilar to a non-dimensional stretch energy. Also, Cauchy stress info.*/
    if (settings.is_energy_densities_printed) {

        outFile << "SCALARS nonDimStretchEnergyDensity double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << stretchEnergyDensities[i] / settings.char_stretch_energy_density_scale << "\n";
        }

        outFile << "SCALARS nonDimBendEnergyDensity double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << bendEnergyDensities[i] / settings.char_stretch_energy_density_scale << "\n";
        }

        outFile << "SCALARS strainMeasure double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << strainMeasures[i] << "\n";
        }


        outFile << "FIELD cauchyStressInfo 4" << "\n";

        outFile << "dimlessCauchyStressEigenval1 1 " << settings.num_triangles << "double" << "\n";
        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << cauchyStressEigenvals[i](0) / (settings.shear_modulus * settings.sheet_thickness) << "\n";
        }
        outFile << "cauchyStressEigenvec1 3 " << settings.num_triangles << "double" << "\n";
        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << cauchyStressEigenvecs[i](0, 0) << " " << cauchyStressEigenvecs[i](1, 0) << " "
                    << cauchyStressEigenvecs[i](2, 0) << "\n";
        }
        outFile << "dimlessCauchyStressEigenval2 1 " << settings.num_triangles << "double" << "\n";
        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << cauchyStressEigenvals[i](1) / (settings.shear_modulus * settings.sheet_thickness) << "\n";
        }
        outFile << "cauchyStressEigenvec2 3 " << settings.num_triangles << "double" << "\n";
        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << cauchyStressEigenvecs[i](0, 1) << " " << cauchyStressEigenvecs[i](1, 1) << " "
                    << cauchyStressEigenvecs[i](2, 1) << "\n";
        }
    }

    // Now print current triangle areas if specified in settings file.
    if (settings.is_triangle_areas_printed) {

        outFile << "SCALARS triArea double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < settings.num_triangles; ++i) {
            outFile << 1.0 / triangles[i].currAreaInv << "\n";
        }
    }

    // Now print radii triangle centroids are at in reference state.
    outFile << "SCALARS triRefRadialCoord double 1" << "\n"
            << "LOOKUP_TABLE default" << "\n";

    for (int i = 0; i < settings.num_triangles; ++i) {
        outFile << triangles[i].refCentroid.norm() << "\n";
    }

    // Now print angle deficits at nodes, if specified in settings file.
    if (settings.is_angle_deficits_printed) {

        outFile << "POINT_DATA " << settings.num_nodes << "\n";

        outFile << "SCALARS angleDeficit double 1" << "\n"
                << "LOOKUP_TABLE default" << "\n";

        for (int i = 0; i < settings.num_nodes; ++i) {
            outFile << angleDeficits[i] << "\n";
        }
    }

    outFile.close();
}
