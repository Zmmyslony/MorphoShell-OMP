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

#define BINARY_EXPORT

#include <Eigen/Dense>

#include <cstddef>
#include <string>
#include <fstream>
#include <iomanip> //for setting output precision etc
#include <vector>
#include <stdexcept>

#include "exportVtk.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "functions/kahanSum.hpp"

// Thanks to https://stackoverflow.com/questions/105252
template<typename T>
void SwapEnd(T &var) {
    char *varArray = reinterpret_cast<char *>(&var);
    for (long i = 0; i < static_cast<long>(sizeof(var) / 2); i++)
        std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
}

using double_vector = std::pair<std::string, std::vector<double>>;

std::stringstream convert_into_stream(const double_vector &data_pair) {
    std::stringstream stream;
    stream << "SCALARS " << data_pair.first << " double 1" << "\n"
           << "LOOKUP_TABLE default" << "\n";

    for (double i: data_pair.second) {
#ifdef BINARY_EXPORT
        SwapEnd(i);
        stream.write((char const *) &i, sizeof(double));
#else
        stream << i << " ";
#endif
    }
    stream << std::endl;
    return stream;
}

long long int
writeVTKDataOutput(const std::vector<Node> &nodes, const std::vector<Triangle> &triangles, const int &stepcount,
                   const double &time, const double &currDialInFactor, const size_t &progTensorSequenceCounter,
                   const std::vector<double> &gaussCurvatures, const std::vector<double> &meanCurvatures,
                   const std::vector<double> &angleDeficits, const std::vector<double> &interiorNodeAngleDeficits,
                   const std::vector<double> &boundaryNodeAngleDeficits, const std::vector<double> &stretchEnergies,
                   const std::vector<double> &bendEnergies, const std::vector<double> &kineticEnergies,
                   const std::vector<double> &strainMeasures, const std::vector<Eigen::Vector2d> &cauchyStressEigenvals,
                   const std::vector<Eigen::Matrix<double, 3, 2> > &cauchyStressEigenvecs, const SettingsNew &settings,
                   const std::string &outputDirName, double energy_loss) {
    auto begin = std::chrono::high_resolution_clock::now();

    std::ofstream outFile(outputDirName + "/step_count.vtk." + std::to_string(stepcount),
                          std::ios::out | std::ios::binary);
    if (!outFile) {
        throw std::runtime_error("Error: Problem creating or opening output file.");
    }
    std::stringstream header_stream;

    header_stream << std::scientific << std::setprecision(5)
                  << "# vtk DataFile Version 4.0" << "\n"
                  << "dial_in_factor = " << currDialInFactor << " dialling programmed tensors _"
                  << progTensorSequenceCounter + 1 << ", time = " << time << ", stepcount = " << stepcount;

    /* If specified in settings file, also print total stretching and bending
    energies here.*/
    if (settings.getCore().isEnergyPrinted()) {
        double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.getStretchEnergyScale();
        double nonDimBendEnergy = kahanSum(bendEnergies) / settings.getStretchEnergyScale();
        double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.getStretchEnergyScale();
        header_stream << ", non-dimensionalised stretch, bend, kinetic, and total energies: "
                      << nonDimStretchEnergy << ", " << nonDimBendEnergy << ", " << nonDimKineticEnergy
                      << ", " << nonDimKineticEnergy << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy;
    }
    header_stream << ", energy_loss: " << energy_loss;
    if (settings.getCore().isAngleDeficitPrinted()) {
        header_stream << ", total interior angle deficit = " << kahanSum(interiorNodeAngleDeficits)
                      << ", total boundary angle deficit = " << kahanSum(boundaryNodeAngleDeficits);
    }


    // Continue with rest of preamble and then output the relevant data to file.
#ifdef BINARY_EXPORT
    header_stream << std::endl << "BINARY" << std::endl;
#else
    header_stream << std::endl << "ASCII" << std::endl;
#endif
    outFile << header_stream.rdbuf();

    std::stringstream mesh_stream;
    mesh_stream << std::scientific << std::setprecision(5);
    mesh_stream << "DATASET POLYDATA" << std::endl
                << "POINTS " << nodes.size() << " double" << std::endl;

    for (const auto &node: nodes) {
#ifdef BINARY_EXPORT
        double x = node.pos(0);
        double y = node.pos(1);
        double z = node.pos(2);
        SwapEnd(x);
        SwapEnd(y);
        SwapEnd(z);

        mesh_stream.write((const char *) &x, sizeof(double));
        mesh_stream.write((const char *) &y, sizeof(double));
        mesh_stream.write((const char *) &z, sizeof(double));
#else
        mesh_stream << node.pos(0) << " " << node.pos(1) << " " << node.pos(2) << "\n";
#endif
    }

    mesh_stream << std::endl;

    mesh_stream << "POLYGONS " << triangles.size() << " " << 4 * triangles.size() << "\n";
    uint32_t element_count = 3;
    SwapEnd(element_count);
    for (const auto &triangle: triangles) {
#ifdef BINARY_EXPORT
        mesh_stream.write((const char *) &element_count, sizeof(uint32_t));
        auto labels = triangle.vertexLabels;
        uint32_t i0 = labels(0);
        uint32_t i1 = labels(1);
        uint32_t i2 = labels(2);
        SwapEnd(i0);
        SwapEnd(i1);
        SwapEnd(i2);

        mesh_stream.write((const char *) &i0, sizeof(uint32_t));
        mesh_stream.write((const char *) &i1, sizeof(uint32_t));
        mesh_stream.write((const char *) &i2, sizeof(uint32_t));
#else
        mesh_stream << "3 " << triangle.vertexLabels(0)
                    << " " << triangle.vertexLabels(1)
                    << " " << triangle.vertexLabels(2) << "\n";
#endif
    }
    mesh_stream << std::endl;
    outFile << mesh_stream.rdbuf();

    outFile << "POINT_DATA " << nodes.size() << "\n";

    double_vector gauss_curvature("gaussCurvature", std::vector<double>(nodes.size(), 0));
    double_vector mean_curvature("meanCurvature", std::vector<double>(nodes.size(), 0));
    double_vector force("force", std::vector<double>(nodes.size(), 0));
    double_vector stretch_energy_density("stretchEnergyDensityNonDim", std::vector<double>(nodes.size(), 0));
    double_vector bend_energy_density("bendEnergyDensityNonDim", std::vector<double>(nodes.size(), 0));
    double_vector node_area("nodeArea", std::vector<double>(nodes.size(), 0));
    double_vector strain_measure("strainMeasure", std::vector<double>(nodes.size(), 0));
    double_vector angle_deficits("angleDeficits", angleDeficits);


#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) {
        unsigned int incident_triangle_count = nodes[i].incidentTriLabels.size();
        for (auto j: nodes[i].incidentTriLabels) {
            gauss_curvature.second[i] += gaussCurvatures[j] / incident_triangle_count;
            mean_curvature.second[i] += meanCurvatures[j] / incident_triangle_count;
            strain_measure.second[i] += strainMeasures[j] / incident_triangle_count;
            const Triangle *triangle = &triangles[j];

            stretch_energy_density.second[i] += triangle->stretchEnergyDensity /
                                                (settings.getStretchEnergyDensityScale() / incident_triangle_count);
            bend_energy_density.second[i] += triangle->bendEnergyDensity /
                                             (settings.getStretchEnergyDensityScale() / incident_triangle_count);
            node_area.second[i] += (1 / triangle->currAreaInv) / incident_triangle_count;
        }
        force.second[i] = nodes[i].force.norm();
    }

    std::vector<double_vector *> data_to_export;
    data_to_export.emplace_back(&gauss_curvature);
    data_to_export.emplace_back(&mean_curvature);
    if (settings.getCore().isEnergyPrinted()) {
        data_to_export.emplace_back(&bend_energy_density);
        data_to_export.emplace_back(&stretch_energy_density);
        data_to_export.emplace_back(&strain_measure);
    }

    if (settings.getCore().isForcePrinted()) { data_to_export.emplace_back(&force); }
    if (settings.getCore().isAngleDeficitPrinted()) { data_to_export.emplace_back(&angle_deficits); }
    if (settings.getCore().isTriangleAreaPrinted()) { data_to_export.emplace_back(&node_area); }

    std::vector<std::stringstream> streams_to_export(data_to_export.size());
#pragma omp parallel for
    for (int i = 0; i < data_to_export.size(); i++) {
        streams_to_export[i] = convert_into_stream(*data_to_export[i]);
    }

    for (auto &stream: streams_to_export) {
        outFile << stream.rdbuf();
    }

    outFile.close();

    auto duration = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
}
