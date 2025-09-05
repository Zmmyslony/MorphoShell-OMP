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

Function that takes a data filename and reads node coordinates and the
node labels of each polygon into prepared std:vectors of the node and triangle
classes. The number of nodes and triangles is also put into the settings struct.
The programmed metric tensors and second fundamental forms are then also read
into the appropriate data structures.
The file is assumed to correspond to a text file containing VTK Legacy
PolyData (triangles specifically). See VTK documentation for explanation of
format. NODE LABELS START AT ZERO. In fact the assumed format is even more
stringent than just vtk (vtk reader would be required otherwise); it must be
exactly as in the example data files accompanying this code.

If a filename is supplied as a third command line argument, a dial-in factor,
counter into the programmed tensor sequence, and a set of node positions are read
in from that file. This is to provide the option of beginning evolution from an
ansatz, which could for example be where a previous unfinished simulation left
off. The format for such ansatz files must therefore be exactly as is used for
the OUTPUT vtk files, except that only the initial preamble (minus time and
stepcount), and the node position data are required.*/

//For testing (spaces can be fiddly with ignore()!) can use:
/*
std::string teststring;
init_DataFile >> teststring;
std::cout << teststring << std::endl;
*/

//Turn Eigen bounds checking off for speed (after running with checks naturally)
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream>
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <stdexcept>
#include <fstream>

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include "readVtk.hpp"
#include "Node.hpp"
#include "Triangle.hpp"


vtkSmartPointer<vtkPolyData> ReadPolyData(const char *fileName) {
    vtkNew<vtkPolyDataReader> reader;
    reader->SetFileName(fileName);
    reader->Update();
    return reader->GetOutput();
}


void readVTKData(std::vector<Node> &nodes, std::vector<Triangle> &triangles,
                 std::vector<std::vector<Eigen::Vector3d> > &sequenceOf_ProgMetricInfo,
                 std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_InvProgMetrics,
                 std::vector<std::vector<double> > &sequenceOf_ProgTaus,
                 std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_ProgSecFFs,
                 bool is_lce_mode_enabled, const std::string &init_data_file_name_str,
                 std::size_t &progTensorSequenceCounterToStartFrom, double &dialInFactorToStartFrom,
                 std::vector<Eigen::Vector3d> &nodeAnsatzPositions, const std::string &ansatz_data_file_name_str,
                 const CoreConfig &config) {
    auto main_vtk_data = ReadPolyData((const char *) init_data_file_name_str.c_str());
    vtkSmartPointer<vtkPolyData> ansatz_vtk_data;
    if (ansatz_data_file_name_str != "no_ansatz_file") {
        ansatz_vtk_data = ReadPolyData((const char *) ansatz_data_file_name_str.c_str());
    }

    unsigned int vertex_count = main_vtk_data->GetNumberOfPoints();
    unsigned int triangle_count = main_vtk_data->GetNumberOfCells();

    // Read in the reference node positions
    nodes.reserve(vertex_count);
    for (unsigned int i = 0; i < vertex_count; i++) {
        double coords[3];
        main_vtk_data->GetPoint(i, coords);
        nodes.emplace_back(i, coords);
    }

    // If ansatz is specified, read in the ansatz nodes too
    if (ansatz_data_file_name_str != "no_ansatz_file") {
        if (vertex_count != ansatz_vtk_data->GetNumberOfPoints() ||
            triangle_count != ansatz_vtk_data->GetNumberOfCells()) {
            throw std::runtime_error(
                "VTK files for programmed quantities and ansatz have different count of vertices or polygons.");
        }

        nodeAnsatzPositions.reserve(vertex_count);
        for (unsigned int i = 0; i < vertex_count; i++) {
            double coords[3];
            ansatz_vtk_data->GetPoint(i, coords);
            nodeAnsatzPositions.emplace_back(Eigen::Vector3d({coords[0], coords[1], coords[2]}));
        }
    }

    // Read in the triangle structure
    triangles.reserve(triangle_count);
    for (unsigned int i = 0; i < triangle_count; i++) {
        vtkNew<vtkIdList> ids;
        main_vtk_data->GetCellPoints(i, ids);
        if (ids->GetNumberOfIds() != 3) {
            throw std::runtime_error("Mesh error: non-triangular polygons are present in the mesh.");
        }
        triangles.emplace_back(i, ids->GetId(0), ids->GetId(1), ids->GetId(2), nodes);
    }


    // Read in the programmed quantities which are defined as cell data.
    auto cell_data = main_vtk_data->GetCellData();
    int metric_information_array_count = cell_data->GetNumberOfArrays();
    if (auto magnetisation_array = cell_data->GetArray("reference_magnetisation_density"); magnetisation_array != nullptr) {
        metric_information_array_count -= 1;
        for (int  i =0; i < triangle_count; i++ ) {
            triangles[i].setReferenceMagnetisationDensity({
                magnetisation_array->GetTuple(i)[0],
                magnetisation_array->GetTuple(i)[1],
                magnetisation_array->GetTuple(i)[2]
            });
        }
   }

    if (metric_information_array_count % 3 != 0) {
        throw std::runtime_error(
            "VTK: too many fields defining the deformation (" + std::to_string(metric_information_array_count) +
            ")- they should be provided in triplets"
            "(programmed_metric_info_i, programmed_bend_info_i, programmed_taus_i), with i starting at 0.");
    }
    if (metric_information_array_count == 0) {
        throw std::runtime_error("VTK: no fields defining the deformation are present in the input file.");
    }

    int progTensorCount = cell_data->GetNumberOfArrays() / 3;
    sequenceOf_ProgMetricInfo.resize(progTensorCount + 1);
    sequenceOf_InvProgMetrics.resize(progTensorCount + 1);
    sequenceOf_ProgTaus.resize(progTensorCount + 1);
    sequenceOf_ProgSecFFs.resize(progTensorCount + 1);


    for (unsigned int i = 1; i < progTensorCount + 1; i++) {
        auto metric_info = cell_data->GetArray(("programmed_metric_info_" + std::to_string(i - 1)).c_str());
        auto bend_info = cell_data->GetArray(("programmed_bend_info_" + std::to_string(i - 1)).c_str());
        auto tau_info = cell_data->GetArray(("programmed_taus_" + std::to_string(i - 1)).c_str());
        if (metric_info == nullptr) {
            throw std::runtime_error("VTK does not contain programmed_metric_info_" + std::to_string(i - 1));
        }
        if (bend_info == nullptr) {
            throw std::runtime_error("VTK does not contain programmed_bend_info_" + std::to_string(i - 1));
        }
        if (tau_info == nullptr) {
            throw std::runtime_error("VTK does not contain programmed_taus_" + std::to_string(i - 1));
        }

        sequenceOf_ProgMetricInfo[i].resize(triangle_count);
        sequenceOf_InvProgMetrics[i].resize(triangle_count);
        sequenceOf_ProgTaus[i].resize(triangle_count);
        sequenceOf_ProgSecFFs[i].resize(triangle_count);

        Eigen::Vector3d tempProgMetricInfo;
        for (unsigned int j = 0; j < triangle_count; j++) {
            tempProgMetricInfo = {
                metric_info->GetTuple(j)[0],
                metric_info->GetTuple(j)[1],
                metric_info->GetTuple(j)[2]
            };
            if (is_lce_mode_enabled) {
                sequenceOf_ProgMetricInfo[i][j] = tempProgMetricInfo;
                // TODO - CHECK IF WE NEED TO ALSO INCLUDE THE INVERSE OF THE METRIC.
                //                /* In case the metric components will be dialled, rather than lambda directly,
                //                set the corresponding programmed inverse metric components accordingly.*/
                //                double dirAng = tempProgMetricInfo(0);
                //                double cosDirAng = cos(dirAng);
                //                double sinDirAng = sin(dirAng);
                //                double lambda = tempProgMetricInfo(1);
                //                double nu = tempProgMetricInfo(2);
                //                double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
                //                double lambdaToThe2Nu = pow(lambda, 2.0 * nu);
                //
                //                sequenceOf_InvProgMetrics[i][j](0, 0) =
                //                        lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
                //                sequenceOf_InvProgMetrics[i][j](0, 1) =
                //                        (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
                //                sequenceOf_InvProgMetrics[i][j](1, 0) = sequenceOf_InvProgMetrics.at([i][j](0, 1);
                //                sequenceOf_InvProgMetrics[i][j](1, 1) =
                //                        lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
            } else {
                Eigen::Matrix<double, 2, 2> tempProgMetric;
                Eigen::FullPivLU<Eigen::Matrix<double, 2, 2> > lu;
                tempProgMetric(0, 0) = tempProgMetricInfo(0);
                tempProgMetric(0, 1) = tempProgMetricInfo(1);
                tempProgMetric(1, 0) = tempProgMetric(0, 1);
                tempProgMetric(1, 1) = tempProgMetricInfo(2);

                /* We actually store the inverse programmed metric for efficiency
                reasons, so invert the metric we just read in (after checking
                that it is actually invertible): */
                lu.compute(tempProgMetric);
                if (!lu.isInvertible()) {
                    throw std::runtime_error(
                        "Error: One of the programmed metrics that was read in was not invertible. Aborting.");
                }
                sequenceOf_InvProgMetrics[i][j] = tempProgMetric.inverse();
            }

            sequenceOf_ProgSecFFs[i][j](0, 0) = bend_info->GetTuple(j)[0];
            sequenceOf_ProgSecFFs[i][j](0, 1) = bend_info->GetTuple(j)[1];
            sequenceOf_ProgSecFFs[i][j](1, 0) = bend_info->GetTuple(j)[1];
            sequenceOf_ProgSecFFs[i][j](1, 1) = bend_info->GetTuple(j)[2];

            sequenceOf_ProgTaus[i][j] = tau_info->GetTuple(j)[0];
        }
    }

    // Read in clamping information
    auto point_data = main_vtk_data->GetPointData();

    auto clamp_and_load_indicators = point_data->GetArray("clamp_And_Load_Indicators");
    for (unsigned int i = 0; i < vertex_count; i++) {
        if (clamp_and_load_indicators->GetTuple(i)[0]) {
            nodes[i].clamp(config);
        }
        nodes[i].isLoadForceEnabled = (bool) clamp_and_load_indicators->GetTuple(i)[1];
    }


    // Read in the dial-in factor and programmed quantity counter.
    if (ansatz_data_file_name_str != "no_ansatz_file") {
        /* Temp variable used to check a value is non-negative
        while assigning it to a std::size_t which is unsigned.*/
        int tempProgTensorSequenceCounter;

        std::ifstream ansatz_DataFile(ansatz_data_file_name_str);
        if (!ansatz_DataFile) {
            throw std::runtime_error("Error: Problem opening ansatz data file given as third command line argument.");
        }
        // Ignore first line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        /* Get the dial-in factor and corresponding point in the programmed
        tensors sequence that evolution of the ansatz should start from.*/
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        ansatz_DataFile >> dialInFactorToStartFrom;
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '_');
        ansatz_DataFile >> tempProgTensorSequenceCounter;
        // Adjust offset so a single set of programmed tensors means dialling starts from the trivial '_0' tensors etc.
        progTensorSequenceCounterToStartFrom = tempProgTensorSequenceCounter - 1;
    }
}
