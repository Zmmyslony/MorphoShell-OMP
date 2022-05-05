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
#include <cstddef>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <Eigen/Dense>
#include <stdexcept>

#include "readVTKData.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "SettingsStruct.hpp"
#include "CustomOutStreamClass.hpp"

void readVTKData(
        std::vector<Node> &nodes,
        std::vector<Triangle> &triangles,
        std::vector<std::vector<Eigen::Vector3d> > &sequenceOf_ProgMetricInfo,
        std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_InvProgMetrics,
        std::vector<std::vector<double> > &sequenceOf_ProgTaus,
        std::vector<std::vector<Eigen::Matrix<double, 2, 2> > > &sequenceOf_ProgSecFFs,
        SettingsStruct &settings,
        const std::string &init_data_file_name_str,
        std::size_t &progTensorSequenceCounterToStartFrom,
        double &dialInFactorToStartFrom,
        std::vector<Eigen::Vector3d> &nodeAnsatzPositions,
        const std::string &ansatz_data_file_name_str,
        CustomOutStreamClass &logStream) {

    /* Variable to hold number of different sets of programmed tensors are
    present in the input file. These are all stored and then activated in
    sequence in the dynamics.*/
    int numProgTensorsInSequence;

    /* Temporary string variable used to help check the data files have the
    correct format. */
    std::string tempString;

    logStream.open();
    logStream << "Beginning reading of input non-ansatz data file." << std::endl;
    logStream.close();


    std::ifstream init_DataFile(init_data_file_name_str);
    if (!init_DataFile) {
        throw std::runtime_error("Error: Problem opening data file given as second command line argument.");
    }

    // Ignore first 4 lines.
    for (int i = 0; i < 4; ++i) {
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Ignore "POINTS".
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');

    // Get number of nodes and resize nodes container accordingly.
    init_DataFile >> settings.NumNodes;
    if (settings.NumNodes <= 0) {
        throw std::runtime_error(
                "Error: Problem with (or before) line giving number of nodes in non-ansatz data file.");
    } else {
        nodes.resize(settings.NumNodes);
    }

    // Ignore rest of line.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Put node label and coordinates in nodes container.
    for (int nodeLabel = 0; nodeLabel < settings.NumNodes; ++nodeLabel) {
        nodes.at(nodeLabel).label = nodeLabel;
        init_DataFile >> nodes.at(nodeLabel).pos(0);
        init_DataFile >> nodes.at(nodeLabel).pos(1);
        init_DataFile >> nodes.at(nodeLabel).pos(2);
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n'); // Ignore delimeter at end of line.
    }


    // Check we've arrived at the expected line of the file: 'POLYGONS'
    init_DataFile >> tempString;
    if (tempString.compare("POLYGONS")) {
        throw std::runtime_error(
                "Error: Problem with non-ansatz data file, at or before 'POLYGONS' line. E.g. may have provided or stated wrong number of nodes.");
    }
    tempString.clear();


    /* Get number of triangles and resize triangles, director angles and
    director twists containers accordingly. */
    init_DataFile >> settings.NumTriangles;
    if (settings.NumTriangles <= 0) {
        throw std::runtime_error("Error: Problem with line giving number of triangles in non-ansatz data file.");
    } else {
        triangles.resize(settings.NumTriangles);
    }

    // Ignore rest of line.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Put triangles' self and vertex labels in triangles container.
    for (int triangleLabel = 0; triangleLabel < settings.NumTriangles; ++triangleLabel) {
        triangles.at(triangleLabel).label = triangleLabel;
        /* Ignore first number of each row, which just says this polygon has
        3 vertices.*/
        init_DataFile.ignore(1);
        init_DataFile >> triangles.at(triangleLabel).vertexLabels(0);
        init_DataFile >> triangles.at(triangleLabel).vertexLabels(1);
        init_DataFile >> triangles.at(triangleLabel).vertexLabels(2);
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    // Check we've arrived at the expected line of the file: 'CELL_DATA'...
    init_DataFile >> tempString;
    if (tempString.compare("CELL_DATA")) {
        throw std::runtime_error(
                "Error: Problem with non-ansatz data file, at or before 'CELL_DATA' line. E.g. may have provided or stated wrong number of triangles.");
    }
    tempString.clear();
    // Ignore rest of line.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Ignore 'FIELD programmed_Quantities'.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    // Record number of separate sets of programmed tensors have been given.
    init_DataFile >> numProgTensorsInSequence;
    if (numProgTensorsInSequence % 3 == 0) {
        numProgTensorsInSequence = numProgTensorsInSequence /
                                   3; //From this line onwards a programmed metric/secFF/tau triple is considered a single entry of the sequence
    } else {
        throw std::runtime_error(
                "Error: Stated no. of programmed metrics + no. of programmed secFFs in sequence (stated in input file after 'FIELD programmed_Tensors') is not even, suggesting mismatch.");
    }
    // Ignore rest of line.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    /* Resize data structures which will hold the sequence of programmed
    tensors. The +1 leaves the first in the sequence to be the trivial tensors
    for the initial flat plane. */
    sequenceOf_ProgMetricInfo.resize(numProgTensorsInSequence + 1);
    sequenceOf_InvProgMetrics.resize(numProgTensorsInSequence + 1);
    sequenceOf_ProgTaus.resize(numProgTensorsInSequence + 1);
    sequenceOf_ProgSecFFs.resize(numProgTensorsInSequence + 1);


    // Now loop to read in the full sequence of programmed tensors.
    for (int sequenceIdx = 1; sequenceIdx < numProgTensorsInSequence + 1; ++sequenceIdx) {

        sequenceOf_ProgMetricInfo.at(sequenceIdx).resize(settings.NumTriangles);
        sequenceOf_InvProgMetrics.at(sequenceIdx).resize(settings.NumTriangles);
        sequenceOf_ProgTaus.at(sequenceIdx).resize(settings.NumTriangles);
        sequenceOf_ProgSecFFs.at(sequenceIdx).resize(settings.NumTriangles);

        // Check we've arrived at the expected line of the file: 'progMetricInfo'...
        init_DataFile >> tempString;
        if (tempString.find("progMetricInfo") == std::string::npos) {
            throw std::runtime_error("Error: Problem with non-ansatz data file, at or before a 'progMetricInfo' line.");
        }
        tempString.clear();
        // Ignore rest of line.
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Read triangles' programmed (energetically favoured) metric info.
        Eigen::Vector3d tempProgMetricInfo;
        Eigen::Matrix<double, 2, 2> tempProgMetric;
        Eigen::FullPivLU<Eigen::Matrix<double, 2, 2> > lu;
        for (int triangleLabel = 0; triangleLabel < settings.NumTriangles; ++triangleLabel) {

            init_DataFile >> tempProgMetricInfo(0);
            init_DataFile >> tempProgMetricInfo(1);
            init_DataFile >> tempProgMetricInfo(2);

            if (settings.isLCEModeEnabled) {

                sequenceOf_ProgMetricInfo.at(sequenceIdx).at(triangleLabel) = tempProgMetricInfo;

                /* In case the metric components will be dialled, rather than lambda directly,
                set the corresponding programmed inverse metric components accordingly.*/
                double dirAng = tempProgMetricInfo(0);
                double cosDirAng = cos(dirAng);
                double sinDirAng = sin(dirAng);
                double lambda = tempProgMetricInfo(1);
                double nu = tempProgMetricInfo(2);
                double lambdaToTheMinus2 = 1.0 / (lambda * lambda);
                double lambdaToThe2Nu = pow(lambda, 2.0 * nu);

                sequenceOf_InvProgMetrics.at(sequenceIdx).at(triangleLabel)(0, 0) =
                        lambdaToTheMinus2 * cosDirAng * cosDirAng + lambdaToThe2Nu * sinDirAng * sinDirAng;
                sequenceOf_InvProgMetrics.at(sequenceIdx).at(triangleLabel)(0, 1) =
                        (lambdaToTheMinus2 - lambdaToThe2Nu) * sinDirAng * cosDirAng;
                sequenceOf_InvProgMetrics.at(sequenceIdx).at(triangleLabel)(1, 0) = sequenceOf_InvProgMetrics.at(
                        sequenceIdx).at(triangleLabel)(0, 1);
                sequenceOf_InvProgMetrics.at(sequenceIdx).at(triangleLabel)(1, 1) =
                        lambdaToThe2Nu * cosDirAng * cosDirAng + lambdaToTheMinus2 * sinDirAng * sinDirAng;
            } else {
                // Note, before I didn't use tempProgMetric, I just used
                // triangles.at(triangleLabel).invProgMetric = triangles.at(triangleLabel).invProgMetric.inverse()
                // but this is WRONG due to the way Eigen works (aliasing)!
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
                sequenceOf_InvProgMetrics.at(sequenceIdx).at(triangleLabel) = tempProgMetric.inverse();
            }

            init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }


        // Check we've arrived at the expected line of the file: 'progSecFFComps'...
        init_DataFile >> tempString;
        if (tempString.find("progSecFFComps") == std::string::npos) {
            throw std::runtime_error(
                    "Error: Problem with non-ansatz data file, at or before a 'progSecFFComps' line. E.g. may have provided or stated wrong number of programmed metrics.");
        }
        tempString.clear();
        // Ignore rest of line.
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        /* Put triangles' programmed (energetically favoured) second fundamental
        forms in their container.*/
        for (int triangleLabel = 0; triangleLabel < settings.NumTriangles; ++triangleLabel) {
            init_DataFile >> sequenceOf_ProgSecFFs.at(sequenceIdx).at(triangleLabel)(0, 0);
            init_DataFile >> sequenceOf_ProgSecFFs.at(sequenceIdx).at(triangleLabel)(0, 1);
            sequenceOf_ProgSecFFs.at(sequenceIdx).at(triangleLabel)(1, 0) = sequenceOf_ProgSecFFs.at(sequenceIdx).at(
                    triangleLabel)(0, 1);
            init_DataFile >> sequenceOf_ProgSecFFs.at(sequenceIdx).at(triangleLabel)(1, 1);
            init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }


        // Check we've arrived at the expected line of the file: 'progTaus'...
        init_DataFile >> tempString;
        if (tempString.find("progTaus") == std::string::npos) {
            throw std::runtime_error(
                    "Error: Problem with non-ansatz data file, at or before a 'progTaus' line. E.g. may have provided or stated wrong number of programmed secFFs.");
        }
        tempString.clear();
        // Ignore rest of line.
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Put triangles' programmed tau values in their container.
        for (int triangleLabel = 0; triangleLabel < settings.NumTriangles; ++triangleLabel) {
            init_DataFile >> sequenceOf_ProgTaus.at(sequenceIdx).at(triangleLabel);
            init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    }

    // Check we've arrived at the expected line of the file: 'POINT_DATA'...
    init_DataFile >> tempString;
    if (tempString.compare("POINT_DATA")) {
        throw std::runtime_error(
                "Error: Problem with non-ansatz data file, at or before 'POINT_DATA'. E.g. may have provided or stated wrong number of programmed tau values.");
    }
    tempString.clear();
    // Ignore rest of line.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    // Ignore next 2 lines preambling point data.
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // Put node clamp and load indicators in nodes container.
    // std::ios_base::boolalpha ensures that '0's and '1's  will be interpreted as bools.
    init_DataFile.unsetf(std::ios_base::boolalpha);
    for (int nodeLabel = 0; nodeLabel < settings.NumNodes; ++nodeLabel) {
        init_DataFile >> nodes.at(nodeLabel).isClamped;
        init_DataFile >> nodes.at(nodeLabel).isLoadForceEnabled;
        init_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    /* Close the data file, as we should be at the end now. We check we
    actually did reach the end first, and for other errors that can be caught
    automatically by a std::ifstream. */
    if (init_DataFile.bad() || init_DataFile.fail()) {
        throw std::runtime_error(
                "Error: Unkown problem with non-ansatz data file. One of .bad(), .fail() was true - investigate! "
                "E.g. you may not have supplied enough node clamp/load indicator data in the file, so the file may have reached its end prematurely.");
    }
    init_DataFile >> tempString; // Attempts to read whatever's left, which should be nothing!
    if (!init_DataFile.eof()) {
        throw std::runtime_error(
                "Error: Did not reach end of non-ansatz data file as expected. Check that it has exactly the correct format.");
    }
    init_DataFile.close();

    logStream.open();
    logStream
            << "Completed reading of non-ansatz data file. Note; it is still possible that something went wrong in the process.\n"
            << std::endl;
    logStream.close();

    //////////////////////////////////////////////////////////////////////
    /* Now read the file with anzatz node positions if given as third command
    line argument. Only the preamble and node positions data is read;
    whatever comes after that is ignored.*/

    if (ansatz_data_file_name_str != "no_ansatz_file") {

        /* Temp variable used to check a value is non-negative
        while assigning it to a std::size_t which is unsigned.*/
        int tempProgTensorSequenceCounter;

        logStream.open();
        logStream << "Beginning reading of ansatz data file." << std::endl;
        logStream.close();

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
        // Ignore rest of line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Check a few things that could have gone wrong.
        if (dialInFactorToStartFrom < 0 || dialInFactorToStartFrom > 1 ||
            tempProgTensorSequenceCounter <= 0 || tempProgTensorSequenceCounter > numProgTensorsInSequence) {
            throw std::runtime_error(
                    "Error: Problem reading ansatz data file. Remember it must have exactly the correct format.");
        }

        // Ignore two lines of preamble.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        // Check number of nodes stated in ansatz file matches settings.NumNodes.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
        int ansatzStatedNumNodes;
        ansatz_DataFile >> ansatzStatedNumNodes;
        if (!(ansatzStatedNumNodes == settings.NumNodes)) {
            throw std::runtime_error(
                    "Error: Problem with (or before) line giving number of nodes in ansatz data file. "
                    "E.g. you may have supplied ansatz and initial data files with inconsistent numbers of nodes.");
        }
        // Ignore rest of line.
        ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');


        // Put node ansatz coordinates into suitable container.
        nodeAnsatzPositions.resize(settings.NumNodes);
        for (int nodeLabel = 0; nodeLabel < settings.NumNodes; ++nodeLabel) {
            ansatz_DataFile >> nodeAnsatzPositions.at(nodeLabel)(0);
            ansatz_DataFile >> nodeAnsatzPositions.at(nodeLabel)(1);
            ansatz_DataFile >> nodeAnsatzPositions.at(nodeLabel)(2);
            ansatz_DataFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }


        /* Close the data file, ignoring any remaining data. We also check
        for other errors that can be caught automatically by a
        std::ifstream. */
        if (ansatz_DataFile.bad() || ansatz_DataFile.fail()) {
            throw std::runtime_error(
                    "Error: Unkown problem with ansatz data file. One of .bad(), .fail() was true - investigate! "
                    " E.g. you may not have supplied enough ansatz node positions, so the file reached its end prematurely.");
        }
        ansatz_DataFile.close();

        logStream.open();
        logStream
                << "Completed reading of ansatz data file. Note, it is still possible that something went wrong in the process.\n"
                << std::endl;
        logStream.close();
    }
}
