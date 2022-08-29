/* Author: Daniel Duffy, University of Cambridge, dld34@cam.ac.uk

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
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.`
/////////////////////////////////////////////////////

Libconfig++ is distributed under the Lesser GPL license (2.1 or later), and copyright is held by Mark A Lindner (at least in large part).
Version 1.7.2 was used for this code, which can be accessed via: http://hyperrealm.github.io/libconfig/

Eigen is distributed under the Mozilla Public License v. 2.0 , and copyright is held by Gael Guennebaud and Benoit Jacob (at least in large part).
Version 3.3.7 was used for this code, which can be accessed via: http://eigen.tuxfamily.org/

//////////////////////////////////////////////////////////////////////////

Main file for Shellmorph, the main code developed for my PhD.

This code simulates a thin shape-morphing 2D sheet with a programmed metric and
second fundamental form.
The main object in the simulations is an unstructured triangulated mesh, with
data stored in three std:vectors: 'nodes', 'triangles' and 'edges', each element
of which holds a class, in turn holding the data for that node/triangle/edge.
So for example nodes[0] is a Node instance holding the data for the first
node, which is labelled zero, matching the indexing of 'nodes'. The initial node
and triangle data is read in from a VTK Legacy file, which is also the output
format. Paraview is recommended for visualising these files.

The code executable takes either two or three command line arguments. The first
is a settings file, and the second is an input data file containing the planar
2D initial (reference) state mesh and the programmed quantities. The optional
third argument is an `ansatz' file (also .vtk), which can be an output file from
a previous simulation. If this is provided, the nodes will be moved to the
positions specified in the ansatz file before dynamics begins from that state.
Thus you can terminate a simulation and restart it from where it left off with
a new settings file (though the velocities start from zero upon restarting).

An elastic energy that penalises stretching and bending deviations from the
programmed metric and secondfundamental form is minimised (using analytical
derivatives) via either linearly-damped Newtonian dynamics or gradient descent.
The former is on the whole recommended, while the latter can be useful for
smoothing out meshes that have significant angles between triangles i.e. small
wavelength oscillations. Equilibrium is deemed achieved when two dimensionless
numbers fall below thresholds specified in the settings file. These numbers are
maxima (over the mesh) of non-dimensionalised node speed and elastic force.
*/

// Turn Eigen bounds checking off for speed (after running with checks naturally).
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#include <iostream> // Used for outputting messages to the terminal
#include <fstream> // Used for outputting to files
#include <cassert> // Used for debugging tools
#include <stdexcept> // Used for exception (error) handling
#include <cmath> // Used for some simple maths functions
#include <vector>// Used for some vectors
#include <iomanip> // For setting time and date format, and std::out precision if needed
#include <chrono>

#include <libconfig.h++> // Used for settings file
#include <Eigen/Dense> // Used for matrices of numbers


#include "CustomOutStreamClass.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "Edge.hpp"
#include "SettingsStruct.hpp"

#include "functions/getRealTime.hpp"
#include "functions/extract_Just_Filename.hpp"
#include "initialDirAndFileHandling.hpp"
#include "readSettingsFile.hpp"
#include "readVTKData.hpp"
#include "calculations/calcTriangleGeometries_and_DialledProgTensors.hpp"
#include "calculations/calcTrianglesIncidentOnNodes.hpp"
#include "calculations/calcTriangleAdjacencies_And_Edges.hpp"
#include "calculations/calcNodeNeighbours.hpp"
#include "calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp"
#include "setRemainingInitCond_and_NodeMasses.hpp"
#include "SimulationStatus.hpp"
#include "functions/perturbInitialPositionsWithRandomNoise.hpp"
#include "functions/zeroForces.hpp"
#include "calculations/calcSecFFsAndRelatedQuantities.hpp"
#include "calculations/calcDeformationForces.hpp"
#include "calculations/calcNonDeformationForces_and_ImposeBCS.hpp"
#include "advanceDynamics.hpp"
#include "writeVTKDataOutput.hpp"
#include "EquilibriumCheck.hpp"
#include "calculations/calcCurvatures.hpp"
#include "calculations/calcEnergiesAndStresses.hpp"
#include "functions/kahanSum.hpp"

//Create useful debug tools.
/*
//#define QUICKCHECK

#ifdef QUICKCHECK
assert( x > 100 );
#endif

//A second useful debug tool.
#define DISPLAY(x) std::cout <<  #x << "=" << x << " at line " << __LINE__ << std::endl;
*/

///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
/* Create a string to contain some things that will later be written to the log
file, as long as all goes to plan and the code gets that far. The first thing to
add to this string is a record of what command line arguments where given, and
at what date and time the run began.*/

    std::string initWriteToLogStr;

    initWriteToLogStr += "Simulation run began at: " + getRealTime() + "\n";
    initWriteToLogStr += "The command that was run was:\n";
    for (int i = 0; i < argc; ++i) {
        initWriteToLogStr += argv[i];
        initWriteToLogStr += " ";
    }
    initWriteToLogStr += "\n\n";


/* Check settings and data files given, with ansatz file optional, and get names
from command line. */
    if ((argc != 3) && (argc != 4)) {
        std::cerr << "Error: Please provide at least one settings file path "
                     "followed by one data file path as command line arguments.\n"
                     "A path to a file with ansatz node coordinates can be given as a third "
                     "argument if desired. All paths should be relative to the current working"
                     "directory, i.e. they should start with ./" << std::endl;
        return -1;
    }
    const char *settings_file_name = argv[1];
    std::string settings_file_name_str(settings_file_name);
    const char *initial_data_file_name = argv[2];
    std::string initial_data_file_name_str(initial_data_file_name);
    std::string ansatz_data_file_name_str("no_ansatz_file"); // Default, since usually no ansatz will be given

/* Will need the settings file name without any preceding directory info,
so we extract that first. Ditto for the initial data file file, and ansatz file.*/
    std::string settings_file_name_str_final_piece = extract_Just_Filename(settings_file_name_str);
    std::string initial_data_file_name_str_final_piece = extract_Just_Filename(initial_data_file_name_str);
    std::string ansatz_data_file_name_str_final_piece("no_ansatz_file");

//Get ansatz file name if given
    if (argc == 4) {
        const char *ansatz_data_file_name = argv[3];
        ansatz_data_file_name_str.assign(ansatz_data_file_name);
        ansatz_data_file_name_str_final_piece = extract_Just_Filename(ansatz_data_file_name_str);
        initWriteToLogStr += "Using ansatz data file " + ansatz_data_file_name_str_final_piece + "\n";
    }


// Do initial handling of directories and input files.
    std::string outputDirName("this_string_will_hold_the_output_dir_name");
    try {
        outputDirName = initialDirAndFileHandling(settings_file_name_str,
                                                  settings_file_name_str_final_piece,
                                                  initial_data_file_name_str,
                                                  initial_data_file_name_str_final_piece,
                                                  ansatz_data_file_name_str,
                                                  ansatz_data_file_name_str_final_piece,
                                                  argc,
                                                  initWriteToLogStr);
    }
    catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
        std::cerr << "Initial directory and file handling failed." << std::endl;
        return -1;
    }

/* Set up log file that will store a copy of the std::cout and std::cerr
output (with the exception of any writes to std::cerr that have already happened,
and the one below if the log file creation failed). The logStream will be passed
by reference to any functions that want to print output to the log file. The
logStream must be opened before writing to it, and it is good practice to close
it immediately afterwards.*/
    CustomOutStreamClass logStream;
    logStream.set_outputFileName(outputDirName + "/log.txt");
// Check we can create a log file with the chosen name.
    try {
        logStream.open();
    }
    catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
        return -1;
    }

/* Set format and precision of logStream, and print some things related to what
has already occurred above.*/
    logStream << std::defaultfloat << std::setprecision(6);
    logStream << initWriteToLogStr << std::endl;
    logStream.close();

/*Read in settings file (using libconfig++ library) and put these settings in
a SettingsStruct that will be passed around between functions when needed.
See SettingsStruct.hpp for details of settings. */
    SettingsStruct settings;
    try { readSettingsFile(settings, settings_file_name); }
    catch (const libconfig::FileIOException &fioex) {
        logStream.open();
        logStream << "I/O error reading settings file. Check file exists in correct directory and has correct name."
                  << std::endl;
        return -1;
    }
    catch (const libconfig::ParseException &pex) {
        logStream << "Parse error reading settings file. Check file has correct format." << std::endl;
        return -1;
    }
    catch (const std::runtime_error &error) {
        logStream << error.what() << std::endl;
        logStream.close();
        return -1;
    }
// Calculate Young's Modulus
    settings.YoungsModulus = 2.0 * settings.ShearModulus * (1.0 + settings.PoissonRatio);


////////////////////////////////////////////////////////////////////////////////

/* Read in *initial* data for each node and triangle
Each element of the node container 'nodes' is a class holding data for that node
etc. seriesOf_InvProgMetrics is a vector of vectors, where each element of the
outer vector holds a vector of (director angle, lambda, nu) vectors, one vector
for each triangle. The different sets held in the outer vector are to allow a
sequence of programmed tensors to be activated in turn, one after the other. For
example a bend first in one direction, with a bend in a second direction only
added afterwards. Ditto for sequenceOf_ProgSecFFs. Between each set of programmed
tensors the full dialling in procedure occurs. Unless an ansatz file is used,
the first round of dialling in will always start from the Identity and Zero for
the programmed metric and secFF respectively. For most uses of this code, only
a single set of programmed tensors will be supplied - multiple sets is a fairly
advanced case.*/
    std::vector <Node> nodes;
    std::vector <Triangle> triangles;
    std::vector <std::vector<Eigen::Vector3d>> sequenceOf_progMetricInfo;
    std::vector <std::vector<Eigen::Matrix<double, 2, 2> >> sequenceOf_InvProgMetrics;
    std::vector <std::vector<double>> sequenceOf_ProgTaus;
    std::vector <std::vector<Eigen::Matrix<double, 2, 2> >> sequenceOf_ProgSecFFs;


/* This variable keeps track of where in the sequence of programmed tensors the
dynamics should start from. The default (zero) is used unless an ansatz file was
given, in which case it may be that evolution from the ansatz should start from
further on in the sequence of programmed tensors. */
    std::size_t progTensorSequenceCounterToStartFrom = 0;
// This is the equivalent variable for currDialInFactor.
    double dialInFactorToStartFrom = 0.0;
// Create container to store node ansatz positions, if given.
    std::vector <Eigen::Vector3d> nodeAnsatzPositions;

    logStream.open();
    logStream << "Now attempting to read data files. An error here likely \n"
                 "implies a problem with a data file, for example a mismatch between the \n"
                 "number of nodes or triangles stated and the number actually present; \n"
                 "or other similar mismatches in numbers of data values; or a format problem. \n"
                 "Remember the input files must have exactly the correct format." << std::endl;
    logStream.close();
    try {
        readVTKData(nodes, triangles, sequenceOf_progMetricInfo, sequenceOf_InvProgMetrics, sequenceOf_ProgTaus,
                    sequenceOf_ProgSecFFs, settings, initial_data_file_name_str, progTensorSequenceCounterToStartFrom,
                    dialInFactorToStartFrom, nodeAnsatzPositions, ansatz_data_file_name_str, logStream);
    }
    catch (const std::out_of_range &out_of_bounds_error) {
        logStream << out_of_bounds_error.what() << std::endl;
        return -1;
    }
    catch (const std::runtime_error &reading_error) {
        logStream << reading_error.what() << std::endl;
        return -1;
    }




// Print numbers of nodes and triangles.
    logStream.open();
    logStream << "Number of nodes = " << settings.NumNodes << std::endl;
    logStream << "Number of triangles = " << settings.NumTriangles << std::endl;
/* Print warning if number of triangles is low - in this case boundary effects
will dominate, and the code should not be trusted, both due to the reduced
accuracy in the treatment of the boundary in e.g. the 2nd F.F. approx, and due
to possible bugs in this rather special case.*/
    if (settings.NumTriangles < 50) {
        logStream << "Your mesh has a small number of triangles. \nBeware that the code "
                     "is likely to be less accurate in this case, \nand unforeseen bugs are more "
                     "likely in extreme cases." << std::endl;
    }
    logStream.close();

/* Set the nodes and triangles debugging display functions to print to the
same output log file as logStream (could send somewhere else in principle).*/
    for (int i = 0; i < settings.NumNodes; ++i) {
        nodes[i].nodeLogStream.set_outputFileName(logStream.get_outputFileName());
    }
    for (int i = 0; i < settings.NumTriangles; ++i) {
        triangles[i].triLogStream.set_outputFileName(logStream.get_outputFileName());
    }


///////////////////////////////////////////////////////////////////////////////
//Here we do some remaining geometry/topology/mesh-related things

// Determine and store the labels of the triangles incident on each node.
    calcTrianglesIncidentOnNodes(nodes, triangles, settings);

/*Now calculate triangle edge-sharing adjacencies (triangle member data), and
set up an 'edges' data structure to store further edge information.*/
    std::vector <Edge> edges;
    try { calcTriangleAdjacencies_And_Edges(nodes, triangles, edges, settings); }
    catch (const std::runtime_error &mesh_adjacencies_or_edges_error) {
        logStream << mesh_adjacencies_or_edges_error.what() << std::endl;
        return -1;
    }

/* Set the edges debugging display functions to print to the same output log
file as logStream (could send somewhere else in principle).*/
    for (int i = 0; i < settings.NumEdges; ++i) {
        edges[i].edgeLogStream.set_outputFileName(logStream.get_outputFileName());
    }

/* Calculate and print the number of boundary and non-boundary edges.
Also label the nodes connected by boundary edges as boundary nodes.
Also calculate the total initial perimeter of the sample, to be used as a
characteristic sample length in estimating characteristic times, time steps
etc.*/
    int numBoundaryEdges = 0;
    std::vector<double> initBoundaryEdgeLengths(settings.NumEdges);

    for (int i = 0; i < settings.NumEdges; ++i) {
        if (edges[i].isOnBoundary) {

            numBoundaryEdges += 1;

            nodes[edges[i].nodeLabels(0)].isOnBoundary = true;
            nodes[edges[i].nodeLabels(1)].isOnBoundary = true;

            //If chosen in settings, clamp whole boundary in addition to clamp
            //indicators from data file.
            if (settings.isBoundaryClamped) {
                nodes[edges[i].nodeLabels(0)].isClamped = true;
                nodes[edges[i].nodeLabels(1)].isClamped = true;
            }

            // Store initial length of this boundary edge
            initBoundaryEdgeLengths[i] = (nodes[edges[i].nodeLabels(0)].pos - nodes[edges[i].nodeLabels(1)].pos).norm();
        }
    }

// Do sum to calculate perimeter, and set the characteristic sample length to it.
    double initPerimeter = kahanSum(initBoundaryEdgeLengths);
    settings.SampleCharLength = initPerimeter;
    logStream.open();
    logStream << "Initial perimeter = " << initPerimeter << std::endl;
    logStream.close();


//A further check that things are ok:
    try {
        if (3 * settings.NumTriangles != 2 * settings.NumEdges - numBoundaryEdges) {
            throw std::runtime_error(
                    "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
        }
    }
    catch (const std::runtime_error &mesh_adjacencies_or_edges_error) {
        logStream << mesh_adjacencies_or_edges_error.what() << std::endl;
        return -1;
    }

    logStream.open();
    logStream << "Number of edges = " << settings.NumEdges << std::endl;
    logStream << "Number of boundary edges = " << numBoundaryEdges << std::endl;
    logStream << "Number of non-boundary edges = " << settings.NumEdges - numBoundaryEdges << std::endl;
    logStream.close();


/* Now set triangle boundary indicators, based on whether they have any *edges*
on the sample boundary.*/
    int numBoundaryTriangles = 0;
    for (int i = 0; i < settings.NumTriangles; ++i) {
        if (triangles[i].isOnBoundary) {
            numBoundaryTriangles += 1;
        }
    }

    logStream.open();
    logStream << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
    logStream << "Number of holes in mesh = " << 1 + settings.NumEdges - settings.NumNodes - settings.NumTriangles
              << std::endl; // From Euler's formula for a planar graph.
    logStream.close();
    try {
        if (1 + settings.NumEdges - settings.NumNodes - settings.NumTriangles < 0) {
            throw std::runtime_error(
                    "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
        }
    }
    catch (const std::runtime_error &mesh_catastrophe_error) {
        logStream << mesh_catastrophe_error.what() << std::endl;
        return -1;
    }


/* Determine and store labels of the node neighbours to each node (those that
share an edge).*/
    calcNodeNeighbours(nodes, edges, settings);

// Calculate and store number of boundary nodes.
    int numBoundaryNodes = 0;
    for (int n = 0; n < settings.NumNodes; ++n) {
        if (nodes[n].isOnBoundary) {
            numBoundaryNodes += 1;
        }
    }
    settings.numBoundaryNodes = numBoundaryNodes;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/* Now calculate a 3x6 matrix for each triangle that can be stored and used
repeatedly to obtain the elements of the Second Fundamental Form for each
triangle during the simulation. The 2nd F.F. is an important extrinsic geometric
object, which completely specifies a surface (up to rigid motions) when combined
with a metric tensor (aka 1st F.F.). This is called the Bonnet Theorem. An
'intended' 2nd F.F. can be imposed with a nematic director that 'twists' through
the thickness of an LCE sheet, and thus the intended and actual 2nd F.F.s enter
into the energy that is being minimised in this simulation. Intuitively, the 2nd
F.F. describes the way the deformed surface is curved in 3D space in terms of
the directional derivative of the surface normal with respect to the 2D
parametrisation coordinates on the surface. Its matrix representation is
symmetric.
A related object is found by 'raising an index' by left-multiplying by the
inverse metric, which gives the 'Shape Operator'. This has determinant and trace
equal to the Gaussian curvature and twice the mean curvature respectively.
The approach taken to approximate the 2nd F.F. for a triangulated mesh in this
code is based on 'patch fitting', which is simple, and has been shown to
converge. Consider a triangle t. It has three vertices, and three other nearby
vertices are selected, subject to the following quadratic fitting being
well-conditioned.
The unique quadratic surface (a function of the 2D parametrisation coords) that
goes through all 6 of these points in the current state is then found and
assigned to triangle t. The 2nd F.F. for this 'patch' could be calculated in
full and averaged over the patch, although this complicates matters. Instead, we
just take the normal to the triangle face and combine that with the 2nd
derivatives of the patch surface following a particular expression for the secFF.
The matrix pre-calculated here is later multiplied onto the current node
positions of the 6 patch nodes to give patch surface derivatives.
*/
    logStream.open();
    logStream << "\n" << "Beginning patch selection and related pre-calculations." << std::endl;
    logStream.close();

    try { calc_nonVertexPatchNodes_and_MatForPatchDerivs(nodes, triangles, settings, logStream); }
    catch (const std::runtime_error &patch_error) {
        logStream << patch_error.what() << std::endl;
        return -1;
    }
    logStream.open();
    logStream << "Successfully completed patch setup." << "\n" << std::endl;
    logStream.close();


/*For the current and initial in-plane coordinate bases to match up, the node
labels for each triangle must be ordered such that the resulting calculated
normals for the initial flat state point along the +z axis, not -z. To do this
I calculate the normals for an arbitrary label order, then loop over triangles,
shuffling the labels for those with negative-z-component normals, before
recalculating the initial normals, areas etc. Note, a more complex procedure
will be required for non-planar initial geometries. The magic numbers in the
argument list are just debugging values that should not matter here.*/
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, waitingForEquilibrium, -12345, 98765,
                                                  sequenceOf_progMetricInfo, sequenceOf_InvProgMetrics,
                                                  sequenceOf_ProgTaus, sequenceOf_ProgSecFFs, settings);
    Eigen::Vector3d tempZAxisVec;
    tempZAxisVec << 0.0, 0.0, 1.0;
    for (int i = 0; i < settings.NumTriangles; ++i) {
        if (tempZAxisVec.dot(triangles[i].currSides.col(0).cross(triangles[i].currSides.col(1))) < 0) {
            int tempLabel = triangles[i].vertexLabels(2);
            triangles[i].vertexLabels(2) = triangles[i].vertexLabels(1);
            triangles[i].vertexLabels(1) = tempLabel;
        }
    }
    std::cout << "CHECK THIS - should shuffle normals before patch matrix calc I think?" << std::endl;
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, waitingForEquilibrium, -12345, 98765,
                                                  sequenceOf_progMetricInfo, sequenceOf_InvProgMetrics,
                                                  sequenceOf_ProgTaus, sequenceOf_ProgSecFFs, settings);

/* Set remaining initial conditions e.g. velocity, store initial triangle
side components and fractional edge lengths, and calculate node masses.
Also set the first set of programmed tensors to be the trivial ones for the
plane, and set up the triangles to have their tensors dialled in from these to
begin with. */
    try {
        setRemainingInitCond_and_NodeMasses(nodes, triangles, edges, sequenceOf_progMetricInfo,
                                            sequenceOf_InvProgMetrics, sequenceOf_ProgTaus, sequenceOf_ProgSecFFs,
                                            settings);
    }
    catch (const std::runtime_error &setInitCondAndInitSidesEtc_error) {
        logStream << setInitCondAndInitSidesEtc_error.what() << std::endl;
        return -1;
    }

/* First find the approximate smallest linear size of mesh element, based on
smallest altitude of each triangle. Find also the smallest value of
this linear size / sqrt(progTau), which is what really matters for the stretching
time step. If progTau varied a a lot over the sequenceOf_ProgTaus, you might
want to be less wasteful and recalculate the time step each time you move on to
the next set of progTaus in the sequence, because the time step then might be
much bigger for other sets in the sequence. If you wanted to get even fancier,
you could change the time step as prog_Tau is dialled in between two sets in the
sequence.*/
    settings.ApproxMinInitElemSize = 2.0 * triangles[0].initArea / triangles[0].currSides.col(0).norm();
    double smallestSizeOverRootTau = settings.ApproxMinInitElemSize / sequenceOf_ProgTaus[0][0];
    int smallestTri = 0;
    for (int i = 0; i < settings.NumTriangles; ++i) {

        // This variable will end up holding the smallest altitude for *this* tri.
        double smallestAltitude =
                2.0 * triangles[i].initArea / (triangles[i].currSides.col(0) - triangles[i].currSides.col(1)).norm();

        for (int s = 0; s < 2; ++s) {
            if (smallestAltitude > 2.0 * triangles[i].initArea / triangles[i].currSides.col(s).norm()) {
                smallestAltitude = 2.0 * triangles[i].initArea / triangles[i].currSides.col(s).norm();
            }
        }

        if (settings.ApproxMinInitElemSize > smallestAltitude) {
            settings.ApproxMinInitElemSize = smallestAltitude;
            smallestTri = i;
        }

        for (auto &sequenceOf_ProgTau: sequenceOf_ProgTaus) {
            if (smallestSizeOverRootTau > smallestAltitude / sqrt(sequenceOf_ProgTau[i])) {
                smallestSizeOverRootTau = smallestAltitude / sqrt(sequenceOf_ProgTau[i]);
            }
        }
    }
    logStream.open();
    logStream << "Sheet thickness = " << settings.SheetThickness << std::endl;
    logStream << "Approx smallest element linear size = " << settings.ApproxMinInitElemSize << " (triangle "
              << smallestTri << ")" << std::endl;
    logStream.close();


/* Calculate 'dialling in' time and damping coefficient based on toy model
stretching and bending analyses. The 'Long times' are approximate characteristic
times for the longest-wavelength modes in the system for stretching and bending.
The damping coefficient is chosen to approximately critically damp the longest-
wavelength bending mode.
*/
    const double PI = 3.14159265358979323846;

    logStream.open();
    double minWavevector = 2.0 * PI / settings.SampleCharLength;
    const double dampingScale =
            2.0 * sqrt(settings.InitDensity * settings.ShearModulus / (6.0 * (1.0 - settings.PoissonRatio))) *
            settings.SheetThickness * minWavevector * minWavevector;
    settings.NumDampFactor = settings.DampingPrefactor1 * dampingScale;
    logStream << "Numerical Damping Coefficient = " << settings.NumDampFactor << std::endl;

    double stretchingLongTime = sqrt(settings.InitDensity / settings.ShearModulus) / minWavevector;
    double bendingLongTime = sqrt(settings.InitDensity * 6.0 * (1.0 - settings.PoissonRatio) / settings.ShearModulus) /
                             (settings.SheetThickness * minWavevector * minWavevector);
    settings.BendingLongTime = bendingLongTime;

    if (stretchingLongTime > bendingLongTime) {
        settings.DialInStepTime = settings.DialInStepTimePrefactor * stretchingLongTime;
        logStream << "Dial-in step time = " << settings.DialInStepTime
                  << ", set based on stretching rather than bending." << std::endl;
    } else {
        settings.DialInStepTime = settings.DialInStepTimePrefactor * bendingLongTime;
        logStream << "Dial-in step time = " << settings.DialInStepTime
                  << ", set based on bending rather than stretching." << std::endl;
    }
    logStream.close();

/* Set settings.TimeBetweenEquilChecks based on the tunable dimensionless value
in the settings file, which relates this time to settings.DialInStepTime.*/
    settings.TimeBetweenEquilChecks = settings.TimeBetweenEquilChecksPrefactor * settings.DialInStepTime;


/* Calculate time step based on toy model stretching and bending analyses (take
 whichever gives shortest characteristic time), and print. */
    double stretchingTimeStep;
    double bendingTimeStep;
    if (!settings.isGradientDescentDynamicsEnabled) {

        stretchingTimeStep = settings.TimeStepPrefactor * sqrt(settings.InitDensity / settings.ShearModulus) *
                             smallestSizeOverRootTau;
        bendingTimeStep = settings.TimeStepPrefactor *
                          sqrt(settings.InitDensity * 6.0 * (1.0 - settings.PoissonRatio) / settings.ShearModulus) *
                          (settings.ApproxMinInitElemSize / settings.SheetThickness) * settings.ApproxMinInitElemSize /
                          (2.0 * PI);
    } else {
        /* Use second damping factor, since we should only be using gradient descent
        for already dialled-in states. Note the damping factor cancels out in the
        dynamics and has no effect - but only if you are consistent with which one
        you use!*/
        double numDampFac2 = settings.DampingPrefactor2 * dampingScale;
        stretchingTimeStep =
                settings.TimeStepPrefactor * (numDampFac2 / settings.ShearModulus) * smallestSizeOverRootTau *
                smallestSizeOverRootTau / (2.0 * PI * 2.0 * PI);
        bendingTimeStep = settings.TimeStepPrefactor * (6.0 * (1.0 - settings.PoissonRatio) * numDampFac2 /
                                                        (settings.SheetThickness * settings.SheetThickness *
                                                         settings.ShearModulus))
                          * smallestSizeOverRootTau * smallestSizeOverRootTau * smallestSizeOverRootTau *
                          smallestSizeOverRootTau / (2.0 * PI * 2.0 * PI * 2.0 * PI * 2.0 * PI);

        logStream.open();
        logStream
                << "\nUSING GRADIENT DESCENT DYNAMICS.\nBE WARNED - this feature is only intended for use on nearly converged states supplied as fully dialled-in ansatzes."
                << "\nThe number of timesteps required for a full simulation with gradient descent is expected to be huge (~10^8).\nAs such, the gradient descent long bending timescale "
                << "is not even used for DialInStepTime,\nas you really shouldn't be doing any dialling in!\nDialInStepTime is instead set such that usual non-gradient-descent settings"
                << " should still give reasonable printout and and equilibrium check frequencies." << std::endl;
        logStream.close();
    }

    logStream.open();
    logStream << "Short stretching and bending timescales: " << stretchingTimeStep << ", " << bendingTimeStep
              << std::endl;

    if (stretchingTimeStep < bendingTimeStep) {
        settings.TimeStep = stretchingTimeStep;
        logStream << "Time step = " << settings.TimeStep << ", set based on stretching rather than bending."
                  << std::endl;
    } else {
        settings.TimeStep = bendingTimeStep;
        logStream << "Time step = " << settings.TimeStep << ", set based on bending rather than stretching."
                  << std::endl;
    }
    logStream.close();

/* Fix up DialInStepTime if gradient descent is being used, just to give
reasonable printout and equilibrium check frequencies without extreme settings.*/
    if (settings.isGradientDescentDynamicsEnabled) {
        settings.DialInStepTime = 1000.0 * settings.TimeStep;
        settings.TimeBetweenEquilChecks = settings.TimeBetweenEquilChecksPrefactor * settings.DialInStepTime;
    }


/* Calculate print frequency based on DialInStepTime / TimeStep, tuned by a
dimensionless parameter in the setting file and rounded to the nearest integer.
This rounding should work fine as long as the PrintFrequency isn't ridiculously
huge. To avoid this we require that PrintFrequency is not in [0.0,1.0). If
PrintFrequency is very large (more likely), so that InversePrintRate rounds to 0,
we set InversePrintRate to 1 to get a printout after every time step. This can
be a useful thing to do sometimes. If PrintFrequency is set to a negative value
we instead set InversePrintRate to -1, which means no output files are regularly
written at all.
*/
    if (settings.PrintFrequency < 0.0) {
        settings.InversePrintRate = -1;
    } else {
        if (settings.PrintFrequency < 1.0) {
            logStream
                    << "Error: To avoid potential divide-by-zero problems we do NOT allow settings.PrintFrequency to lie in [0.0, 1.0). This is overkill somewhat, and could be relaxed if need be."
                    << std::endl;
            return -1;
        } else {
            settings.InversePrintRate = lround(settings.DialInStepTime / (settings.TimeStep * settings.PrintFrequency));
            if (settings.InversePrintRate == 0) {
                settings.InversePrintRate = 1;
            }
        }
    }



// Print total load force.
    int numLoadedNodes = 0;
    for (int i = 0; i < settings.NumNodes; ++i) {
        if (nodes[i].isLoadForceEnabled) {
            numLoadedNodes += 1;
        }
    }
    logStream.open();
    logStream << "Total load force applied = " <<
              numLoadedNodes * settings.LoadStrength * settings.ShearModulus * settings.ApproxMinInitElemSize *
              settings.SheetThickness << std::endl;
    logStream.close();


// Calculate and store characteristic force, energy, and energy density scales.
    settings.charForceScale = settings.ShearModulus * settings.ApproxMinInitElemSize * settings.SheetThickness;
    settings.charStretchEnergyDensityScale = settings.ShearModulus * settings.SheetThickness;
// settings.charBendEnergyDensityScale = settings.ShearModulus * settings.SheetThickness * settings.SheetThickness * settings.SheetThickness / (settings.SampleCharLength * settings.SampleCharLength);
    double totInitArea = 0;
    std::vector<double> initAreas(settings.NumTriangles);
    for (int i = 0; i < settings.NumTriangles; ++i) {
        initAreas[i] = triangles[i].initArea;
    }
    totInitArea = kahanSum(initAreas);
    settings.charStretchEnergyScale = settings.charStretchEnergyDensityScale * totInitArea;
// settings.charBendEnergyScale = settings.charBendEnergyDensityScale * totInitArea;



/* Variables used to control the programmed properties' time variation.
This is not necessarily directly physical, especially for a fast optical
activation instead of slow heating; the idea here is to gradually `dial in'
the programmed (energetically favoured) metric and second fundamental form.
This prevents anything too explosive happening the code, and can be used to
make the simulation quasistatic (which can help avoid undesired isometries
for example), though it may not be necessary.
The range between 0, and 1 is split into a sequence of discrete values based
on the chosen settings.DialInResolution. The dial-in factor will then
evolve linearly from one such value to the next, over a chosen
settings.DialInStepTime. After each such step, the value is held constant
until equilibrium is reached (to within the chosen tolerance). In the
'waiting' stage, the check for equilibrium occurs at a rate determined by
the settings.TimeBetweenEquilChecks. The dialling in process occurs
between each pair of the sequence of programmed tensors, if more than one is
given in the input file.
*/
    double timeSinceLastEquilCheck = -98765.0; // Reset to zero every time equil is checked, or every time a new DialInFactor value is reached.
    double timeSinceCurrDiallingInPhaseStarted = -7654321.0; // Reset to zero every time equil is reached and a new `dialling in' phase starts.
    SimulationStatus status; // Indicates whether the simulation is currently in a 'dialling in' phase, or is not dialling and is instead waiting for equilibrium, or whether equilibrium has been reached but the next dialling phase has not yet begun.
    std::size_t DialInFactorCounter = 123456; // Keeps track of which DialInFactor value in DialInFactorValuesToHoldAt was last dialled *from* (not held at, which is the next value along in the list)
    double currDialInFactor = -54321.0; // Keeps track of current value

    std::vector<double> DialInFactorValuesToHoldAt;
    double tempDialInFactor = 0.0;
    if ((settings.DialInResolution > 0) && settings.DialInStepTime >= 0) {
        while (tempDialInFactor < 1.0) {
            DialInFactorValuesToHoldAt.push_back(tempDialInFactor);
            tempDialInFactor += settings.DialInResolution;
        }
        DialInFactorValuesToHoldAt.push_back(1.0);
    } else {
        logStream.open();
        logStream << "settings.DialInResolution and settings.DialInStepTime must "
                     "both be >0 and >=0 respectively, and they aren't currently. Aborting."
                  << std::endl;
        logStream.close();
        return -1;
    }

/* Create std::vectors (with one element for each triangle) that will be passed
by reference to functions calculating the Gauss and mean curvatures, and
stretching and bending energies and energy densities. */
    std::vector<double> gaussCurvatures(settings.NumTriangles,
                                        -98765.4321); //Recognisable initialisation for debugging.
    std::vector<double> meanCurvatures(settings.NumTriangles, -98765.4321);
    std::vector<double> stretchEnergies(settings.NumTriangles, -98765.4321);
    std::vector<double> bendEnergies(settings.NumTriangles, -98765.4321);
    std::vector<double> stretchEnergyDensities(settings.NumTriangles, -98765.4321);
    std::vector<double> bendEnergyDensities(settings.NumTriangles, -98765.4321);
    std::vector<double> kineticEnergies(settings.NumNodes, -98765.4321);
    std::vector<double> strainMeasures(settings.NumTriangles, -98765.4321);
    std::vector <Eigen::Vector2d> cauchyStressEigenvals(settings.NumTriangles);
    std::vector <Eigen::Matrix<double, 3, 2>> cauchyStressEigenvecs(settings.NumTriangles);


/* Create further std::vectors to store angle deficits if specified in settings
file.*/
    std::vector<double> angleDeficits(settings.NumNodes, -98765.4321);
//std::vector<double> interiorNodeAngleDeficits(settings.NumNodes - settings.numBoundaryNodes, -98765.4321);
//std::vector<double> boundaryNodeAngleDeficits(settings.numBoundaryNodes, -98765.4321);
    std::cout << "CHECK THIS" << std::endl;
    std::vector<double> interiorNodeAngleDeficits(settings.NumNodes, -98765.4321);
    std::vector<double> boundaryNodeAngleDeficits(settings.NumNodes, -98765.4321);

/*
std::cout<< "STORING NODE REF POSITIONS TO HELP MAKE EXACT CONE." << std::endl;
std::vector< Eigen::Vector3d > nodeRefPosits(settings.NumNodes);
for(int n = 0; n < settings.NumNodes; ++n){
    nodeRefPosits[n] = nodes[n].pos;
}*/

/* Perturb node positions with small random noise if desired, to allow 'breaking
away' from flat plane initial condition, for example. This will have no effect
if an ansatz is used. */
    if (settings.isPerturbationOfInitialPositionsEnabled) {
        perturbInitialPositionsWithRandomNoise(nodes, settings);
    }

////////////////////////////////////////////////////////////////////////////////

/* Initialise time and step counter. Set logStream to more useful
format and precision for e.g. energy printouts. */
    double time = 0;
    int stepCount = 0;
    logStream << std::scientific << std::setprecision(8);

    std::ofstream forceDistFile;


/* Loop over the sequence of programmed tensors, dialling-in and waiting for
equilibrium between each pair in the sequence. This loop is redundant in most use
cases for this code, where only a single set of programmed tensors is supplied.*/
    for (std::size_t progTensorSequenceCounter = progTensorSequenceCounterToStartFrom;
         progTensorSequenceCounter <= sequenceOf_progMetricInfo.size() - 2; ++progTensorSequenceCounter) {

        /* A third input file specifying an ansatz for the node positions may have
        been read in if it was given as a command line argument. This could for
        example correspond to an output file from this code, to carry on where some
        previous simulation left off. In this case, the nodes are moved to their
        ansatz positions here, and other relevant variables are set up. */
        if (ansatz_data_file_name_str != "no_ansatz_file" &&
            progTensorSequenceCounter == progTensorSequenceCounterToStartFrom) {

            for (int i = 0; i < settings.NumNodes; ++i) {
                nodes[i].pos = nodeAnsatzPositions[i];
            }

            timeSinceLastEquilCheck = 0.0;
            status = dialling;

            if (!settings.isDialingFromAnsatzEnabled) {

                currDialInFactor = dialInFactorToStartFrom;
                DialInFactorCounter = 0;
                for (std::size_t i = 0; i < DialInFactorValuesToHoldAt.size(); ++i) {
                    if (DialInFactorValuesToHoldAt[i] < currDialInFactor) {
                        DialInFactorCounter = i;
                    }
                }

                timeSinceCurrDiallingInPhaseStarted =
                        ((currDialInFactor - DialInFactorValuesToHoldAt[DialInFactorCounter]) /
                         (DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                          - DialInFactorValuesToHoldAt[DialInFactorCounter])) * settings.DialInStepTime;

                /* Here we implement another feature, in which the progMetric and
                progSecFF that we would normally be dialling FROM at this point, are
                just set to the values they would be dialling TO, so that no dialling
                of these quantities happens between the current (starting) entry in
                the progTensor sequence and the next. progTau however IS dialled in as
                normal. Upon reaching the next entry in the sequence, normal behaviour
                resumes so there will be a phase of waiting for equilibrium,
                followed potentially by normal dialling in of further progTensors in
                the sequence, if there are any more in the sequence.
                This feature's intended use case is where you want to dial in some
                large progTaus on the boundary triangles, to reduce stretch-bend
                tradeoff boundary effects, but you don't want to do any dialling of
                the progMetric or progSecFF because you think the shape is already
                close to these in its achieved metric and secFF, probably because
                you already did a simulation with them (but no large progTau), and
                are using the output of that run as an ansatz. This does mean that
                you'll want to CHANGE THE currDialInFactor IN THE ANSATZ FILE TO = 0.
                Otherwise the whole dialling in phase will be skipped entirely, and
                the progTaus will just jump to their new large values, which may well
                produce violent behaviour. You might be tempted to just use the
                isDialingFromAnsatzEnabled feature instead of this one, since
                then the new progTaus are jumped to, but the progMetric and secFF are
                dialled in from the initial actual metric and actual secFF values,
                ensuring gentle dynamics. This is probably non-ideal though, because
                the progMetric and progSecFF along that dialling path may then
                deviate significantly from the form that you want to arrive at, and
                that you thought you started close to, pumping unnecessary energy in
                to the system and wasting time.
                There may be other uses for this feature too.*/
                if (settings.ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF) {
                    for (int i = 0; i < settings.NumTriangles; ++i) {
                        sequenceOf_ProgSecFFs[progTensorSequenceCounterToStartFrom][i] = sequenceOf_ProgSecFFs[
                                progTensorSequenceCounterToStartFrom + 1][i];
                        sequenceOf_progMetricInfo[progTensorSequenceCounterToStartFrom][i] = sequenceOf_progMetricInfo[
                                progTensorSequenceCounterToStartFrom + 1][i];
                    }
                }
            }

                /* The isDialingFromAnsatzEnabled setting is very handy if you can
                only make an approximate guess at a solution, and want to start from
                there and see where things go without explosive dynamics.
                If settings.isDialingFromAnsatzEnabled == true, we instead
                arrange to ignore the currDialInFactor, DialInFactorCounter and
                timeSinceCurrDiallingInPhaseStarted variables that are implied by the
                the ansatz data file, and instead just start a new dialling in phase
                from whatever the *calculated* metric and secFF are in the ansatz state.
                Tau is not a calculable geometric quantity so we just set it to the
                programmed value immediately, which will not cause explosive dynamics
                even if e.g. the ansatz shape came from a simulation with a very different
                programmed tau, because the stretch and bend energies both still start
                at zero etc.
                */
            else {
                logStream.open();
                logStream
                        << "\nsettings.isDialingFromAnsatzEnabled == true, FIRST DIALLING IN PHASE WILL BE FROM ANSATZ STATE."
                        << std::endl;
                logStream.close();

                currDialInFactor = 0.0;
                DialInFactorCounter = 0;
                timeSinceCurrDiallingInPhaseStarted = 0.0;

                // Calculate all necessary geometry for the ansatz state.
                calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, waitingForEquilibrium, currDialInFactor,
                                                              progTensorSequenceCounter, sequenceOf_progMetricInfo,
                                                              sequenceOf_InvProgMetrics, sequenceOf_ProgTaus,
                                                              sequenceOf_ProgSecFFs, settings);
                calcSecFFsAndRelatedQuantities(triangles, settings);


                // Temp LU decomp of triangle metric, used to check invertibility.
                Eigen::FullPivLU<Eigen::Matrix<double, 2, 2> > tempMetricDecomp;


                /* Alter sequenceOf_InvProgMetrics[progTensorSequenceCounterToStartFrom] and similar to change where
                the programmed quantities are dialling from.*/

                for (int i = 0; i < settings.NumTriangles; ++i) {

                    // Test for invertibility of metric before taking inverse.
                    try {
                        tempMetricDecomp.compute(
                                (triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse());
                        if (!tempMetricDecomp.isInvertible()) {
                            throw std::runtime_error(
                                    "At least one triangle had a non-invertible metric in the ansatz state. \n"
                                    "This should not occur in a reasonable mesh. Aborting.");
                        } else {
                            sequenceOf_InvProgMetrics[progTensorSequenceCounterToStartFrom][i] = (
                                    triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                            if (settings.isDialingDisabled) {
                                sequenceOf_InvProgMetrics[progTensorSequenceCounterToStartFrom + 1][i] = (
                                        triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                            }
                        }
                    }
                    catch (const std::runtime_error &nonInvertibleMetric_error) {
                        logStream.open();
                        logStream << nonInvertibleMetric_error.what() << std::endl;
                        logStream.close();
                        return -1;
                    }

                    // Other quantities require no inverse, so are easier.
                    sequenceOf_ProgTaus[progTensorSequenceCounterToStartFrom][i] = sequenceOf_ProgTaus[
                            progTensorSequenceCounterToStartFrom + 1][i];
                    sequenceOf_ProgSecFFs[progTensorSequenceCounterToStartFrom][i] = triangles[i].secFF;
                    if (settings.isDialingDisabled) {
                        sequenceOf_ProgTaus[progTensorSequenceCounterToStartFrom + 1][i] = sequenceOf_ProgTaus[
                                progTensorSequenceCounterToStartFrom + 1][i];
                        sequenceOf_ProgSecFFs[progTensorSequenceCounterToStartFrom + 1][i] = triangles[i].secFF;
                    }
                }
            }

            /////////////////////////////////////////////////////////////////////

            // If set in settings file, clamp all nodes within a certain vertical distance
            // above the node with the lowest initial z value (only works with ansazt clearly).
            if (settings.ThicknessesAboveLowestNodeToClampUpTo > 0) {
                double minNodeZCoord = nodes[0].pos(2);
                for (int n = 0; n < settings.NumNodes; ++n) {
                    if (minNodeZCoord > nodes[n].pos(2)) {
                        minNodeZCoord = nodes[n].pos(2);
                    }
                }

                for (int n = 0; n < settings.NumNodes; ++n) {
                    if (nodes[n].pos(2) - minNodeZCoord <
                        settings.ThicknessesAboveLowestNodeToClampUpTo * settings.SheetThickness) {
                        nodes[n].isClamped = true;
                    }
                }
            }

            ///////////////////////////////////////////////////////////////////
        }

            // If no ansatz is being used:
        else {
            // Reset the variables that control each dialling/waiting process.
            timeSinceLastEquilCheck = 0.0;
            status = dialling;
            currDialInFactor = 0.0;
            DialInFactorCounter = 0;
            timeSinceCurrDiallingInPhaseStarted = 0.0;
        }


        /* Handle cases where settings.DialInStepTime is zero or very small. This is
        a slightly hacky way to ensure that no dialling actually occurs, and the
        simulation jumps to a status = waitingForEquilibrium state. The magic number is
        just chosen to be recognisable for debugging.*/
        if (settings.DialInStepTime < settings.TimeStep) {
            settings.DialInStepTime = 0.0;
            timeSinceCurrDiallingInPhaseStarted = 1.23456789;
        }

        ///////////////////////////////////////////////////////////////////////////

        // Begin dynamical evolution of node positions and velocities.
        logStream.open();
        logStream << "\nBeginning dynamical evolution.\n" << std::endl;
        logStream.close();

        logStream.open();
        logStream << "\nCREATING VECTOR TO STORE UNSTRESSED CONE NODE POSITIONS.\n" << std::endl;
        logStream.close();
        std::vector <Eigen::Vector3d> nodeUnstressedConePosits(settings.NumNodes);
        double s1 = 12345678.9;
        //double sTest = 98765432.1;

        std::vector < std::vector < std::pair <
        int, int>>> correspondingTrianglesForNodes = getCorrespondingTrianglesForNodes(
                triangles, nodes);

        while (DialInFactorCounter <= DialInFactorValuesToHoldAt.size() - 2) {
            int highestNode = -99;
            int lowestNode = -99;
            if (stepCount == 0) {
                settings.initSlideZCoord_lower = nodes[0].pos(2);
                settings.initSlideZCoord_upper = nodes[0].pos(2);
                for (int n = 0; n < settings.NumNodes; ++n) {
                    if (settings.initSlideZCoord_lower > nodes[n].pos(2)) {
                        settings.initSlideZCoord_lower = nodes[n].pos(2);
                        lowestNode = n;
                    }
                    if (settings.initSlideZCoord_upper < nodes[n].pos(2)) {
                        settings.initSlideZCoord_upper = nodes[n].pos(2);
                        highestNode = n;
                    }
                    // Can instead choose initial slide separation directly from settings file, to help avoid 'jumping' when starting a
                    // squashing run part way through from a previous run's output.
                    if (settings.SpecifyInitSlideZCoord_upper > -99.0) {
                        settings.initSlideZCoord_upper = settings.SpecifyInitSlideZCoord_upper;
                    }
                    if (settings.SpecifyInitSlideZCoord_lower > -99.0) {
                        settings.initSlideZCoord_lower = settings.SpecifyInitSlideZCoord_lower;
                    }

                    // To do constant-weight slide instead
                    if (settings.constSlideWeightFac > 0) {
                        settings.initSlideZCoord_upper = settings.initSlideZCoord_lower + settings.SpacerHeight *
                                                                                          settings.SampleCharLength;// This used to be settings.SpacerHeight * settings.SampleCharLength instead, which is wrong
                    }

                    settings.upperSlideDisplacement = 0.0;
                    settings.upperSlideVel = 0.0;
                    settings.slideJustReachedEquil = 0;
                    if (settings.isControlledForceEnabled) {
                        settings.upperSlideWeight = settings.InitialSlideWeightForCtrldForceInUnitsOfMuTsq *
                                                    (settings.ShearModulus * settings.SheetThickness *
                                                     settings.SheetThickness);
                    }
                }
                // Uncomment the following line for single ridge experiment with point(ish) load to applied tip.
                //settings.initSlideZCoord_upper = settings.initSlideZCoord_lower + 2.0 * settings.SheetThickness; // to apply point(ish) load to tip

                // Readjust for the case of glass cones instead of glass slides. The
                // slide Z coords correspond to the tips of the glass cones.
                if (settings.GlassCones) {
                    settings.ConeAngle = 1.02327019;
                    settings.initSlideZCoord_upper += -tan(settings.ConeAngle) *
                                                      sqrt(nodes[highestNode].pos(0) * nodes[highestNode].pos(0) +
                                                           nodes[highestNode].pos(1) * nodes[highestNode].pos(1));
                    settings.initSlideZCoord_lower += -tan(settings.ConeAngle) *
                                                      sqrt(nodes[lowestNode].pos(0) * nodes[lowestNode].pos(0) +
                                                           nodes[lowestNode].pos(1) * nodes[lowestNode].pos(1));
                    logStream.open();
                    logStream << "USING TWO GLASS CONES FOR SQUASHING." << std::endl;
                    logStream.close();

                    // Hijack nodes[i].isOnBoundary to instead label nodes whose
                    // forces we will modify to kill any components not tangential to
                    // a perfect cone base state. We do this to nodes within an intermediate
                    // distance vertically from the ends of the cone.
                    double intermLengthScaleUpper = 2.0 * sqrt(settings.SheetThickness * 0.18);
                    double intermLengthScaleLower = 2.0 * sqrt(settings.SheetThickness * 1.8);
                    logStream.open();
                    logStream
                            << "Apply normal-force-killer within the following vertical distances of the top and bottom: "
                            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
                    logStream.close();
                    for (int n = 0; n < settings.NumNodes; ++n) {
                        if (((nodes[highestNode].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
                            ((nodes[n].pos(2) - nodes[lowestNode].pos(2)) < intermLengthScaleLower)) {
                            nodes[n].isOnBoundary = true;
                        }
                    }
                }

                // Instead set up imposed-Seide-deformations idea.
                if (settings.isSeideDeformationsEnabled) {
                    settings.lambda = 0.9;
                    settings.ConeAngle = asin(pow(settings.lambda, 1.5));
                    logStream.open();
                    logStream << "IMPOSING SEIDE DEFORMATIONS." << std::endl;
                    logStream.close();

                    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
                    // positions we will force to be those of Seide's setup (found
                    // by findng the membrane stress solution for his base state, and
                    // integrating the strains etc).
                    double intermLengthScaleUpper = 3.0 * sqrt(settings.SheetThickness * 0.18);
                    double intermLengthScaleLower = 3.0 * sqrt(settings.SheetThickness * 1.8);
                    logStream.open();
                    logStream
                            << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
                            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
                    logStream.close();
                    for (int n = 0; n < settings.NumNodes; ++n) {
                        if (((nodes[highestNode].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
                            ((nodes[n].pos(2) - nodes[lowestNode].pos(2)) < intermLengthScaleLower)) {
                            nodes[n].isSeideDisplacementEnabled = true;
                        }
                    }

                    // STORE PERFECT CONE ANSATZ
                    for (int n = 0; n < settings.NumNodes; ++n) {
                        nodeUnstressedConePosits[n] = nodes[n].pos;
                    }
                    s1 = sqrt(nodes[highestNode].pos(0) * nodes[highestNode].pos(0) +
                              nodes[highestNode].pos(1) * nodes[highestNode].pos(1)) / sin(settings.ConeAngle);

                    Eigen::Vector3d testTriCurrCentroid;
                    testTriCurrCentroid = (nodes[triangles[settings.testTriangle].vertexLabels(0)].pos +
                                           nodes[triangles[settings.testTriangle].vertexLabels(1)].pos +
                                           nodes[triangles[settings.testTriangle].vertexLabels(2)].pos) / 3;
                    //sTest = sqrt(testTriCurrCentroid(0)*testTriCurrCentroid(0) + testTriCurrCentroid(1)*testTriCurrCentroid(1)) / sin(settings.ConeAngle);
                }
            }


            if (!settings.isControlledForceEnabled) {
                settings.upperSlideDisplacement =
                        time * settings.slideSpeedPrefactor * settings.SampleCharLength / bendingLongTime;
            } else { // settings.isControlledForceEnabled == true instead
                //settings.upperSlideWeight = (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) * (time * settings.slideSpeedPrefactor / bendingLongTime);
                settings.slideDampingParam =
                        0.4 * settings.ShearModulus * settings.SheetThickness * settings.SheetThickness /
                        (settings.slideSpeedPrefactor * settings.SampleCharLength / bendingLongTime);
                if (fabs(settings.upperTotSlideForce + settings.upperSlideWeight) /
                    (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) <
                    settings.totalSlideForceToMuTSqRatioEquilThreshold
                    && settings.constSlideWeightFac < 0) {
                    if (timeSinceLastEquilCheck > settings.TimeBetweenEquilChecks) {
                        if (equilibriumCheck(nodes, triangles, settings, logStream) == equilibriumReached) {
                            settings.slideJustReachedEquil = 1;
                            settings.upperSlideWeight +=
                                    settings.slideWeightDialSpeedFac * (settings.TimeStep / bendingLongTime) *
                                    (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness);
                        }
                        timeSinceLastEquilCheck = 0.0;
                    }
                }

                // To do constant-weight slide instead
                if (settings.constSlideWeightFac > 0) {
                    settings.upperSlideWeight =
                            settings.constSlideWeightFac * settings.ShearModulus * settings.SheetThickness *
                            settings.SheetThickness;
                }
            }
            settings.currSlideZCoord_upper = settings.initSlideZCoord_upper - settings.upperSlideDisplacement;
            std::pair<double, double> upperAndLowerTotSlideForces;

            // Impose Seide deformations.
            if (settings.isSeideDeformationsEnabled) {
                double pInit = 3.0 * (settings.ShearModulus * settings.SheetThickness *
                                      settings.SheetThickness); // So we don't have to start all the way from p=0, chosen based on previous sims.
                settings.p = pInit + time * settings.pSpeedPrefactor * settings.ShearModulus * settings.SheetThickness *
                                     settings.SheetThickness / bendingLongTime;

                for (int n = 0; n < settings.NumNodes; ++n) {
                    if (nodes[n].isSeideDisplacementEnabled || stepCount == 0) {

                        double tCone = settings.SheetThickness / sqrt(settings.lambda);

                        double polarAng = atan2(nodeUnstressedConePosits[n](1), nodeUnstressedConePosits[n](0));
                        double s = sqrt(nodeUnstressedConePosits[n](0) * nodeUnstressedConePosits[n](0) +
                                        nodeUnstressedConePosits[n](1) * nodeUnstressedConePosits[n](1)) /
                                   sin(settings.ConeAngle);
                        Eigen::Vector3d uHat;
                        uHat << sin(settings.ConeAngle) * cos(polarAng), sin(settings.ConeAngle) * sin(polarAng), -cos(
                                settings.ConeAngle);
                        Eigen::Vector3d wHat;
                        wHat << -cos(settings.ConeAngle) * cos(polarAng), -cos(settings.ConeAngle) *
                                                                          sin(polarAng), -sin(settings.ConeAngle);
                        double u = -settings.p * log(s / s1) /
                                   (2.0 * PI * settings.YoungsModulus * tCone * sin(settings.ConeAngle) *
                                    cos(settings.ConeAngle));
                        double w = -settings.p * (log(s / s1) + settings.PoissonRatio) /
                                   (2.0 * PI * settings.YoungsModulus * tCone * cos(settings.ConeAngle) *
                                    cos(settings.ConeAngle));
                        nodes[n].pos = nodeUnstressedConePosits[n] + u * uHat + w * wHat;
                    }
                }
            }


            /* The forces are set each timestep using a += procedure, and therefore
            must be set to zero each time before this is done.*/
            zeroForces(nodes);

            // Check if still in dialling in phase or whether it is time to wait for
            //equilibrium.
            if (timeSinceCurrDiallingInPhaseStarted >= settings.DialInStepTime && status == dialling) {

                /* If the current dialling phase has indeed finished, set
                currDialInFactor to the value that was being dialled up to in
                that phase, and set the programmed tensors accordingly.*/
                currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter + 1];
                calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, status, currDialInFactor,
                                                              progTensorSequenceCounter, sequenceOf_progMetricInfo,
                                                              sequenceOf_InvProgMetrics, sequenceOf_ProgTaus,
                                                              sequenceOf_ProgSecFFs, settings);
                status = waitingForEquilibrium;
                // NB an EquilCheck has not actually just occurred, but this has the
                //desired effect of ensuring that each DialInFactor value is held
                //for at least one TimeBetweenEquilChecks.
                timeSinceLastEquilCheck = 0.0;

                settings.NumDampFactor =
                        settings.DampingPrefactor2 * dampingScale; // Set damping factor to waiting phase value.

                logStream.open();
                logStream << "Reached Dial-In Factor of " << DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                          << ", now waiting for equilibrium" << std::endl;
                logStream.close();
            }

            /* If not waiting for equilibrium, set the current value of the Dial-In
            Factor, based on linear dialling in between the previously calculated
            'checkpoint' values. */
            if (status == dialling) {

                currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter] +
                                   (DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                                    - DialInFactorValuesToHoldAt[DialInFactorCounter]) *
                                   timeSinceCurrDiallingInPhaseStarted / settings.DialInStepTime;
            }

            /* Calculate sides, areas...etc of each triangle, as well as the current
            dialled-in programmed (inverse) metric and second fundamental form.*/

            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
            calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, status, currDialInFactor,
                                                          progTensorSequenceCounter, sequenceOf_progMetricInfo,
                                                          sequenceOf_InvProgMetrics, sequenceOf_ProgTaus,
                                                          sequenceOf_ProgSecFFs, settings);

            std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();

            /* Calculate secFF estimates for triangles, and related quantities such
            as the derivative of the bending energy wrt the secFF components.*/
            calcSecFFsAndRelatedQuantities(triangles, settings);
            std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();

            // Calculate current strain and bending force on each node.
            calcDeformationForces(nodes, triangles, settings, correspondingTrianglesForNodes);
            std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();

            /* Add force contributions from e.g. damping, loads, 'prod' perturbation, and
             account for BCs e.g. clamping. */
            upperAndLowerTotSlideForces = calcNonDeformationForces_and_ImposeBCS(nodes, time, settings);
            std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
            settings.upperTotSlideForce = upperAndLowerTotSlideForces.first;


            // Hijack upperSlideDisplacement to hold end displacement in pulling experiment.
            /*
            if( settings.LoadStrength > 0.0 ){
                upperSlideDisplacement = nodes[0].pos(0);
                for(int n = 0; n < settings.NumNodes; ++n){
                    if(upperSlideDisplacement < nodes[n].pos(0)){
                        upperSlideDisplacement = nodes[n].pos(0);
                    }
                }
            }*/

            /* Write output data regularly.
            Can be switched off with settings.PrintFrequency < 0.0.
            Doing the write-out at this point in the loop means the node positions
            and the triangle geometry data match in the output, which is desirable!*/
            if ((stepCount % settings.InversePrintRate == 0 && settings.InversePrintRate > 0) ||
                settings.slideJustReachedEquil == 1) {
                try {
                    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
                    if (settings.isEnergyDensitiesPrinted) {
                        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
                    }

                    // Here we hack the angleDeficits to instead tell us whether a node has a Seide displacement imposed or not.
                    for (int n = 0; n < settings.NumNodes; ++n) {
                        if (nodes[n].isSeideDisplacementEnabled) {
                            angleDeficits[n] = 1;
                        } else {
                            angleDeficits[n] = 0;
                        }
                    }
                    writeVTKDataOutput(nodes, triangles, stepCount, time, currDialInFactor, progTensorSequenceCounter,
                                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
                }
                catch (const std::runtime_error &writing_error) {
                    logStream << writing_error.what() << std::endl;
                    return -1;
                }
                logStream.open();
                logStream << std::fixed << "Wrote VTK output at " << getRealTime() << ", stepCount = " << stepCount
                          << ", simulation time = " << time + settings.TimeStep << ", current dial-in factor = "
                          << currDialInFactor << std::scientific;
                logStream << ", last step's execution time "
                          << std::chrono::duration_cast<std::chrono::microseconds>(end4 - begin).count() << " us"
                          << std::endl;

//                logStream << "\tcalcTriangleGeometries execution time "
//                          << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin).count() << " us"
//                          << std::endl;
//                logStream << "\tcalcSecFFsAndRelatedQuantities execution time "
//                          << std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1).count() << " us"
//                          << std::endl;
//                logStream << "\tcalcDeformationForces execution time "
//                          << std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2).count() << " us"
//                          << std::endl;
//                logStream << "\tcalcNonDeformationForces_and_ImposeBCS execution time "
//                          << std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3).count() << " us"
//                          << std::endl;

                logStream.close();

                forceDistFile.open(outputDirName + "/force_displacement_vals.txt", std::ofstream::app);
                forceDistFile << std::scientific << std::setprecision(14)
                              << upperAndLowerTotSlideForces.first /
                                 (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) << " "
                              << upperAndLowerTotSlideForces.second /
                                 (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) << " "
                              << settings.upperSlideDisplacement / settings.SampleCharLength << " "
                              << settings.currSlideZCoord_upper << " "
                              << settings.initSlideZCoord_lower << " "
                              << settings.slideJustReachedEquil << " "
                              << settings.upperSlideWeight /
                                 (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness)

                              //<< settings.p / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " "
                              //<< cauchyStressEigenvals.at(settings.testTriangle)(0) * sTest / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " "
                              //<< cauchyStressEigenvals.at(settings.testTriangle)(1) * sTest / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " " << std::endl;
                              << std::endl;
                forceDistFile.close();
                settings.slideJustReachedEquil = 0;
            }


            if (!settings.isControlledForceEnabled) {


                /* If last check for equilibrium or last reaching of a new DialInFactor
                 value was more than TimeBetweenEquilChecks ago, check for equilibrium.
                 Also print total stretching and bending energies.*/
                if (timeSinceLastEquilCheck > settings.TimeBetweenEquilChecks && status == waitingForEquilibrium) {
                    logStream.open();
                    logStream << "Checking for equilibrium at " << getRealTime() << ", stepCount = " << stepCount
                              << ", simulation time = " << time << ", current dial-in factor = " << currDialInFactor
                              << std::endl;
                    logStream.close();

                    status = equilibriumCheck(nodes, triangles, settings, logStream);
                    timeSinceLastEquilCheck = 0.0;

                    calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                            stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                            cauchyStressEigenvals, cauchyStressEigenvecs, settings);
                    double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.charStretchEnergyScale;
                    double nonDimBendEnergy = kahanSum(bendEnergies) / settings.charStretchEnergyScale;
                    double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.charStretchEnergyScale;

                    logStream.open();
                    logStream << "Non-dimensionalised energies:" << std::endl;
                    logStream << "\t total \t\t" << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy
                              << std::endl;
                    logStream << "\t stretch \t" << nonDimStretchEnergy << std::endl;
                    logStream << "\t bend \t\t" << nonDimBendEnergy << std::endl;
                    logStream << "\t kinetic \t" << nonDimKineticEnergy << std::endl;
                    logStream.close();
                }

                /* If equilibrium reached, write output data to file, and move to next
                'dialling in' phase. */
                if (status == equilibriumReached) {
                    try {
                        calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                                       interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
                        if (settings.isEnergyDensitiesPrinted) {
                            calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                                    stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                                    cauchyStressEigenvals, cauchyStressEigenvecs, settings);
                        }
                        writeVTKDataOutput(nodes, triangles, stepCount, time, currDialInFactor,
                                           progTensorSequenceCounter, gaussCurvatures, meanCurvatures, angleDeficits,
                                           interiorNodeAngleDeficits, boundaryNodeAngleDeficits, stretchEnergyDensities,
                                           bendEnergyDensities, stretchEnergies, bendEnergies, kineticEnergies,
                                           strainMeasures, cauchyStressEigenvals, cauchyStressEigenvecs, settings,
                                           outputDirName);
                    }
                    catch (const std::runtime_error &writing_error) {
                        logStream.open();
                        logStream << writing_error.what() << std::endl;
                        logStream.close();
                        return -1;
                    }
                    logStream.open();
                    logStream << "Equilibrium reached. Wrote VTK output at " << getRealTime() << ", stepCount = "
                              << stepCount << ", simulation time = " << time << ", current dial-in factor = "
                              << currDialInFactor << std::endl;
                    status = dialling;
                    timeSinceCurrDiallingInPhaseStarted = 0.0;
                    DialInFactorCounter += 1;
                    settings.NumDampFactor = settings.DampingPrefactor1 *
                                             dampingScale; // Set damping factor back to dialling phase value.
                    if (DialInFactorCounter <= DialInFactorValuesToHoldAt.size() - 2) {
                        logStream << "New dialling in phase beginning, from value "
                                  << DialInFactorValuesToHoldAt[DialInFactorCounter] << ", to "
                                  << DialInFactorValuesToHoldAt[DialInFactorCounter + 1] << std::endl;
                    }
                    logStream.close();
                }


            }


            /* Advance node positions and velocities using dynamical solver. Error
            caught here if any node force is suspiciously high, indicating probable
            'blowing, up', and the code aborts in that case. */
            try { advanceDynamics(nodes, triangles, settings, logStream); }
            catch (const std::runtime_error &error) {
                logStream.open();
                logStream << "At " << getRealTime() << ", stepCount = " << stepCount << ", time = " << time
                          << ", dial-in factor = " << currDialInFactor
                          << " a force was suspiciously large, there is probably a problem. Writing VTK output and then aborting."
                          << std::endl;
                logStream.close();
                try {
                    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
                    if (settings.isEnergyDensitiesPrinted) {
                        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
                    }
                    writeVTKDataOutput(nodes, triangles, stepCount, time, currDialInFactor, progTensorSequenceCounter,
                                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
                }
                catch (const std::runtime_error &writing_error) {
                    logStream.open();
                    logStream << writing_error.what() << std::endl;
                    logStream.close();
                    return -1;
                }
                return -1;
            }

            // Advance times and step counter.
            time += settings.TimeStep;
            timeSinceLastEquilCheck += settings.TimeStep;
            timeSinceCurrDiallingInPhaseStarted += settings.TimeStep;
            stepCount += 1;
        }

        if (sequenceOf_progMetricInfo.size() > 2 && progTensorSequenceCounter < sequenceOf_progMetricInfo.size() - 2) {
            logStream.open();
            logStream << "Moving on to next set of programmed tensors in sequence." << std::endl;
            logStream.close();
        }
    }

// Print some helpful final things.
    logStream.open();
    logStream << "Reached simulation time = " << time << " using " << stepCount << " time steps" << std::endl;
    logStream.close();

    return 0;
}
