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

#define _USE_MATH_DEFINES

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


//#include "CustomOutStreamClass.hpp"
#include "Node.hpp"
#include "Triangle.hpp"
#include "Edge.hpp"

#include "functions/getRealTime.hpp"
#include "functions/extract_Just_Filename.hpp"
#include "initialDirAndFileHandling.hpp"
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
#include "calculations/calcDeformationForces.hpp"
#include "advanceDynamics.hpp"
#include "writeVTKDataOutput.hpp"
#include "EquilibriumCheck.hpp"
#include "calculations/calcCurvatures.hpp"
#include "calculations/calcEnergiesAndStresses.hpp"
#include "functions/kahanSum.hpp"
#include "simulation.h"

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
    Simulation simulation(argc, argv);
    return simulation.run_simulation();
}


//int legacy_main(int argc, char *argv[]) {
///* Create a string to contain some things that will later be written to the log
//file, as long as all goes to plan and the code gets that far. The first thing to
//add to this string is a record of what command line arguments where given, and
//at what date and time the run began.*/
//
//    std::string initWriteToLogStr;
//
//    initWriteToLogStr += "Simulation run began at: " + getRealTime() + "\n";
//    initWriteToLogStr += "The command that was run was:\n";
//    for (int i = 0; i < argc; ++i) {
//        initWriteToLogStr += argv[i];
//        initWriteToLogStr += " ";
//    }
//    initWriteToLogStr += "\n\n";
//
//
///* Check settings and data files given, with ansatz file optional, and get names
//from command line. */
//    if ((argc != 3) && (argc != 4)) {
//        std::cerr << "Error: Please provide at least one settings file path "
//                     "followed by one data file path as command line arguments.\n"
//                     "A path to a file with ansatz node coordinates can be given as a third "
//                     "argument if desired. All paths should be relative to the current working"
//                     "directory, i.e. they should start with ./" << std::endl;
//        return -1;
//    }
//    const char *settings_file_name = argv[1];
//    std::string settings_file_name_str(settings_file_name);
//    const char *initial_data_file_name = argv[2];
//    std::string initial_data_file_name_str(initial_data_file_name);
//    std::string ansatz_data_file_name_str("no_ansatz_file"); // Default, since usually no ansatz will be given
//
///* Will need the settings file name without any preceding directory info,
//so we extract that first. Ditto for the initial data  file, and ansatz file.*/
//    std::string settings_file_name_str_final_piece = extract_Just_Filename(settings_file_name_str);
//    std::string initial_data_file_name_str_final_piece = extract_Just_Filename(initial_data_file_name_str);
//    std::string ansatz_data_file_name_str_final_piece("no_ansatz_file");
//
////Get ansatz file name if given
//    if (argc == 4) {
//        const char *ansatz_data_file_name = argv[3];
//        ansatz_data_file_name_str.assign(ansatz_data_file_name);
//        ansatz_data_file_name_str_final_piece = extract_Just_Filename(ansatz_data_file_name_str);
//        initWriteToLogStr += "Using ansatz data file " + ansatz_data_file_name_str_final_piece + "\n";
//    }
//
//
//// Do initial handling of directories and input files.
//    std::string outputDirName("this_string_will_hold_the_output_dir_name");
//    try {
//        outputDirName = directory_setup(settings_file_name_str,
//                                                  settings_file_name_str_final_piece,
//                                                  initial_data_file_name_str,
//                                                  initial_data_file_name_str_final_piece,
//                                                  ansatz_data_file_name_str,
//                                                  ansatz_data_file_name_str_final_piece,
//                                                  argc,
//                                                  initWriteToLogStr);
//    }
//    catch (const std::runtime_error &error) {
//        std::cerr << error.what() << std::endl;
//        std::cerr << "Initial directory and file handling failed." << std::endl;
//        return -1;
//    }
//
///* Set up log file that will store a copy of the std::cout and std::cerr
//output (with the exception of any writes to std::cerr that have already happened,
//and the one below if the log file creation failed). The logStream will be passed
//by reference to any functions that want to print output to the log file. The
//logStream must be opened before writing to it, and it is good practice to close
//it immediately afterwards.*/
//    CustomOutStreamClass logStream;
//    logStream.setOutputFileName(outputDirName + "/log.txt");
//// Check we can create a log file with the chosen name.
//    try {
//        logStream.open();
//    }
//    catch (const std::runtime_error &error) {
//        std::cerr << error.what() << std::endl;
//        return -1;
//    }
//
///* Set format and precision of logStream, and print some things related to what
//has already occurred above.*/
//    logStream << std::defaultfloat << std::setprecision(6);
//    logStream << initWriteToLogStr << std::endl;
//    logStream.close();
//
///*Read in settings file (using libconfig++ library) and put these settings in
//a Settings that will be passed around between functions when needed.
//See Settings.hpp for details of settings. */
//    Settings settings;
//    try { readSettingsFile(settings, settings_file_name); }
//    catch (const libconfig::FileIOException &fioex) {
//        logStream.open();
//        logStream << "I/O error reading settings file. Check file exists in correct directory and has correct name."
//                  << std::endl;
//        return -1;
//    }
//    catch (const libconfig::ParseException &pex) {
//        logStream << "Parse error reading settings file. Check file has correct format." << std::endl;
//        return -1;
//    }
//    catch (const std::runtime_error &error) {
//        logStream << error.what() << std::endl;
//        logStream.close();
//        return -1;
//    }
//// Calculate Young's Modulus
//    settings.youngs_modulus = 2.0 * settings.shear_modulus * (1.0 + settings.poisson_ratio);
//
//
//////////////////////////////////////////////////////////////////////////////////
//
///* Read in *initial* data for each node and triangle
//Each element of the node container 'nodes' is a class holding data for that node
//etc. seriesOf_InvProgMetrics is a vector of vectors, where each element of the
//outer vector holds a vector of (director angle, lambda, nu) vectors, one vector
//for each triangle. The different sets held in the outer vector are to allow a
//sequence of programmed tensors to be activated in turn, one after the other. For
//example a bend first in one direction, with a bend in a second direction only
//added afterwards. Ditto for programmed_second_fundamental_forms. Between each set of programmed
//tensors the full dialling in procedure occurs. Unless an ansatz file is used,
//the first round of dialling in will always start from the Identity and Zero for
//the programmed metric and secFF respectively. For most uses of this code, only
//a single set of programmed tensors will be supplied - multiple sets is a fairly
//advanced case.*/
//    std::vector<Node> nodes;
//    std::vector<Triangle> triangles;
//    std::vector<std::vector<Eigen::Vector3d>> programmed_metric_infos;
//    std::vector<std::vector<Eigen::Matrix<double, 2, 2> >> inverted_programmed_metrics;
//    std::vector<std::vector<double>> programmed_taus;
//    std::vector<std::vector<Eigen::Matrix<double, 2, 2> >> programmed_second_fundamental_forms;
//
//
///* This variable keeps track of where in the sequence of programmed tensors the
//dynamics should start from. The default (zero) is used unless an ansatz file was
//given, in which case it may be that evolution from the ansatz should start from
//further on in the sequence of programmed tensors. */
//    std::size_t initial_stage = 0;
//// This is the equivalent variable for dial_in_factor.
//    double dialInFactorToStartFrom = 0.0;
//// Create container to store node ansatz positions, if given.
//    std::vector<Eigen::Vector3d> nodeAnsatzPositions;
//
//    logStream.open();
//    logStream << "Now attempting to read data files. An error here likely \n"
//                 "implies a problem with a data file, for example a mismatch between the \n"
//                 "number of nodes or triangles stated and the number actually present; \n"
//                 "or other similar mismatches in numbers of data values; or a format problem. \n"
//                 "Remember the input files must have exactly the correct format." << std::endl;
//    logStream.close();
//    try {
//        readVTKData(nodes, triangles, programmed_metric_infos, inverted_programmed_metrics, programmed_taus,
//                    programmed_second_fundamental_forms, settings, initial_data_file_name_str,
//                    initial_stage,
//                    dialInFactorToStartFrom, nodeAnsatzPositions, ansatz_data_file_name_str, logStream);
//    }
//    catch (const std::out_of_range &out_of_bounds_error) {
//        logStream << out_of_bounds_error.what() << std::endl;
//        return -1;
//    }
//    catch (const std::runtime_error &reading_error) {
//        logStream << reading_error.what() << std::endl;
//        return -1;
//    }
//
//
//
//
//// Print numbers of nodes and triangles.
//    logStream.open();
//    logStream << "Number of nodes = " << settings.num_nodes << std::endl;
//    logStream << "Number of triangles = " << settings.num_triangles << std::endl;
///* Print warning if number of triangles is low - in this case boundary effects
//will dominate, and the code should not be trusted, both due to the reduced
//accuracy in the treatment of the boundary in e.g. the 2nd F.F. approx, and due
//to possible bugs in this rather special case.*/
//    if (settings.num_triangles < 50) {
//        logStream << "Your mesh has a small number of triangles. \nBeware that the code "
//                     "is likely to be less accurate in this case, \nand unforeseen bugs are more "
//                     "likely in extreme cases." << std::endl;
//    }
//    logStream.close();
//
///* Set the nodes and triangles debugging display functions to print to the
//same output log file as logStream (could send somewhere else in principle).*/
//    for (int i = 0; i < settings.num_nodes; ++i) {
//        nodes[i].nodeLogStream.setOutputFileName(logStream.getOutputFileName());
//    }
//    for (int i = 0; i < settings.num_triangles; ++i) {
//        triangles[i].triLogStream.setOutputFileName(logStream.getOutputFileName());
//    }
//
//
/////////////////////////////////////////////////////////////////////////////////
////Here we do some remaining geometry/topology/mesh-related things
//
//// Determine and store the labels of the triangles incident on each node.
//    calcTrianglesIncidentOnNodes(nodes, triangles);
//
///*Now calculate triangle edge-sharing adjacencies (triangle member data), and
//set up an 'edges' data structure to store further edge information.*/
//    std::vector<Edge> edges;
//    try { calcTriangleAdjacencies_And_Edges(nodes, triangles, edges); }
//    catch (const std::runtime_error &mesh_adjacencies_or_edges_error) {
//        logStream << mesh_adjacencies_or_edges_error.what() << std::endl;
//        return -1;
//    }
//
///* Set the edges debugging display functions to print to the same output log
//file as logStream (could send somewhere else in principle).*/
//    for (int i = 0; i < settings.num_edges; ++i) {
//        edges[i].edgeLogStream.setOutputFileName(logStream.getOutputFileName());
//    }
//
///* Calculate and print the number of boundary and non-boundary edges.
//Also label the nodes connected by boundary edges as boundary nodes.
//Also calculate the total initial perimeter of the sample, to be used as a
//characteristic sample length in estimating characteristic times, time steps
//etc.*/
//    int numBoundaryEdges = 0;
//    std::vector<double> initBoundaryEdgeLengths(settings.num_edges);
//
//    for (int i = 0; i < settings.num_edges; ++i) {
//        if (edges[i].isOnBoundary) {
//
//            numBoundaryEdges += 1;
//
//            nodes[edges[i].nodeLabels(0)].isOnBoundary = true;
//            nodes[edges[i].nodeLabels(1)].isOnBoundary = true;
//
//            //If chosen in settings, clamp whole boundary in addition to clamp
//            //indicators from data file.
//            if (settings.is_boundary_clamped) {
//                nodes[edges[i].nodeLabels(0)].isClamped = true;
//                nodes[edges[i].nodeLabels(1)].isClamped = true;
//            }
//
//            // Store initial length of this boundary edge
//            initBoundaryEdgeLengths[i] = (nodes[edges[i].nodeLabels(0)].pos - nodes[edges[i].nodeLabels(1)].pos).norm();
//        }
//    }
//
//// Do sum to calculate perimeter, and set the characteristic sample length to it.
//    double initPerimeter = kahanSum(initBoundaryEdgeLengths);
//    settings.sample_char_length = initPerimeter;
//    logStream.open();
//    logStream << "Initial perimeter = " << initPerimeter << std::endl;
//    logStream.close();
//
//
////A further check that things are ok:
//    try {
//        if (3 * settings.num_triangles != 2 * settings.num_edges - numBoundaryEdges) {
//            throw std::runtime_error(
//                    "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
//        }
//    }
//    catch (const std::runtime_error &mesh_adjacencies_or_edges_error) {
//        logStream << mesh_adjacencies_or_edges_error.what() << std::endl;
//        return -1;
//    }
//
//    logStream.open();
//    logStream << "Number of edges = " << settings.num_edges << std::endl;
//    logStream << "Number of boundary edges = " << numBoundaryEdges << std::endl;
//    logStream << "Number of non-boundary edges = " << settings.num_edges - numBoundaryEdges << std::endl;
//    logStream.close();
//
//
///* Now set triangle boundary indicators, based on whether they have any *edges*
//on the sample boundary.*/
//    int numBoundaryTriangles = 0;
//    for (int i = 0; i < settings.num_triangles; ++i) {
//        if (triangles[i].isOnBoundary) {
//            numBoundaryTriangles += 1;
//        }
//    }
//
//    logStream.open();
//    logStream << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
//    logStream << "Number of holes in mesh = " << 1 + settings.num_edges - settings.num_nodes - settings.num_triangles
//              << std::endl; // From Euler's formula for a planar graph.
//    logStream.close();
//    try {
//        if (1 + settings.num_edges - settings.num_nodes - settings.num_triangles < 0) {
//            throw std::runtime_error(
//                    "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
//        }
//    }
//    catch (const std::runtime_error &mesh_catastrophe_error) {
//        logStream << mesh_catastrophe_error.what() << std::endl;
//        return -1;
//    }
//
//
///* Determine and store labels of the node neighbours to each node (those that
//share an edge).*/
//    configureNodeAdjacency(nodes, edges);
//
//// Calculate and store number of boundary nodes.
//    int numBoundaryNodes = 0;
//    for (int n = 0; n < settings.num_nodes; ++n) {
//        if (nodes[n].isOnBoundary) {
//            numBoundaryNodes += 1;
//        }
//    }
//    settings.num_boundary_nodes = numBoundaryNodes;
//
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//
///* Now calculate a 3x6 matrix for each triangle that can be stored and used
//repeatedly to obtain the elements of the Second Fundamental Form for each
//triangle during the simulation. The 2nd F.F. is an important extrinsic geometric
//object, which completely specifies a surface (up to rigid motions) when combined
//with a metric tensor (aka 1st F.F.). This is called the Bonnet Theorem. An
//'intended' 2nd F.F. can be imposed with a nematic director that 'twists' through
//the thickness of an LCE sheet, and thus the intended and actual 2nd F.F.s enter
//into the energy that is being minimised in this simulation. Intuitively, the 2nd
//F.F. describes the way the deformed surface is curved in 3D space in terms of
//the directional derivative of the surface normal with respect to the 2D
//parametrisation coordinates on the surface. Its matrix representation is
//symmetric.
//A related object is found by 'raising an index' by left-multiplying by the
//inverse metric, which gives the 'Shape Operator'. This has determinant and trace
//equal to the Gaussian curvature and twice the mean curvature respectively.
//The approach taken to approximate the 2nd F.F. for a triangulated mesh in this
//code is based on 'patch fitting', which is simple, and has been shown to
//converge. Consider a triangle t. It has three vertices, and three other nearby
//vertices are selected, subject to the following quadratic fitting being
//well-conditioned.
//The unique quadratic surface (a function of the 2D parametrisation coords) that
//goes through all 6 of these points in the current state is then found and
//assigned to triangle t. The 2nd F.F. for this 'patch' could be calculated in
//full and averaged over the patch, although this complicates matters. Instead, we
//just take the normal to the triangle face and combine that with the 2nd
//derivatives of the patch surface following a particular expression for the secFF.
//The matrix pre-calculated here is later multiplied onto the current node
//positions of the 6 patch nodes to give patch surface derivatives.
//*/
//    logStream.open();
//    logStream << "\n" << "Beginning patch selection and related pre-calculations." << std::endl;
//    logStream.close();
//
//    try { calc_nonVertexPatchNodes_and_MatForPatchDerivs(nodes, triangles, logStream, 0); }
//    catch (const std::runtime_error &patch_error) {
//        logStream << patch_error.what() << std::endl;
//        return -1;
//    }
//    logStream.open();
//    logStream << "Successfully completed patch setup." << "\n" << std::endl;
//    logStream.close();
//
//
///*For the current and initial in-plane coordinate bases to match up, the node
//labels for each triangle must be ordered such that the resulting calculated
//normals for the initial flat state point along the +z axis, not -z. To do this
//I calculate the normals for an arbitrary label order, then loop over triangles,
//shuffling the labels for those with negative-z-component normals, before
//recalculating the initial normals, areas etc. Note, a more complex procedure
//will be required for non-planar initial geometries. The magic numbers in the
//argument list are just debugging values that should not matter here.*/
//    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, WaitingForEquilibrium, -12345, 98765,
//                                                  programmed_metric_infos, inverted_programmed_metrics,
//                                                  programmed_taus, programmed_second_fundamental_forms, settings);
//    Eigen::Vector3d tempZAxisVec;
//    tempZAxisVec << 0.0, 0.0, 1.0;
//    for (int i = 0; i < settings.num_triangles; ++i) {
//        if (tempZAxisVec.dot(triangles[i].currSides.col(0).cross(triangles[i].currSides.col(1))) < 0) {
//            int tempLabel = triangles[i].vertexLabels(2);
//            triangles[i].vertexLabels(2) = triangles[i].vertexLabels(1);
//            triangles[i].vertexLabels(1) = tempLabel;
//        }
//    }
//    std::cout << "CHECK THIS - should shuffle normals before patch matrix calc I think?" << std::endl;
//    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, WaitingForEquilibrium, -12345, 98765,
//                                                  programmed_metric_infos, inverted_programmed_metrics,
//                                                  programmed_taus, programmed_second_fundamental_forms, settings);
//
///* Set remaining initial conditions e.g. velocity, store initial triangle
//side components and fractional edge lengths, and calculate node masses.
//Also set the first set of programmed tensors to be the trivial ones for the
//plane, and set up the triangles to have their tensors dialled in from these to
//begin with. */
//    try {
//        setRemainingInitCond_and_NodeMasses(nodes, triangles, edges, programmed_metric_infos,
//                                            inverted_programmed_metrics, programmed_taus,
//                                            programmed_second_fundamental_forms,
//                                            settings);
//    }
//    catch (const std::runtime_error &setInitCondAndInitSidesEtc_error) {
//        logStream << setInitCondAndInitSidesEtc_error.what() << std::endl;
//        return -1;
//    }
//
////    settings.SetupSmallestElements(logStream, triangles, programmed_taus);
//
//    /* First find the approximate smallest linear size of mesh element, based on
//    smallest altitude of each triangle. Find also the smallest value of
//    this linear size / sqrt(progTau), which is what really matters for the stretching
//    time step. If progTau varied a a lot over the programmed_taus, you might
//    want to be less wasteful and recalculate the time step each time you move on to
//    the next set of progTaus in the sequence, because the time step then might be
//    much bigger for other sets in the sequence. If you wanted to get even fancier,
//    you could change the time step as prog_Tau is dialled in between two sets in the
//    sequence.*/
//    settings.approx_min_init_elem_size = DBL_MAX;
//    settings.smallest_size_over_root_tau = DBL_MAX;
////    int smallestTri = 0;
//    for (int i = 0; i < triangles.size(); ++i) {
//
//        // This variable will end up holding the smallest altitude for *this* tri.
//        double smallestAltitude =
//                2 * triangles[i].initArea / (triangles[i].currSides.col(0) - triangles[i].currSides.col(1)).norm();
//
//        for (int s = 0; s < 2; ++s) {
//            if (smallestAltitude > 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm()) {
//                smallestAltitude = 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm();
//            }
//        }
//
//        if (settings.approx_min_init_elem_size > smallestAltitude) {
//            settings.approx_min_init_elem_size = smallestAltitude;
//        }
//
//        for (auto &sequenceOf_ProgTau: programmed_taus) {
//            if (settings.smallest_size_over_root_tau > smallestAltitude / sqrt(sequenceOf_ProgTau[i])) {
//                settings.smallest_size_over_root_tau = smallestAltitude / sqrt(sequenceOf_ProgTau[i]);
//            }
//        }
//    }
//    logStream.open();
//    logStream << "Sheet thickness = " << settings.sheet_thickness << std::endl;
//    logStream << "Approx smallest element linear size = " << settings.approx_min_init_elem_size << std::endl;
//    logStream.close();
//    settings.SetupDialIn(logStream);
//
///* Set settings.TimeBetweenEquilChecks based on the tunable dimensionless value
//in the settings file, which relates this time to settings.DialInStepTime.*/
//    settings.time_between_equil_checks = settings.time_between_equil_checks_prefactor * settings.dial_in_step_time;
//    settings.SetupStepTime(logStream);
//    settings.SetupPrintFrequency(logStream);
//
//
//// Print total load force.
//    int numLoadedNodes = 0;
//    for (int i = 0; i < settings.num_nodes; ++i) {
//        if (nodes[i].isLoadForceEnabled) {
//            numLoadedNodes += 1;
//        }
//    }
//    logStream.open();
//    logStream << "Total load force applied = " <<
//                                               numLoadedNodes * settings.load_strength * settings.shear_modulus * settings.approx_min_init_elem_size *
//                                               settings.sheet_thickness << std::endl;
//    logStream.close();
//
//
//// Calculate and store characteristic force, energy, and energy density scales.
//    settings.char_force_scale = settings.shear_modulus * settings.approx_min_init_elem_size * settings.sheet_thickness;
//    settings.char_stretch_energy_density_scale = settings.shear_modulus * settings.sheet_thickness;
//// settings.charBendEnergyDensityScale = settings.ShearModulus * settings.SheetThickness * settings.SheetThickness * settings.SheetThickness / (settings.SampleCharLength * settings.SampleCharLength);
//    double totInitArea = 0;
//    std::vector<double> initAreas(settings.num_triangles);
//    for (int i = 0; i < settings.num_triangles; ++i) {
//        initAreas[i] = triangles[i].initArea;
//    }
//    totInitArea = kahanSum(initAreas);
//    settings.char_stretch_energy_scale = settings.char_stretch_energy_density_scale * totInitArea;
//// settings.charBendEnergyScale = settings.charBendEnergyDensityScale * totInitArea;
//
//
//
///* Variables used to control the programmed properties' time variation.
//This is not necessarily directly physical, especially for a fast optical
//activation instead of slow heating; the idea here is to gradually `dial in'
//the programmed (energetically favoured) metric and second fundamental form.
//This prevents anything too explosive happening the code, and can be used to
//make the simulation quasistatic (which can help avoid undesired isometries
//for example), though it may not be necessary.
//The range between 0, and 1 is split into a sequence of discrete values based
//on the chosen settings.DialInResolution. The dial-in factor will then
//evolve linearly from one such value to the next, over a chosen
//settings.DialInStepTime. After each such step, the value is held constant
//until equilibrium is reached (to within the chosen tolerance). In the
//'waiting' stage, the check for equilibrium occurs at a rate determined by
//the settings.TimeBetweenEquilChecks. The dialling in process occurs
//between each pair of the sequence of programmed tensors, if more than one is
//given in the input file.
//*/
//    auto timeSinceLastEquilCheck = DBL_MAX; // Reset to zero every time equil is checked, or every time a new DialInFactor value is reached.
//    auto timeSinceCurrDiallingInPhaseStarted = DBL_MAX; // Reset to zero every time equil is reached and a new `dialling in' phase starts.
//    SimulationStatus status; // Indicates whether the simulation is currently in a 'dialling in' phase, or is not dialling and is instead waiting for equilibrium, or whether equilibrium has been reached but the next dialling phase has not yet begun.
//    std::size_t phase_counter = SIZE_MAX; // Keeps track of which DialInFactor value in dial_in_phases was last dialled *from* (not held at, which is the next value along in the list)
//    auto dial_in_factor = DBL_MAX; // Keeps track of current value
//
//    std::vector<double> dial_in_phases;
//    double tempDialInFactor = 0.0;
//    if (settings.dial_in_resolution > 0 && settings.dial_in_step_time >= 0) {
//        while (tempDialInFactor < 1.0) {
//            dial_in_phases.push_back(tempDialInFactor);
//            tempDialInFactor += settings.dial_in_resolution;
//        }
//        dial_in_phases.push_back(1.0);
//    } else {
//        logStream.open();
//        logStream << "settings.DialInResolution and settings.DialInStepTime must "
//                     "both be >0 and >=0 respectively, and they aren't currently. Aborting."
//                  << std::endl;
//        logStream.close();
//        return -1;
//    }
//
///* Create std::vectors (with one element for each triangle) that will be passed
//by reference to functions calculating the Gauss and mean curvatures, and
//stretching and bending energies and energy densities. */
//    std::vector<double> gaussCurvatures(settings.num_triangles,
//                                        DBL_MAX); //Recognisable initialisation for debugging.
//    std::vector<double> meanCurvatures(settings.num_triangles, DBL_MAX);
//    std::vector<double> stretchEnergies(settings.num_triangles, DBL_MAX);
//    std::vector<double> bendEnergies(settings.num_triangles, DBL_MAX);
//    std::vector<double> stretchEnergyDensities(settings.num_triangles, DBL_MAX);
//    std::vector<double> bendEnergyDensities(settings.num_triangles, DBL_MAX);
//    std::vector<double> kineticEnergies(settings.num_nodes, DBL_MAX);
//    std::vector<double> strainMeasures(settings.num_triangles, DBL_MAX);
//    std::vector<Eigen::Vector2d> cauchyStressEigenvals(settings.num_triangles);
//    std::vector<Eigen::Matrix<double, 3, 2>> cauchyStressEigenvecs(settings.num_triangles);
//
//
///* Create further std::vectors to store angle deficits if specified in settings
//file.*/
//    std::vector<double> angleDeficits(settings.num_nodes, DBL_MAX);
////std::vector<double> interiorNodeAngleDeficits(settings.NumNodes - settings.numBoundaryNodes, DBL_MAX);
////std::vector<double> boundaryNodeAngleDeficits(settings.numBoundaryNodes, DBL_MAX);
//    std::cout << "CHECK THIS" << std::endl;
//    std::vector<double> interiorNodeAngleDeficits(settings.num_nodes, DBL_MAX);
//    std::vector<double> boundaryNodeAngleDeficits(settings.num_nodes, DBL_MAX);
//
///*
//std::cout<< "STORING NODE REF POSITIONS TO HELP MAKE EXACT CONE." << std::endl;
//std::vector< Eigen::Vector3d > nodeRefPosits(settings.NumNodes);
//for(int n = 0; n < settings.NumNodes; ++n){
//    nodeRefPosits[n] = nodes[n].pos;
//}*/
//
///* Perturb node positions with small random noise if desired, to allow 'breaking
//away' from flat plane initial condition, for example. This will have no effect
//if an ansatz is used. */
//    if (settings.is_perturbation_of_initial_positions_enabled) {
//        perturbInitialPositionsWithRandomNoise(nodes, 0);
//    }
//
//////////////////////////////////////////////////////////////////////////////////
//
///* Initialise time and step counter. Set logStream to more useful
//format and precision for e.g. energy printouts. */
//    double time = 0;
//    int stepCount = 0;
//    logStream << std::scientific << std::setprecision(8);
//
//    std::ofstream forceDistFile;
//
//
///* Loop over the sequence of programmed tensors, dialling-in and waiting for
//equilibrium between each pair in the sequence. This loop is redundant in most use
//cases for this code, where only a single set of programmed tensors is supplied.*/
//    for (std::size_t progTensorSequenceCounter = initial_stage;
//         progTensorSequenceCounter <= programmed_metric_infos.size() - 2; ++progTensorSequenceCounter) {
//
//        /* A third input file specifying an ansatz for the node positions may have
//        been read in if it was given as a command line argument. This could for
//        example correspond to an output file from this code, to carry on where some
//        previous simulation left off. In this case, the nodes are moved to their
//        ansatz positions here, and other relevant variables are set up. */
//        if (ansatz_data_file_name_str != "no_ansatz_file" &&
//            progTensorSequenceCounter == initial_stage) {
//
//            for (int i = 0; i < settings.num_nodes; ++i) {
//                nodes[i].pos = nodeAnsatzPositions[i];
//            }
//
//            timeSinceLastEquilCheck = 0.0;
//            status = Dialling;
//
//            if (!settings.is_dialing_from_ansatz_enabled) {
//
//                dial_in_factor = dialInFactorToStartFrom;
//                phase_counter = 0;
//                for (std::size_t i = 0; i < dial_in_phases.size(); ++i) {
//                    if (dial_in_phases[i] < dial_in_factor) {
//                        phase_counter = i;
//                    }
//                }
//
//                timeSinceCurrDiallingInPhaseStarted =
//                        ((dial_in_factor - dial_in_phases[phase_counter]) /
//                         (dial_in_phases[phase_counter + 1]
//                          - dial_in_phases[phase_counter])) * settings.dial_in_step_time;
//
//                /* Here we implement another feature, in which the progMetric and
//                progSecFF that we would normally be dialling FROM at this point, are
//                just set to the values they would be dialling TO, so that no dialling
//                of these quantities happens between the current (starting) entry in
//                the progTensor sequence and the next. progTau however IS dialled in as
//                normal. Upon reaching the next entry in the sequence, normal behaviour
//                resumes so there will be a phase of waiting for equilibrium,
//                followed potentially by normal dialling in of further progTensors in
//                the sequence, if there are any more in the sequence.
//                This feature's intended use case is where you want to dial in some
//                large progTaus on the boundary triangles, to reduce stretch-bend
//                tradeoff boundary effects, but you don't want to do any dialling of
//                the progMetric or progSecFF because you think the shape is already
//                close to these in its achieved metric and secFF, probably because
//                you already did a simulation with them (but no large progTau), and
//                are using the output of that run as an ansatz. This does mean that
//                you'll want to CHANGE THE dial_in_factor IN THE ANSATZ FILE TO = 0.
//                Otherwise the whole dialling in phase will be skipped entirely, and
//                the progTaus will just jump to their new large values, which may well
//                produce violent behaviour. You might be tempted to just use the
//                isDialingFromAnsatzEnabled feature instead of this one, since
//                then the new progTaus are jumped to, but the progMetric and secFF are
//                dialled in from the initial actual metric and actual secFF values,
//                ensuring gentle dynamics. This is probably non-ideal though, because
//                the progMetric and progSecFF along that dialling path may then
//                deviate significantly from the form that you want to arrive at, and
//                that you thought you started close to, pumping unnecessary energy in
//                to the system and wasting time.
//                There may be other uses for this feature too.*/
//                if (settings.for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff) {
//                    for (int i = 0; i < settings.num_triangles; ++i) {
//                        programmed_second_fundamental_forms[initial_stage][i] = programmed_second_fundamental_forms[
//                                initial_stage + 1][i];
//                        programmed_metric_infos[initial_stage][i] = programmed_metric_infos[
//                                initial_stage + 1][i];
//                    }
//                }
//            }
//
//                /* The isDialingFromAnsatzEnabled setting is very handy if you can
//                only make an approximate guess at a solution, and want to start from
//                there and see where things go without explosive dynamics.
//                If settings.isDialingFromAnsatzEnabled == true, we instead
//                arrange to ignore the dial_in_factor, phase_counter and
//                time_phase variables that are implied by the
//                the ansatz data file, and instead just start a new dialling in phase
//                from whatever the *calculated* metric and secFF are in the ansatz state.
//                Tau is not a calculable geometric quantity so we just set it to the
//                programmed value immediately, which will not cause explosive dynamics
//                even if e.g. the ansatz shape came from a simulation with a very different
//                programmed tau, because the stretch and bend energies both still start
//                at zero etc.
//                */
//            else {
//                logStream.open();
//                logStream
//                        << "\nsettings.isDialingFromAnsatzEnabled == true, FIRST DIALLING IN PHASE WILL BE FROM ANSATZ STATE."
//                        << std::endl;
//                logStream.close();
//
//                dial_in_factor = 0.0;
//                phase_counter = 0;
//                timeSinceCurrDiallingInPhaseStarted = 0.0;
//
//                // Calculate all necessary geometry for the ansatz state.
//                calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, WaitingForEquilibrium, dial_in_factor,
//                                                              progTensorSequenceCounter, programmed_metric_infos,
//                                                              inverted_programmed_metrics, programmed_taus,
//                                                              programmed_second_fundamental_forms, settings);
//                updateSecondFundamentalForms(triangles, settings);
//
//
//                // Temp LU decomp of triangle metric, used to check invertibility.
//                Eigen::FullPivLU<Eigen::Matrix<double, 2, 2> > tempMetricDecomp;
//
//
//                /* Alter inverted_programmed_metrics[initial_stage] and similar to change where
//                the programmed quantities are dialling from.*/
//
//                for (int i = 0; i < settings.num_triangles; ++i) {
//
//                    // Test for invertibility of metric before taking inverse.
//                    try {
//                        tempMetricDecomp.compute(
//                                (triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse());
//                        if (!tempMetricDecomp.isInvertible()) {
//                            throw std::runtime_error(
//                                    "At least one triangle had a non-invertible metric in the ansatz state. \n"
//                                    "This should not occur in a reasonable mesh. Aborting.");
//                        } else {
//                            inverted_programmed_metrics[initial_stage][i] = (
//                                    triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
//                            if (settings.is_dialing_disabled) {
//                                inverted_programmed_metrics[initial_stage + 1][i] = (
//                                        triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
//                            }
//                        }
//                    }
//                    catch (const std::runtime_error &nonInvertibleMetric_error) {
//                        logStream.open();
//                        logStream << nonInvertibleMetric_error.what() << std::endl;
//                        logStream.close();
//                        return -1;
//                    }
//
//                    // Other quantities require no inverse, so are easier.
//                    programmed_taus[initial_stage][i] = programmed_taus[
//                            initial_stage + 1][i];
//                    programmed_second_fundamental_forms[initial_stage][i] = triangles[i].secFF;
//                    if (settings.is_dialing_disabled) {
//                        programmed_taus[initial_stage + 1][i] = programmed_taus[
//                                initial_stage + 1][i];
//                        programmed_second_fundamental_forms[initial_stage +
//                                                            1][i] = triangles[i].secFF;
//                    }
//                }
//            }
//
//            /////////////////////////////////////////////////////////////////////
//
//            // If set in settings file, clamp all nodes within a certain vertical distance
//            // above the node with the lowest initial z value (only works with ansazt clearly).
//            if (settings.thicknesses_above_lowest_node_to_clamp_up_to > 0) {
//                double minNodeZCoord = nodes[0].pos(2);
//                for (int n = 0; n < settings.num_nodes; ++n) {
//                    if (minNodeZCoord > nodes[n].pos(2)) {
//                        minNodeZCoord = nodes[n].pos(2);
//                    }
//                }
//
//                for (int n = 0; n < settings.num_nodes; ++n) {
//                    if (nodes[n].pos(2) - minNodeZCoord <
//                        settings.thicknesses_above_lowest_node_to_clamp_up_to * settings.sheet_thickness) {
//                        nodes[n].isClamped = true;
//                    }
//                }
//            }
//
//            ///////////////////////////////////////////////////////////////////
//        }
//
//            // If no ansatz is being used:
//        else {
//            // Reset the variables that control each dialling/waiting process.
//            timeSinceLastEquilCheck = 0.0;
//            status = Dialling;
//            dial_in_factor = 0.0;
//            phase_counter = 0;
//            timeSinceCurrDiallingInPhaseStarted = 0.0;
//        }
//
//
//        /* Handle cases where settings.DialInStepTime is zero or very small. This is
//        a slightly hacky way to ensure that no dialling actually occurs, and the
//        simulation jumps to a status = waitingForEquilibrium state. The magic number is
//        just chosen to be recognisable for debugging.*/
//        if (settings.dial_in_step_time < settings.time_step) {
//            settings.dial_in_step_time = 0.0;
//            timeSinceCurrDiallingInPhaseStarted = 1.23456789;
//        }
//
//        ///////////////////////////////////////////////////////////////////////////
//
//        // Begin dynamical evolution of node positions and velocities.
//        logStream.open();
//        logStream << "\nBeginning dynamical evolution.\n" << std::endl;
//        logStream.close();
//
//        logStream.open();
//        logStream << "\nCREATING VECTOR TO STORE UNSTRESSED CONE NODE POSITIONS.\n" << std::endl;
//        logStream.close();
//        std::vector<Eigen::Vector3d> nodeUnstressedConePosits(settings.num_nodes);
//        double s1 = 12345678.9;
//        //double sTest = 98765432.1;
//
//        std::vector<std::vector<std::pair<
//                int, int>>> correspondingTrianglesForNodes = getCorrespondingTrianglesForNodes(
//                triangles, nodes);
//
//        while (phase_counter <= dial_in_phases.size() - 2) {
//            int highestNode = -99;
//            int lowestNode = -99;
//            if (stepCount == 0) {
//                settings.init_slide_z_coord_lower = nodes[0].pos(2);
//                settings.init_slide_z_coord_upper = nodes[0].pos(2);
//                for (int n = 0; n < settings.num_nodes; ++n) {
//                    if (settings.init_slide_z_coord_lower > nodes[n].pos(2)) {
//                        settings.init_slide_z_coord_lower = nodes[n].pos(2);
//                        lowestNode = n;
//                    }
//                    if (settings.init_slide_z_coord_upper < nodes[n].pos(2)) {
//                        settings.init_slide_z_coord_upper = nodes[n].pos(2);
//                        highestNode = n;
//                    }
//                    // Can instead choose initial slide separation directly from settings file, to help avoid 'jumping' when starting a
//                    // squashing run part way through from a previous run's output.
//                    if (settings.specify_init_slide_z_coord_upper > -99.0) {
//                        settings.init_slide_z_coord_upper = settings.specify_init_slide_z_coord_upper;
//                    }
//                    if (settings.specify_init_slide_z_coord_lower > -99.0) {
//                        settings.init_slide_z_coord_lower = settings.specify_init_slide_z_coord_lower;
//                    }
//
//                    // To do constant-weight slide instead
//                    if (settings.const_slide_weight_fac > 0) {
//                        settings.init_slide_z_coord_upper = settings.init_slide_z_coord_lower + settings.spacer_height *
//                                                                                                settings.sample_char_length;// This used to be settings.SpacerHeight * settings.SampleCharLength instead, which is wrong
//                    }
//
//                    settings.upper_slide_displacement = 0.0;
//                    settings.upper_slide_vel = 0.0;
//                    settings.is_slide_just_equilibrated = 0;
//                    if (settings.is_controlled_force_enabled) {
//                        settings.upper_slide_weight = settings.initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq *
//                                                      (settings.shear_modulus * settings.sheet_thickness *
//                                                       settings.sheet_thickness);
//                    }
//                }
//                // Uncomment the following line for single ridge experiment with point(ish) load to applied tip.
//                //settings.initSlideZCoord_upper = settings.initSlideZCoord_lower + 2.0 * settings.SheetThickness; // to apply point(ish) load to tip
//
//                // Readjust for the case of glass cones instead of glass slides. The
//                // slide Z coords correspond to the tips of the glass cones.
//                if (settings.glass_cones) {
//                    settings.cone_angle = 1.02327019;
//                    settings.init_slide_z_coord_upper += -tan(settings.cone_angle) *
//                                                         sqrt(nodes[highestNode].pos(0) * nodes[highestNode].pos(0) +
//                                                           nodes[highestNode].pos(1) * nodes[highestNode].pos(1));
//                    settings.init_slide_z_coord_lower += -tan(settings.cone_angle) *
//                                                         sqrt(nodes[lowestNode].pos(0) * nodes[lowestNode].pos(0) +
//                                                           nodes[lowestNode].pos(1) * nodes[lowestNode].pos(1));
//                    logStream.open();
//                    logStream << "USING TWO GLASS CONES FOR SQUASHING." << std::endl;
//                    logStream.close();
//
//                    // Hijack nodes[i].isOnBoundary to instead label nodes whose
//                    // forces we will modify to kill any components not tangential to
//                    // a perfect cone base state. We do this to nodes within an intermediate
//                    // distance vertically from the ends of the cone.
//                    double intermLengthScaleUpper = 2.0 * sqrt(settings.sheet_thickness * 0.18);
//                    double intermLengthScaleLower = 2.0 * sqrt(settings.sheet_thickness * 1.8);
//                    logStream.open();
//                    logStream
//                            << "Apply normal-force-killer within the following vertical distances of the top and bottom: "
//                            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
//                    logStream.close();
//                    for (int n = 0; n < settings.num_nodes; ++n) {
//                        if (((nodes[highestNode].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
//                            ((nodes[n].pos(2) - nodes[lowestNode].pos(2)) < intermLengthScaleLower)) {
//                            nodes[n].isOnBoundary = true;
//                        }
//                    }
//                }
//
//                // Instead set up imposed-Seide-deformations idea.
//                if (settings.is_seide_deformations_enabled) {
//                    settings.lambda = 0.9;
//                    settings.cone_angle = asin(pow(settings.lambda, 1.5));
//                    logStream.open();
//                    logStream << "IMPOSING SEIDE DEFORMATIONS." << std::endl;
//                    logStream.close();
//
//                    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
//                    // positions we will force to be those of Seide's setup (found
//                    // by findng the membrane stress solution for his base state, and
//                    // integrating the strains etc).
//                    double intermLengthScaleUpper = 3.0 * sqrt(settings.sheet_thickness * 0.18);
//                    double intermLengthScaleLower = 3.0 * sqrt(settings.sheet_thickness * 1.8);
//                    logStream.open();
//                    logStream
//                            << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
//                            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
//                    logStream.close();
//                    for (int n = 0; n < settings.num_nodes; ++n) {
//                        if (((nodes[highestNode].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
//                            ((nodes[n].pos(2) - nodes[lowestNode].pos(2)) < intermLengthScaleLower)) {
//                            nodes[n].isSeideDisplacementEnabled = true;
//                        }
//                    }
//
//                    // STORE PERFECT CONE ANSATZ
//                    for (int n = 0; n < settings.num_nodes; ++n) {
//                        nodeUnstressedConePosits[n] = nodes[n].pos;
//                    }
//                    s1 = sqrt(nodes[highestNode].pos(0) * nodes[highestNode].pos(0) +
//                              nodes[highestNode].pos(1) * nodes[highestNode].pos(1)) / sin(settings.cone_angle);
//
//                    Eigen::Vector3d testTriCurrCentroid;
//                    testTriCurrCentroid = (nodes[triangles[settings.test_triangle].vertexLabels(0)].pos +
//                                           nodes[triangles[settings.test_triangle].vertexLabels(1)].pos +
//                                           nodes[triangles[settings.test_triangle].vertexLabels(2)].pos) / 3;
//                    //sTest = sqrt(testTriCurrCentroid(0)*testTriCurrCentroid(0) + testTriCurrCentroid(1)*testTriCurrCentroid(1)) / sin(settings.ConeAngle);
//                }
//            }
//
//
//            if (!settings.is_controlled_force_enabled) {
//                settings.upper_slide_displacement =
//                        time * settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time;
//            } else { // settings.isControlledForceEnabled == true instead
//                //settings.upperSlideWeight = (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) * (time * settings.slideSpeedPrefactor / bending_long_time);
//                settings.slide_damping_param =
//                        0.4 * settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness /
//                        (settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time);
//                if (fabs(settings.upper_tot_slide_force + settings.upper_slide_weight) /
//                    (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) <
//                    settings.total_slide_force_to_mu_t_sq_ratio_equil_threshold
//                    && settings.const_slide_weight_fac < 0) {
//                    if (timeSinceLastEquilCheck > settings.time_between_equil_checks) {
//                        if (equilibriumCheck(nodes, settings, triangles, logStream) == EquilibriumReached) {
//                            settings.is_slide_just_equilibrated = 1;
//                            settings.upper_slide_weight +=
//                                    settings.slide_weight_dial_speed_fac *
//                                    (settings.time_step / settings.bending_long_time) *
//                                    (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness);
//                        }
//                        timeSinceLastEquilCheck = 0.0;
//                    }
//                }
//
//                // To do constant-weight slide instead
//                if (settings.const_slide_weight_fac > 0) {
//                    settings.upper_slide_weight =
//                            settings.const_slide_weight_fac * settings.shear_modulus * settings.sheet_thickness *
//                            settings.sheet_thickness;
//                }
//            }
//            settings.curr_slide_z_coord_upper = settings.init_slide_z_coord_upper - settings.upper_slide_displacement;
//            std::pair<double, double> upperAndLowerTotSlideForces;
//
//            // Impose Seide deformations.
//            if (settings.is_seide_deformations_enabled) {
//                double pInit = 3.0 * (settings.shear_modulus * settings.sheet_thickness *
//                                      settings.sheet_thickness); // So we don't have to start all the way from p=0, chosen based on previous sims.
//                settings.p = pInit + time * settings.p_speed_prefactor * settings.shear_modulus * settings.sheet_thickness *
//                                     settings.sheet_thickness / settings.bending_long_time;
//
//                for (int n = 0; n < settings.num_nodes; ++n) {
//                    if (nodes[n].isSeideDisplacementEnabled || stepCount == 0) {
//
//                        double tCone = settings.sheet_thickness / sqrt(settings.lambda);
//
//                        double polarAng = atan2(nodeUnstressedConePosits[n](1), nodeUnstressedConePosits[n](0));
//                        double s = sqrt(nodeUnstressedConePosits[n](0) * nodeUnstressedConePosits[n](0) +
//                                        nodeUnstressedConePosits[n](1) * nodeUnstressedConePosits[n](1)) /
//                                   sin(settings.cone_angle);
//                        Eigen::Vector3d uHat;
//                        uHat << sin(settings.cone_angle) * cos(polarAng), sin(settings.cone_angle) * sin(polarAng), -cos(
//                                settings.cone_angle);
//                        Eigen::Vector3d wHat;
//                        wHat << -cos(settings.cone_angle) * cos(polarAng), -cos(settings.cone_angle) *
//                                                                           sin(polarAng), -sin(settings.cone_angle);
//                        double u = -settings.p * log(s / s1) /
//                                   (2.0 * M_PI * settings.youngs_modulus * tCone * sin(settings.cone_angle) *
//                                    cos(settings.cone_angle));
//                        double w = -settings.p * (log(s / s1) + settings.poisson_ratio) /
//                                   (2.0 * M_PI * settings.youngs_modulus * tCone * cos(settings.cone_angle) *
//                                    cos(settings.cone_angle));
//                        nodes[n].pos = nodeUnstressedConePosits[n] + u * uHat + w * wHat;
//                    }
//                }
//            }
//
//
//            /* The forces are set each timestep using a += procedure, and therefore
//            must be set to zero each time before this is done.*/
//            zeroForces(nodes);
//
//            // Check if still in dialling in phase or whether it is time to wait for
//            //equilibrium.
//            if (timeSinceCurrDiallingInPhaseStarted >= settings.dial_in_step_time && status == Dialling) {
//
//                /* If the current dialling phase has indeed finished, set
//                dial_in_factor to the value that was being dialled up to in
//                that phase, and set the programmed tensors accordingly.*/
//                dial_in_factor = dial_in_phases[phase_counter + 1];
//                calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, status, dial_in_factor,
//                                                              progTensorSequenceCounter, programmed_metric_infos,
//                                                              inverted_programmed_metrics, programmed_taus,
//                                                              programmed_second_fundamental_forms, settings);
//                status = WaitingForEquilibrium;
//                // NB an EquilCheck has not actually just occurred, but this has the
//                //desired effect of ensuring that each DialInFactor value is held
//                //for at least one TimeBetweenEquilChecks.
//                timeSinceLastEquilCheck = 0.0;
//
//                settings.num_damp_factor =
//                        settings.equilibriation_damping *
//                        settings.damping_scale; // Set damping factor to waiting phase value.
//
//                logStream.open();
//                logStream << "Reached Dial-In Factor of " << dial_in_phases[phase_counter + 1]
//                          << ", now waiting for equilibrium" << std::endl;
//                logStream.close();
//            }
//
//            /* If not waiting for equilibrium, set the current value of the Dial-In
//            Factor, based on linear dialling in between the previously calculated
//            'checkpoint' values. */
//            if (status == Dialling) {
//
//                dial_in_factor = dial_in_phases[phase_counter] +
//                                   (dial_in_phases[phase_counter + 1]
//                                    - dial_in_phases[phase_counter]) *
//                                   timeSinceCurrDiallingInPhaseStarted / settings.dial_in_step_time;
//            }
//
//            /* Calculate sides, areas...etc of each triangle, as well as the current
//            dialled-in programmed (inverse) metric and second fundamental form.*/
//
//            std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
//            calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, status, dial_in_factor,
//                                                          progTensorSequenceCounter, programmed_metric_infos,
//                                                          inverted_programmed_metrics, programmed_taus,
//                                                          programmed_second_fundamental_forms, settings);
//
//            std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
//
//            /* Calculate secFF estimates for triangles, and related quantities such
//            as the derivative of the bending energy wrt the secFF components.*/
//            updateSecondFundamentalForms(triangles, settings);
//            std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
//
//            // Calculate current strain and bending force on each node.
//            update_elastic_forces(nodes, triangles, settings, correspondingTrianglesForNodes);
//            std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
//
//            /* Add force contributions from e.g. damping, loads, 'prod' perturbation, and
//             account for BCs e.g. clamping. */
//            upperAndLowerTotSlideForces = calcNonDeformationForces_and_ImposeBCS(nodes, time, settings);
//            std::chrono::steady_clock::time_point end4 = std::chrono::steady_clock::now();
//            settings.upper_tot_slide_force = upperAndLowerTotSlideForces.first;
//
//
//            // Hijack upperSlideDisplacement to hold end displacement in pulling experiment.
//            /*
//            if( settings.LoadStrength > 0.0 ){
//                upperSlideDisplacement = nodes[0].pos(0);
//                for(int n = 0; n < settings.NumNodes; ++n){
//                    if(upperSlideDisplacement < nodes[n].pos(0)){
//                        upperSlideDisplacement = nodes[n].pos(0);
//                    }
//                }
//            }*/
//
//            /* Write output data regularly.
//            Can be switched off with settings.PrintFrequency < 0.0.
//            Doing the write-out at this point in the loop means the node positions
//            and the triangle geometry data match in the output, which is desirable!*/
//            if ((stepCount % settings.inverse_print_rate == 0 && settings.inverse_print_rate > 0) ||
//                settings.is_slide_just_equilibrated == 1) {
//                try {
//                    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
//                                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
//                    if (settings.is_energy_densities_printed) {
//                        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
//                                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
//                    }
//
//                    // Here we hack the angleDeficits to instead tell us whether a node has a Seide displacement imposed or not.
//                    for (int n = 0; n < settings.num_nodes; ++n) {
//                        if (nodes[n].isSeideDisplacementEnabled) {
//                            angleDeficits[n] = 1;
//                        } else {
//                            angleDeficits[n] = 0;
//                        }
//                    }
//                    writeVTKDataOutput(nodes, triangles, stepCount, time, dial_in_factor, progTensorSequenceCounter,
//                                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
//                                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
//                                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
//                }
//                catch (const std::runtime_error &writing_error) {
//                    logStream << writing_error.what() << std::endl;
//                    return -1;
//                }
//                logStream.open();
//                logStream << std::fixed << "Wrote VTK output at " << getRealTime() << ", stepCount = " << stepCount
//                          << ", simulation time = " << time + settings.time_step << ", current dial-in factor = "
//                          << dial_in_factor << std::scientific;
//                logStream << ", last step's execution time "
//                          << std::chrono::duration_cast<std::chrono::microseconds>(end4 - begin).count() << " us"
//                          << std::endl;
//
////                logStream << "\tcalcTriangleGeometries execution time "
////                          << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin).count() << " us"
////                          << std::endl;
////                logStream << "\tcalcSecFFsAndRelatedQuantities execution time "
////                          << std::chrono::duration_cast<std::chrono::microseconds>(end2 - end1).count() << " us"
////                          << std::endl;
////                logStream << "\tcalcDeformationForces execution time "
////                          << std::chrono::duration_cast<std::chrono::microseconds>(end3 - end2).count() << " us"
////                          << std::endl;
////                logStream << "\tcalcNonDeformationForces_and_ImposeBCS execution time "
////                          << std::chrono::duration_cast<std::chrono::microseconds>(end4 - end3).count() << " us"
////                          << std::endl;
//
//                logStream.close();
//
//                forceDistFile.open(outputDirName + "/force_displacement_vals.txt", std::ofstream::app);
//                forceDistFile << std::scientific << std::setprecision(14)
//                              << upperAndLowerTotSlideForces.first /
//                                 (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) << " "
//                              << upperAndLowerTotSlideForces.second /
//                                 (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) << " "
//                              << settings.upper_slide_displacement / settings.sample_char_length << " "
//                              << settings.curr_slide_z_coord_upper << " "
//                              << settings.init_slide_z_coord_lower << " "
//                              << settings.is_slide_just_equilibrated << " "
//                              << settings.upper_slide_weight /
//                                 (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness)
//
//                              //<< settings.p / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " "
//                              //<< cauchyStressEigenvals.at(settings.testTriangle)(0) * sTest / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " "
//                              //<< cauchyStressEigenvals.at(settings.testTriangle)(1) * sTest / (settings.ShearModulus*settings.SheetThickness*settings.SheetThickness) << " " << std::endl;
//                              << std::endl;
//                forceDistFile.close();
//                settings.is_slide_just_equilibrated = 0;
//            }
//
//
//            if (!settings.is_controlled_force_enabled) {
//
//
//                /* If last check for equilibrium or last reaching of a new DialInFactor
//                 value was more than TimeBetweenEquilChecks ago, check for equilibrium.
//                 Also print total stretching and bending energies.*/
//                if (timeSinceLastEquilCheck > settings.time_between_equil_checks && status == WaitingForEquilibrium) {
//                    logStream.open();
//                    logStream << "Checking for equilibrium at " << getRealTime() << ", stepCount = " << stepCount
//                              << ", simulation time = " << time << ", current dial-in factor = " << dial_in_factor
//                              << std::endl;
//                    logStream.close();
//
//                    status = equilibriumCheck(nodes, settings, triangles, logStream);
//                    timeSinceLastEquilCheck = 0.0;
//
//                    calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
//                                            stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                            cauchyStressEigenvals, cauchyStressEigenvecs, settings);
//                    double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.char_stretch_energy_scale;
//                    double nonDimBendEnergy = kahanSum(bendEnergies) / settings.char_stretch_energy_scale;
//                    double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.char_stretch_energy_scale;
//
//                    logStream.open();
//                    logStream << "Non-dimensionalised energies:" << std::endl;
//                    logStream << "\t total \t\t" << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy
//                              << std::endl;
//                    logStream << "\t stretch \t" << nonDimStretchEnergy << std::endl;
//                    logStream << "\t bend \t\t" << nonDimBendEnergy << std::endl;
//                    logStream << "\t kinetic \t" << nonDimKineticEnergy << std::endl;
//                    logStream.close();
//                }
//
//                /* If equilibrium reached, write output data to file, and move to next
//                'dialling in' phase. */
//                if (status == EquilibriumReached) {
//                    try {
//                        calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
//                                       interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
//                        if (settings.is_energy_densities_printed) {
//                            calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
//                                                    stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                                    cauchyStressEigenvals, cauchyStressEigenvecs, settings);
//                        }
//                        writeVTKDataOutput(nodes, triangles, stepCount, time, dial_in_factor,
//                                           progTensorSequenceCounter, gaussCurvatures, meanCurvatures, angleDeficits,
//                                           interiorNodeAngleDeficits, boundaryNodeAngleDeficits, stretchEnergyDensities,
//                                           bendEnergyDensities, stretchEnergies, bendEnergies, kineticEnergies,
//                                           strainMeasures, cauchyStressEigenvals, cauchyStressEigenvecs, settings,
//                                           outputDirName);
//                    }
//                    catch (const std::runtime_error &writing_error) {
//                        logStream.open();
//                        logStream << writing_error.what() << std::endl;
//                        logStream.close();
//                        return -1;
//                    }
//                    logStream.open();
//                    logStream << "Equilibrium reached. Wrote VTK output at " << getRealTime() << ", stepCount = "
//                              << stepCount << ", simulation time = " << time << ", current dial-in factor = "
//                              << dial_in_factor << std::endl;
//                    status = Dialling;
//                    timeSinceCurrDiallingInPhaseStarted = 0.0;
//                    phase_counter += 1;
//                    settings.num_damp_factor = settings.dial_in_damping *
//                                               settings.damping_scale; // Set damping factor back to dialling phase value.
//                    if (phase_counter <= dial_in_phases.size() - 2) {
//                        logStream << "New dialling in phase beginning, from value "
//                                  << dial_in_phases[phase_counter] << ", to "
//                                  << dial_in_phases[phase_counter + 1] << std::endl;
//                    }
//                    logStream.close();
//                }
//
//
//            }
//
//
//            /* Advance node positions and velocities using dynamical solver. Error
//            caught here if any node force is suspiciously high, indicating probable
//            'blowing, up', and the code aborts in that case. */
//            try { advanceDynamics(nodes, triangles, settings, logStream); }
//            catch (const std::runtime_error &error) {
//                logStream.open();
//                logStream << "At " << getRealTime() << ", stepCount = " << stepCount << ", time = " << time
//                          << ", dial-in factor = " << dial_in_factor
//                          << " a force was suspiciously large, there is probably a problem. Writing VTK output and then aborting."
//                          << std::endl;
//                logStream.close();
//                try {
//                    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
//                                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
//                    if (settings.is_energy_densities_printed) {
//                        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
//                                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
//                    }
//                    writeVTKDataOutput(nodes, triangles, stepCount, time, dial_in_factor, progTensorSequenceCounter,
//                                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
//                                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
//                                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
//                                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
//                }
//                catch (const std::runtime_error &writing_error) {
//                    logStream.open();
//                    logStream << writing_error.what() << std::endl;
//                    logStream.close();
//                    return -1;
//                }
//                return -1;
//            }
//
//            // Advance times and step counter.
//            time += settings.time_step;
//            timeSinceLastEquilCheck += settings.time_step;
//            timeSinceCurrDiallingInPhaseStarted += settings.time_step;
//            stepCount += 1;
//        }
//
//        if (programmed_metric_infos.size() > 2 && progTensorSequenceCounter < programmed_metric_infos.size() - 2) {
//            logStream.open();
//            logStream << "Moving on to next set of programmed tensors in sequence." << std::endl;
//            logStream.close();
//        }
//    }
//
//// Print some helpful final things.
//    logStream.open();
//    logStream << "Reached simulation time = " << time << " using " << stepCount << " time steps" << std::endl;
//    logStream.close();
//
//    return 0;
//}
