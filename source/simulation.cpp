//
// Created by Michał Zmyślony on 12/02/2024.
//

#define _USE_MATH_DEFINES

#include "simulation.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <libconfig.h++>
#include <chrono>

#include "functions/getRealTime.hpp"
#include "functions/extract_Just_Filename.hpp"
#include "initialDirAndFileHandling.hpp"
#include "readSettingsFile.hpp"
#include "readVTKData.hpp"
#include "calculations/calcTrianglesIncidentOnNodes.hpp"
#include "calculations/calcTriangleAdjacencies_And_Edges.hpp"
#include "functions/kahanSum.hpp"
#include "calculations/calcNodeNeighbours.hpp"
#include "calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp"
#include "calculations/calcTriangleGeometries_and_DialledProgTensors.hpp"
#include "setRemainingInitCond_and_NodeMasses.hpp"
#include "functions/perturbInitialPositionsWithRandomNoise.hpp"
#include "calculations/calcEnergiesAndStresses.hpp"
#include "EquilibriumCheck.hpp"
#include "calculations/calcCurvatures.hpp"
#include "advanceDynamics.hpp"
#include "writeVTKDataOutput.hpp"
#include "calculations/calcDeformationForces.hpp"
#include "calculations/calcNonDeformationForces_and_ImposeBCS.hpp"
#include "functions/zeroForces.hpp"


void Simulation::setup_logstream() {
    log_stream.setOutputFileName(outputDirName + "/log.txt");
// Check we can create a log file with the chosen name.
    try {
        log_stream.open();
    }
    catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
        throw std::runtime_error("Unable to create logfile.");
    }
    log_stream << std::defaultfloat << std::setprecision(6);
    log_stream.close();
    log_filename = log_stream.getOutputFileName();
}


void Simulation::setup_filenames(int argc, char *argv[]) {
    /* Create a string to contain some things that will later be written to the log
file, as long as all goes to plan and the code gets that far. The first thing to
add to this string is a record of what command line arguments where given, and
at what date and time the run began.*/

    init_string += "Simulation run began at: " + getRealTime() + "\n";
    init_string += "The command that was run was:\n";
    for (int i = 0; i < argc; ++i) {
        init_string += argv[i];
        init_string += " ";
    }
    init_string += "\n\n";

/* Check settings and data files given, with ansatz file optional, and get names
from command line. */
    if ((argc != 3) && (argc != 4)) {
        throw std::runtime_error("Error: Please provide at least one settings file path "
                                 "followed by one data file path as command line arguments.\n"
                                 "A path to a file with ansatz node coordinates can be given as a third "
                                 "argument if desired. All paths should be relative to the current working"
                                 "directory, i.e. they should start with ./");
    }
    settings_filename = argv[1];
    initialisation_filename = argv[2];

    std::string settings_file_name_str_final_piece = extract_Just_Filename(settings_filename);
    std::string initial_data_file_name_str_final_piece = extract_Just_Filename(initialisation_filename);
    std::string ansatz_data_file_name_str_final_piece("no_ansatz_file");

    if (argc == 4) {
        ansatz_filename = argv[3];
        ansatz_data_file_name_str_final_piece = extract_Just_Filename(ansatz_filename);
        init_string += "Using ansatz data file " + ansatz_data_file_name_str_final_piece + "\n";
    }


    try {
        outputDirName = initialDirAndFileHandling(settings_filename,
                                                  settings_file_name_str_final_piece,
                                                  initialisation_filename,
                                                  initial_data_file_name_str_final_piece,
                                                  settings_filename,
                                                  ansatz_data_file_name_str_final_piece,
                                                  argc,
                                                  init_string);
    }
    catch (const std::runtime_error &error) {
        std::cerr << error.what() << std::endl;
        std::cerr << "Initial directory and file handling failed." << std::endl;
        throw std::runtime_error("Initial directory and file handling failed.");
    }
    setup_logstream();

    log_stream.open();
    log_stream << init_string << std::endl;
    log_stream.close();
}


void Simulation::read_settings() {
    try { readSettingsFile(settings, settings_filename); }
    catch (const libconfig::FileIOException &fioex) {
        log_stream.open();
        log_stream << "I/O error reading settings file. Check file exists in correct directory and has correct name."
                   << std::endl;
        return;
    }
    catch (const libconfig::ParseException &pex) {
        log_stream << "Parse error reading settings file. Check file has correct format." << std::endl;
        return;
    }
    catch (const std::runtime_error &error) {
        log_stream << error.what() << std::endl;
        log_stream.close();
        return;
    }
}


void Simulation::read_vtk_data() {
    log_stream.open();
    log_stream << "Now attempting to read data files. An error here likely \n"
                  "implies a problem with a data file, for example a mismatch between the \n"
                  "number of nodes or triangles stated and the number actually present; \n"
                  "or other similar mismatches in numbers of data values; or a format problem. \n"
                  "Remember the input files must have exactly the correct format." << std::endl;
    log_stream.close();
    try {
        readVTKData(nodes, triangles, programmed_metric_infos, inverted_programmed_metrics, programmed_taus,
                    programmed_second_fundamental_forms, settings, initialisation_filename,
                    progTensorSequenceCounterToStartFrom,
                    dialInFactorToStartFrom, nodeAnsatzPositions, ansatz_filename, log_stream);
    }
    catch (const std::out_of_range &out_of_bounds_error) {
        log_stream << out_of_bounds_error.what() << std::endl;
        return;
    }
    catch (const std::runtime_error &reading_error) {
        log_stream << reading_error.what() << std::endl;
        return;
    }
}


void Simulation::configure_nodes() {
    log_stream.open();
    log_stream << "Number of nodes = " << settings.NumNodes << std::endl;
    log_stream << "Number of triangles = " << settings.NumTriangles << std::endl;

    /* Print warning if number of triangles is low - in this case boundary effects
    will dominate, and the code should not be trusted, both due to the reduced
    accuracy in the treatment of the boundary in e.g. the 2nd F.F. approx, and due
    to possible bugs in this rather special case.*/
    if (settings.NumTriangles < 50) {
        log_stream << "Your mesh has a small number of triangles. \nBeware that the code "
                      "is likely to be less accurate in this case, \nand unforeseen bugs are more "
                      "likely in extreme cases." << std::endl;
    }
    log_stream.close();


    for (int i = 0; i < settings.NumNodes; ++i) {
        nodes[i].nodeLogStream.setOutputFileName(log_filename);
    }
}


void Simulation::configure_topological_properties() {
    calcTrianglesIncidentOnNodes(nodes, triangles, settings);
    calcTriangleAdjacencies_And_Edges(nodes, triangles, edges, settings);

    for (int i = 0; i < settings.NumEdges; ++i) {
        edges[i].edgeLogStream.setOutputFileName(log_filename);
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
    log_stream.open();
    log_stream << "Initial perimeter = " << initPerimeter << std::endl;
    log_stream.close();


//A further check that things are ok:
    if (3 * settings.NumTriangles != 2 * settings.NumEdges - numBoundaryEdges) {
        throw std::runtime_error(
                "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
    }

    log_stream.open();
    log_stream << "Number of edges = " << settings.NumEdges << std::endl;
    log_stream << "Number of boundary edges = " << numBoundaryEdges << std::endl;
    log_stream << "Number of non-boundary edges = " << settings.NumEdges - numBoundaryEdges << std::endl;
    log_stream.close();
}


void Simulation::configure_triangles() {
    for (int i = 0; i < settings.NumTriangles; ++i) {
        triangles[i].triLogStream.setOutputFileName(log_filename);
    }

    int numBoundaryTriangles = 0;
    for (int i = 0; i < settings.NumTriangles; ++i) {
        if (triangles[i].isOnBoundary) {
            numBoundaryTriangles += 1;
        }
    }

    log_stream.open();
    log_stream << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
    log_stream << "Number of holes in mesh = " << 1 + settings.NumEdges - settings.NumNodes - settings.NumTriangles
               << std::endl; // From Euler's formula for a planar graph.
    log_stream.close();

    if (1 + settings.NumEdges - settings.NumNodes - settings.NumTriangles < 0) {
        throw std::runtime_error(
                "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
    }
}


void Simulation::count_boundary_nodes() {
    int numBoundaryNodes = 0;
    for (int n = 0; n < settings.NumNodes; ++n) {
        if (nodes[n].isOnBoundary) {
            numBoundaryNodes += 1;
        }
    }
    settings.numBoundaryNodes = numBoundaryNodes;
}


void Simulation::set_node_patches() {
    log_stream.open();
    log_stream << "\n" << "Beginning patch selection and related pre-calculations." << std::endl;
    log_stream.close();

    calc_nonVertexPatchNodes_and_MatForPatchDerivs(nodes, triangles, settings, log_stream);

    log_stream.open();
    log_stream << "Successfully completed patch setup." << "\n" << std::endl;
    log_stream.close();
}


void Simulation::orient_node_labels() {
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, simulation_status, -12345, 98765,
                                                  programmed_metric_infos, inverted_programmed_metrics,
                                                  programmed_taus, programmed_second_fundamental_forms, settings);
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
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, simulation_status, -12345, 98765,
                                                  programmed_metric_infos, inverted_programmed_metrics,
                                                  programmed_taus, programmed_second_fundamental_forms, settings);
}


void Simulation::set_initial_conditions() {
    setRemainingInitCond_and_NodeMasses(nodes, triangles, edges, programmed_metric_infos,
                                        inverted_programmed_metrics, programmed_taus,
                                        programmed_second_fundamental_forms,
                                        settings);
}

void Simulation::find_smallest_element() {
    settings.ApproxMinInitElemSize = DBL_MAX;
    settings.smallestSizeOverRootTau = DBL_MAX;
    for (int i = 0; i < triangles.size(); ++i) {

        // This variable will end up holding the smallest altitude for *this* tri.
        double smallestAltitude =
                2 * triangles[i].initArea / (triangles[i].currSides.col(0) - triangles[i].currSides.col(1)).norm();

        for (int s = 0; s < 2; ++s) {
            if (smallestAltitude > 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm()) {
                smallestAltitude = 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm();
            }
        }

        if (settings.ApproxMinInitElemSize > smallestAltitude) {
            settings.ApproxMinInitElemSize = smallestAltitude;
        }

        for (auto &sequenceOf_ProgTau: programmed_taus) {
            if (settings.smallestSizeOverRootTau > smallestAltitude / sqrt(sequenceOf_ProgTau[i])) {
                settings.smallestSizeOverRootTau = smallestAltitude / sqrt(sequenceOf_ProgTau[i]);
            }
        }
    }
    log_stream.open();
    log_stream << "Sheet thickness = " << settings.SheetThickness << std::endl;
    log_stream << "Approx smallest element linear size = " << settings.ApproxMinInitElemSize << std::endl;
    log_stream.close();
}

void Simulation::print_total_load_force() {
    int numLoadedNodes = 0;
    for (int i = 0; i < settings.NumNodes; ++i) {
        if (nodes[i].isLoadForceEnabled) {
            numLoadedNodes += 1;
        }
    }
    log_stream.open();
    log_stream << "Total load force applied = " <<
               numLoadedNodes * settings.LoadStrength * settings.ShearModulus * settings.ApproxMinInitElemSize *
               settings.SheetThickness << std::endl;
    log_stream.close();
}

void Simulation::setup_characteristic_scales() {
    settings.TimeBetweenEquilChecks = settings.TimeBetweenEquilChecksPrefactor * settings.DialInStepTime;

    settings.charForceScale = settings.ShearModulus * settings.ApproxMinInitElemSize * settings.SheetThickness;
    settings.charStretchEnergyDensityScale = settings.ShearModulus * settings.SheetThickness;
// settings.charBendEnergyDensityScale = settings.ShearModulus * settings.SheetThickness * settings.SheetThickness * settings.SheetThickness / (settings.SampleCharLength * settings.SampleCharLength);
    double totInitArea;
    std::vector<double> initAreas(settings.NumTriangles);
    for (int i = 0; i < settings.NumTriangles; ++i) {
        initAreas[i] = triangles[i].initArea;
    }
    totInitArea = kahanSum(initAreas);
    settings.charStretchEnergyScale = settings.charStretchEnergyDensityScale * totInitArea;
// settings.charBendEnergyScale = settings.charBendEnergyDensityScale * totInitArea;
}


void Simulation::setup_equilibrium_dial_in_factors() {
    double tempDialInFactor = 0.0;
    if (settings.DialInResolution > 0 && settings.DialInStepTime >= 0) {
        while (tempDialInFactor < 1.0) {
            DialInFactorValuesToHoldAt.push_back(tempDialInFactor);
            tempDialInFactor += settings.DialInResolution;
        }
        DialInFactorValuesToHoldAt.push_back(1.0);
    } else {
        log_stream.open();
        log_stream << "settings.DialInResolution and settings.DialInStepTime must "
                      "both be >0 and >=0 respectively, and they aren't currently. Aborting."
                   << std::endl;
        log_stream.close();
        throw std::runtime_error("settings.DialInResolution and settings.DialInStepTime must "
                                 "both be >0 and >=0 respectively, and they aren't currently. Aborting.");
    }
}

void Simulation::initialise_simulation_vectors() {
    gaussCurvatures = std::vector<double>(settings.NumTriangles, DBL_MAX);
    meanCurvatures = std::vector<double>(settings.NumTriangles, DBL_MAX);
    stretchEnergies = std::vector<double>(settings.NumTriangles, DBL_MAX);
    bendEnergies = std::vector<double>(settings.NumTriangles, DBL_MAX);
    stretchEnergyDensities = std::vector<double>(settings.NumTriangles, DBL_MAX);
    bendEnergyDensities = std::vector<double>(settings.NumTriangles, DBL_MAX);
    kineticEnergies = std::vector<double>(settings.NumNodes, DBL_MAX);
    strainMeasures = std::vector<double>(settings.NumTriangles, DBL_MAX);
    cauchyStressEigenvals = std::vector<Eigen::Vector2d>(settings.NumTriangles);
    cauchyStressEigenvecs = std::vector<Eigen::Matrix<double, 3, 2>>(settings.NumTriangles);

    angleDeficits = std::vector<double>(settings.NumNodes, DBL_MAX);
    interiorNodeAngleDeficits = std::vector<double>(settings.NumNodes, DBL_MAX);
    boundaryNodeAngleDeficits = std::vector<double>(settings.NumNodes, DBL_MAX);
}


void Simulation::init(int argc, char *argv[]) {
    setup_filenames(argc, argv);
    read_settings();
    read_vtk_data();

    configure_nodes();
    configure_topological_properties();
    configure_triangles();
    configureNodeAdjacency(nodes, edges, settings);
    count_boundary_nodes();
    set_node_patches();
    orient_node_labels();

    set_initial_conditions();
    find_smallest_element();

    settings.SetupDialIn(log_stream);
    settings.SetupStepTime(log_stream);
    settings.SetupPrintFrequency(log_stream);

    print_total_load_force();
    setup_characteristic_scales();
    setup_equilibrium_dial_in_factors();

    initialise_simulation_vectors();

    if (settings.isPerturbationOfInitialPositionsEnabled) {
        perturbInitialPositionsWithRandomNoise(nodes, settings);
    }
}


void Simulation::run_ansatz(int counter) {
    for (int i = 0; i < settings.NumNodes; ++i) {
        nodes[i].pos = nodeAnsatzPositions[i];
    }

    timeSinceLastEquilCheck = 0;
    simulation_status = dialling;
    DialInFactorCounter = 0;

    if (!settings.isDialingFromAnsatzEnabled) {

        currDialInFactor = dialInFactorToStartFrom;
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
                programmed_second_fundamental_forms[progTensorSequenceCounterToStartFrom][i] = programmed_second_fundamental_forms[
                        progTensorSequenceCounterToStartFrom + 1][i];
                programmed_metric_infos[progTensorSequenceCounterToStartFrom][i] = programmed_metric_infos[
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
        log_stream.open();
        log_stream
                << "\nsettings.isDialingFromAnsatzEnabled == true, FIRST DIALLING IN PHASE WILL BE FROM ANSATZ STATE."
                << std::endl;
        log_stream.close();

        currDialInFactor = 0.0;
        timeSinceCurrDiallingInPhaseStarted = 0.0;

        // Calculate all necessary geometry for the ansatz state.
        calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, waitingForEquilibrium, currDialInFactor,
                                                      counter, programmed_metric_infos,
                                                      inverted_programmed_metrics, programmed_taus,
                                                      programmed_second_fundamental_forms, settings);
        updateSecondFundamentalForms(triangles, settings);


        // Temp LU decomp of triangle metric, used to check invertibility.
        Eigen::FullPivLU<Eigen::Matrix<double, 2, 2>> tempMetricDecomp;


        /* Alter inverted_programmed_metrics[progTensorSequenceCounterToStartFrom] and similar to change where
        the programmed quantities are dialling from.*/

        for (int i = 0; i < settings.NumTriangles; ++i) {

            // Test for invertibility of metric before taking inverse.
            tempMetricDecomp.compute(
                    (triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse());
            if (!tempMetricDecomp.isInvertible()) {
                throw std::runtime_error(
                        "At least one triangle had a non-invertible metric in the ansatz state. \n"
                        "This should not occur in a reasonable mesh. Aborting.");
            } else {
                inverted_programmed_metrics[progTensorSequenceCounterToStartFrom][i] = (
                        triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                if (settings.isDialingDisabled) {
                    inverted_programmed_metrics[progTensorSequenceCounterToStartFrom + 1][i] = (
                            triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                }
            }


            // Other quantities require no inverse, so are easier.
            programmed_taus[progTensorSequenceCounterToStartFrom][i] = programmed_taus[
                    progTensorSequenceCounterToStartFrom + 1][i];
            programmed_second_fundamental_forms[progTensorSequenceCounterToStartFrom][i] = triangles[i].secFF;
            if (settings.isDialingDisabled) {
                programmed_taus[progTensorSequenceCounterToStartFrom + 1][i] = programmed_taus[
                        progTensorSequenceCounterToStartFrom + 1][i];
                programmed_second_fundamental_forms[progTensorSequenceCounterToStartFrom +
                                                    1][i] = triangles[i].secFF;
            }
        }
    }

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
}

void Simulation::setup_imposed_seide_deformations(double &s1, int highest_node, int lowest_node,
                                                  std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    settings.lambda = 0.9;
    settings.ConeAngle = asin(pow(settings.lambda, 1.5));
    log_stream.open();
    log_stream << "IMPOSING SEIDE DEFORMATIONS." << std::endl;
    log_stream.close();

    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
    // positions we will force to be those of Seide's setup (found
    // by findng the membrane stress solution for his base state, and
    // integrating the strains etc).
    double intermLengthScaleUpper = 3.0 * sqrt(settings.SheetThickness * 0.18);
    double intermLengthScaleLower = 3.0 * sqrt(settings.SheetThickness * 1.8);
    log_stream.open();
    log_stream
            << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
    log_stream.close();
    for (int n = 0; n < settings.NumNodes; ++n) {
        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
            nodes[n].isSeideDisplacementEnabled = true;
        }
    }

    // STORE PERFECT CONE ANSATZ
    for (int n = 0; n < settings.NumNodes; ++n) {
        nodeUnstressedConePosits[n] = nodes[n].pos;
    }
    s1 = sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
              nodes[highest_node].pos(1) * nodes[highest_node].pos(1)) / sin(settings.ConeAngle);

//    Eigen::Vector3d testTriCurrCentroid;
//    testTriCurrCentroid = (nodes[triangles[settings.testTriangle].vertexLabels(0)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(1)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(2)].pos) / 3;
    //sTest = sqrt(testTriCurrCentroid(0)*testTriCurrCentroid(0) + testTriCurrCentroid(1)*testTriCurrCentroid(1)) / sin(settings.ConeAngle);

}

void Simulation::setup_glass_cones(int highest_node, int lowest_node) {
    settings.ConeAngle = 1.02327019;
    settings.initSlideZCoord_upper += -tan(settings.ConeAngle) *
                                      sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
                                           nodes[highest_node].pos(1) * nodes[highest_node].pos(1));
    settings.initSlideZCoord_lower += -tan(settings.ConeAngle) *
                                      sqrt(nodes[lowest_node].pos(0) * nodes[lowest_node].pos(0) +
                                           nodes[lowest_node].pos(1) * nodes[lowest_node].pos(1));
    log_stream.open();
    log_stream << "USING TWO GLASS CONES FOR SQUASHING." << std::endl;
    log_stream.close();

    // Hijack nodes[i].isOnBoundary to instead label nodes whose
    // forces we will modify to kill any components not tangential to
    // a perfect cone base state. We do this to nodes within an intermediate
    // distance vertically from the ends of the cone.
    double intermLengthScaleUpper = 2.0 * sqrt(settings.SheetThickness * 0.18);
    double intermLengthScaleLower = 2.0 * sqrt(settings.SheetThickness * 1.8);
    log_stream.open();
    log_stream
            << "Apply normal-force-killer within the following vertical distances of the top and bottom: "
            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
    log_stream.close();
    for (int n = 0; n < settings.NumNodes; ++n) {
        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
            nodes[n].isOnBoundary = true;
        }
    }
}

void Simulation::update_slide_properties() {
    if (!settings.isControlledForceEnabled) {
        settings.upperSlideDisplacement =
                time * settings.slideSpeedPrefactor * settings.SampleCharLength / settings.bending_long_time;
    } else { // settings.isControlledForceEnabled == true instead
        //settings.upperSlideWeight = (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) * (time * settings.slideSpeedPrefactor / bending_long_time);
        settings.slideDampingParam =
                0.4 * settings.ShearModulus * settings.SheetThickness * settings.SheetThickness /
                (settings.slideSpeedPrefactor * settings.SampleCharLength / settings.bending_long_time);
        if (fabs(settings.upperTotSlideForce + settings.upperSlideWeight) /
            (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) <
            settings.totalSlideForceToMuTSqRatioEquilThreshold
            && settings.constSlideWeightFac < 0) {
            if (timeSinceLastEquilCheck > settings.TimeBetweenEquilChecks) {
                if (equilibriumCheck(nodes, triangles, settings, log_stream) == equilibriumReached) {
                    settings.slideJustReachedEquil = 1;
                    settings.upperSlideWeight +=
                            settings.slideWeightDialSpeedFac *
                            (settings.TimeStep / settings.bending_long_time) *
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
}


void Simulation::impose_seide_deformation(double s1, const std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    if (settings.isSeideDeformationsEnabled) {
        double pInit = 3.0 * (settings.ShearModulus * settings.SheetThickness *
                              settings.SheetThickness); // So we don't have to start all the way from p=0, chosen based on previous sims.
        settings.p = pInit + time * settings.pSpeedPrefactor * settings.ShearModulus * settings.SheetThickness *
                             settings.SheetThickness / settings.bending_long_time;

        for (int n = 0; n < settings.NumNodes; ++n) {
            if (nodes[n].isSeideDisplacementEnabled || step_count == 0) {

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
                           (2.0 * M_PI * settings.YoungsModulus * tCone * sin(settings.ConeAngle) *
                            cos(settings.ConeAngle));
                double w = -settings.p * (log(s / s1) + settings.PoissonRatio) /
                           (2.0 * M_PI * settings.YoungsModulus * tCone * cos(settings.ConeAngle) *
                            cos(settings.ConeAngle));
                nodes[n].pos = nodeUnstressedConePosits[n] + u * uHat + w * wHat;
            }
        }
    }
}

void Simulation::first_step_configuration(double &seide_quotient,
                                          std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    int highest_node = INT_MIN;
    int lowest_node = INT_MIN;
    settings.initSlideZCoord_lower = nodes[0].pos(2);
    settings.initSlideZCoord_upper = nodes[0].pos(2);
    for (int n = 0; n < settings.NumNodes; ++n) {
        if (settings.initSlideZCoord_lower > nodes[n].pos(2)) {
            settings.initSlideZCoord_lower = nodes[n].pos(2);
            lowest_node = n;
        }
        if (settings.initSlideZCoord_upper < nodes[n].pos(2)) {
            settings.initSlideZCoord_upper = nodes[n].pos(2);
            highest_node = n;
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
        setup_glass_cones(highest_node, lowest_node);
    }

    // Instead set up imposed-Seide-deformations idea.
    if (settings.isSeideDeformationsEnabled) {
        setup_imposed_seide_deformations(seide_quotient, highest_node, lowest_node, nodeUnstressedConePosits);
    }
}

void Simulation::begin_equilibrium_search(int counter) {
    currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter + 1];
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, simulation_status, currDialInFactor,
                                                  counter, programmed_metric_infos,
                                                  inverted_programmed_metrics, programmed_taus,
                                                  programmed_second_fundamental_forms, settings);
    simulation_status = waitingForEquilibrium;
    // NB an EquilCheck has not actually just occurred, but this has the
    //desired effect of ensuring that each DialInFactor value is held
    //for at least one TimeBetweenEquilChecks.
    timeSinceLastEquilCheck = 0.0;

    settings.NumDampFactor =
            settings.DampingPrefactor2 *
            settings.dampingScale; // Set damping factor to waiting phase value.

    log_stream.open();
    log_stream << "Reached dial-in factor of " << DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
               << ". Waiting for equilibrium." << std::endl;
    log_stream.close();
}

void Simulation::progress_single_step(int counter, std::pair<double, double> upperAndLowerTotSlideForces,
                                      std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes) {
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, simulation_status, currDialInFactor,
                                                  counter, programmed_metric_infos,
                                                  inverted_programmed_metrics, programmed_taus,
                                                  programmed_second_fundamental_forms, settings);



    /* Calculate secFF estimates for triangles, and related quantities such
    as the derivative of the bending energy wrt the secFF components.*/
    updateSecondFundamentalForms(triangles, settings);

    // Calculate current strain and bending force on each node.
    calcDeformationForces(nodes, triangles, settings, correspondingTrianglesForNodes);

    /* Add force contributions from e.g. damping, loads, 'prod' perturbation, and
     account for BCs e.g. clamping. */
    upperAndLowerTotSlideForces = calcNonDeformationForces_and_ImposeBCS(nodes, time, settings);
    settings.upperTotSlideForce = upperAndLowerTotSlideForces.first;
}

void Simulation::update_dial_in_factor() {
    currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter] +
                       (DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                        - DialInFactorValuesToHoldAt[DialInFactorCounter]) *
                       timeSinceCurrDiallingInPhaseStarted / settings.DialInStepTime;
}


void Simulation::save_and_print_details(int counter, double duration_us,
                                        std::pair<double, double> upperAndLowerTotSlideForces) {
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
    writeVTKDataOutput(nodes, triangles, step_count, time, currDialInFactor, counter,
                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);

    log_stream.open();
    log_stream << std::fixed << "Wrote VTK output at " << getRealTime() << ", stepCount = " << step_count
               << ", simulation time = " << time + settings.TimeStep << ", current dial-in factor = "
               << currDialInFactor << std::scientific;
    log_stream << ", last step's execution time "
               << duration_us << " us"
               << std::endl;

    log_stream.close();

    std::ofstream forceDistFile;
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
                  << std::endl;
    forceDistFile.close();
}

std::stringstream Simulation::log_prefix() {
    std::stringstream ss;
    ss << getRealTime() << ", step = " << step_count << ", time = " << time;
    return ss;
}

void Simulation::error_large_force(int counter) {
    log_stream.open();
    log_stream << "At " << getRealTime() << ", stepCount = " << step_count << ", time = " << time
               << ", dial-in factor = " << currDialInFactor
               << " a force was suspiciously large, there is probably a problem. Writing VTK output and then aborting."
               << std::endl;
    log_stream.close();

    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
    if (settings.isEnergyDensitiesPrinted) {
        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
    }
    writeVTKDataOutput(nodes, triangles, step_count, time, currDialInFactor, counter,
                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
    throw std::runtime_error("Unexpectedly large force.");
}

void Simulation::check_for_equilibrium() {
    log_stream.open();
    log_stream << "Checking for equilibrium at " << getRealTime() << ", stepCount = " << step_count
               << ", simulation time = " << time << ", current dial-in factor = " << currDialInFactor
               << std::endl;
    log_stream.close();

    simulation_status = equilibriumCheck(nodes, triangles, settings, log_stream);
    timeSinceLastEquilCheck = 0.0;

    calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                            stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                            cauchyStressEigenvals, cauchyStressEigenvecs, settings);
//    double nonDimStretchEnergy = kahanSum(stretchEnergies) / settings.charStretchEnergyScale;
//    double nonDimBendEnergy = kahanSum(bendEnergies) / settings.charStretchEnergyScale;
//    double nonDimKineticEnergy = kahanSum(kineticEnergies) / settings.charStretchEnergyScale;

//    log_stream.open();
//    log_stream << "Non-dimensionalised energies:" << std::endl;
//    log_stream << "\t total \t\t" << nonDimStretchEnergy + nonDimBendEnergy + nonDimKineticEnergy
//               << std::endl;
//    log_stream << "\t stretch \t" << nonDimStretchEnergy << std::endl;
//    log_stream << "\t bend \t\t" << nonDimBendEnergy << std::endl;
//    log_stream << "\t kinetic \t" << nonDimKineticEnergy << std::endl;
//    log_stream.close();
}

void Simulation::setup_reached_equilibrium() {
    log_stream.open();
    simulation_status = dialling;
    timeSinceCurrDiallingInPhaseStarted = 0.0;
    DialInFactorCounter += 1;
    settings.NumDampFactor = settings.DampingPrefactor1 *
                             settings.dampingScale; // Set damping factor back to dialling phase value.
    if (DialInFactorCounter <= DialInFactorValuesToHoldAt.size() - 2) {
        log_stream << "New dialling in phase beginning, from value "
                   << DialInFactorValuesToHoldAt[DialInFactorCounter] << ", to "
                   << DialInFactorValuesToHoldAt[DialInFactorCounter + 1] << std::endl;
    }
    log_stream.close();
}

void Simulation::run_tensor_increment(int counter) {
    /* A third input file specifying an ansatz for the node positions may have
        been read in if it was given as a command line argument. This could for
        example correspond to an output file from this code, to carry on where some
        previous simulation left off. In this case, the nodes are moved to their
        ansatz positions here, and other relevant variables are set up. */
    if (ansatz_filename != "no_ansatz_file" &&
        counter == progTensorSequenceCounterToStartFrom) {
        run_ansatz(counter);
    } else {
        // Reset the variables that control each dialling/waiting process.
        timeSinceLastEquilCheck = 0;
        simulation_status = dialling;
        currDialInFactor = 0;
        DialInFactorCounter = 0;
        timeSinceCurrDiallingInPhaseStarted = 0;
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
    log_stream.open();
    log_stream << "\nBeginning dynamical evolution.\n" << std::endl;
    log_stream.close();

    log_stream.open();
    log_stream << "\nCREATING VECTOR TO STORE UNSTRESSED CONE NODE POSITIONS.\n" << std::endl;
    log_stream.close();
    std::vector<Eigen::Vector3d> nodeUnstressedConePosits(settings.NumNodes);
    double seide_quotient = DBL_MAX;

    std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes = getCorrespondingTrianglesForNodes(
            triangles, nodes);

    while (DialInFactorCounter <= DialInFactorValuesToHoldAt.size() - 2) {
        if (step_count == 0) { first_step_configuration(seide_quotient, nodeUnstressedConePosits); }

        update_slide_properties();

        if (settings.isSeideDeformationsEnabled) {
            impose_seide_deformation(seide_quotient, nodeUnstressedConePosits);
        }

        zeroForces(nodes);

        // Check if still in dialling in phase or whether it is time to wait for
        // equilibrium.
        if (timeSinceCurrDiallingInPhaseStarted >= settings.DialInStepTime && simulation_status == dialling) {
            begin_equilibrium_search(counter);
        }

        /* If not waiting for equilibrium, set the current value of the Dial-In
        Factor, based on linear dialling in between the previously calculated
        'checkpoint' values. */
        if (simulation_status == dialling) { update_dial_in_factor(); }

        std::pair<double, double> upperAndLowerTotSlideForces;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        progress_single_step(counter, upperAndLowerTotSlideForces, correspondingTrianglesForNodes);
        auto duration = std::chrono::steady_clock::now() - begin;
        double duration_us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();


        /* Write output data regularly.
        Can be switched off with settings.PrintFrequency < 0.0.
        Doing the write-out at this point in the loop means the node positions
        and the triangle geometry data match in the output, which is desirable!*/
        if ((step_count % settings.InversePrintRate == 0 && settings.InversePrintRate > 0) ||
            settings.slideJustReachedEquil == 1) {
            save_and_print_details(counter, duration_us, upperAndLowerTotSlideForces);
            settings.slideJustReachedEquil = 0;
        }


        if (!settings.isControlledForceEnabled) {
            /* If last check for equilibrium or last reaching of a new DialInFactor
             value was more than TimeBetweenEquilChecks ago, check for equilibrium.
             Also print total stretching and bending energies.*/
            if (timeSinceLastEquilCheck > settings.TimeBetweenEquilChecks &&
                simulation_status == waitingForEquilibrium) {
                check_for_equilibrium();
            }

            /* If equilibrium reached, write output data to file, and move to next
            'dialling in' phase. */
            if (simulation_status == equilibriumReached) {
                save_and_print_details(counter, duration_us, upperAndLowerTotSlideForces);
                setup_reached_equilibrium();
            }
        }

        /* Advance node positions and velocities using dynamical solver. Error
        caught here if any node force is suspiciously high, indicating probable
        'blowing, up', and the code aborts in that case. */
        try { advanceDynamics(nodes, triangles, settings, log_stream); }
        catch (const std::runtime_error &error) {
            error_large_force(counter);
        }

        // Advance times and step counter.
        time += settings.TimeStep;
        timeSinceLastEquilCheck += settings.TimeStep;
        timeSinceCurrDiallingInPhaseStarted += settings.TimeStep;
        step_count += 1;
    }
}


void Simulation::run_simulation() {
    log_stream << std::scientific << std::setprecision(8);
    std::ofstream forceDistFile;

    /* Loop over the sequence of programmed tensors, dialling-in and waiting for
    equilibrium between each pair in the sequence. This loop is redundant in most use
    cases for this code, where only a single set of programmed tensors is supplied.*/
    for (std::size_t progTensorSequenceCounter = progTensorSequenceCounterToStartFrom;
         progTensorSequenceCounter <= programmed_metric_infos.size() - 2; ++progTensorSequenceCounter) {
        run_tensor_increment(progTensorSequenceCounter);

        if (programmed_metric_infos.size() > 2 && progTensorSequenceCounter < programmed_metric_infos.size() - 2) {
            log_stream.open();
            log_stream << "Moving on to next set of programmed tensors in sequence." << std::endl;
            log_stream.close();
        }
    }

// Print some helpful final things.
    log_stream.open();
    log_stream << "Reached simulation time = " << time << " using " << step_count << " time steps" << std::endl;
    log_stream.close();
}

Simulation::Simulation(int argc, char **argv) {
    init(argc, argv);
}
