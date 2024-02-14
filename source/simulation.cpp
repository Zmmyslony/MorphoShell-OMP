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
    readSettingsFile(settings, settings_filename);
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
    log_stream << "Number of nodes = " << settings.num_nodes << std::endl;
    log_stream << "Number of triangles = " << settings.num_triangles << std::endl;

    /* Print warning if number of triangles is low - in this case boundary effects
    will dominate, and the code should not be trusted, both due to the reduced
    accuracy in the treatment of the boundary in e.g. the 2nd F.F. approx, and due
    to possible bugs in this rather special case.*/
    if (settings.num_triangles < 50) {
        log_stream << "Your mesh has a small number of triangles. \nBeware that the code "
                      "is likely to be less accurate in this case, \nand unforeseen bugs are more "
                      "likely in extreme cases." << std::endl;
    }
    log_stream.close();


    for (int i = 0; i < settings.num_nodes; ++i) {
        nodes[i].nodeLogStream.setOutputFileName(log_filename);
    }
}


void Simulation::configure_topological_properties() {
    calcTrianglesIncidentOnNodes(nodes, triangles, settings);
    calcTriangleAdjacencies_And_Edges(nodes, triangles, edges, settings);

    for (int i = 0; i < settings.num_edges; ++i) {
        edges[i].edgeLogStream.setOutputFileName(log_filename);
    }

    /* Calculate and print the number of boundary and non-boundary edges.
    Also label the nodes connected by boundary edges as boundary nodes.
    Also calculate the total initial perimeter of the sample, to be used as a
    characteristic sample length in estimating characteristic times, time steps
    etc.*/
    int numBoundaryEdges = 0;
    std::vector<double> initBoundaryEdgeLengths(settings.num_edges);

    for (int i = 0; i < settings.num_edges; ++i) {
        if (edges[i].isOnBoundary) {

            numBoundaryEdges += 1;

            nodes[edges[i].nodeLabels(0)].isOnBoundary = true;
            nodes[edges[i].nodeLabels(1)].isOnBoundary = true;

            //If chosen in settings, clamp whole boundary in addition to clamp
            //indicators from data file.
            if (settings.is_boundary_clamped) {
                nodes[edges[i].nodeLabels(0)].isClamped = true;
                nodes[edges[i].nodeLabels(1)].isClamped = true;
            }

            // Store initial length of this boundary edge
            initBoundaryEdgeLengths[i] = (nodes[edges[i].nodeLabels(0)].pos - nodes[edges[i].nodeLabels(1)].pos).norm();
        }
    }

// Do sum to calculate perimeter, and set the characteristic sample length to it.
    double initPerimeter = kahanSum(initBoundaryEdgeLengths);
    settings.sample_char_length = initPerimeter;
    log_stream.open();
    log_stream << "Initial perimeter = " << initPerimeter << std::endl;
    log_stream.close();


//A further check that things are ok:
    if (3 * settings.num_triangles != 2 * settings.num_edges - numBoundaryEdges) {
        throw std::runtime_error(
                "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
    }

    log_stream.open();
    log_stream << "Number of edges = " << settings.num_edges << std::endl;
    log_stream << "Number of boundary edges = " << numBoundaryEdges << std::endl;
    log_stream << "Number of non-boundary edges = " << settings.num_edges - numBoundaryEdges << std::endl;
    log_stream.close();
}


void Simulation::configure_triangles() {
    for (int i = 0; i < settings.num_triangles; ++i) {
        triangles[i].triLogStream.setOutputFileName(log_filename);
    }

    int numBoundaryTriangles = 0;
    for (int i = 0; i < settings.num_triangles; ++i) {
        if (triangles[i].isOnBoundary) {
            numBoundaryTriangles += 1;
        }
    }

    log_stream.open();
    log_stream << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
    log_stream << "Number of holes in mesh = " << 1 + settings.num_edges - settings.num_nodes - settings.num_triangles
               << std::endl; // From Euler's formula for a planar graph.
    log_stream.close();

    if (1 + settings.num_edges - settings.num_nodes - settings.num_triangles < 0) {
        throw std::runtime_error(
                "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
    }
}


void Simulation::count_boundary_nodes() {
    int numBoundaryNodes = 0;
    for (int n = 0; n < settings.num_nodes; ++n) {
        if (nodes[n].isOnBoundary) {
            numBoundaryNodes += 1;
        }
    }
    settings.num_boundary_nodes = numBoundaryNodes;
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
    for (int i = 0; i < settings.num_triangles; ++i) {
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
    settings.approx_min_init_elem_size = DBL_MAX;
    settings.smallest_size_over_root_tau = DBL_MAX;
    for (int i = 0; i < triangles.size(); ++i) {

        // This variable will end up holding the smallest altitude for *this* tri.
        double smallestAltitude =
                2 * triangles[i].initArea / (triangles[i].currSides.col(0) - triangles[i].currSides.col(1)).norm();

        for (int s = 0; s < 2; ++s) {
            if (smallestAltitude > 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm()) {
                smallestAltitude = 2 * triangles[i].initArea / triangles[i].currSides.col(s).norm();
            }
        }

        if (settings.approx_min_init_elem_size > smallestAltitude) {
            settings.approx_min_init_elem_size = smallestAltitude;
        }

        for (auto &sequenceOf_ProgTau: programmed_taus) {
            if (settings.smallest_size_over_root_tau > smallestAltitude / sqrt(sequenceOf_ProgTau[i])) {
                settings.smallest_size_over_root_tau = smallestAltitude / sqrt(sequenceOf_ProgTau[i]);
            }
        }
    }
    log_stream.open();
    log_stream << "Sheet thickness = " << settings.sheet_thickness << std::endl;
    log_stream << "Approx smallest element linear size = " << settings.approx_min_init_elem_size << std::endl;
    log_stream.close();
}

void Simulation::print_total_load_force() {
    int numLoadedNodes = 0;
    for (int i = 0; i < settings.num_nodes; ++i) {
        if (nodes[i].isLoadForceEnabled) {
            numLoadedNodes += 1;
        }
    }
    log_stream.open();
    log_stream << "Total load force applied = " <<
               numLoadedNodes * settings.load_strength * settings.shear_modulus * settings.approx_min_init_elem_size *
               settings.sheet_thickness << std::endl;
    log_stream.close();
}

void Simulation::setup_characteristic_scales() {
    settings.time_between_equil_checks = settings.time_between_equil_checks_prefactor * settings.dial_in_step_time;

    settings.char_force_scale = settings.shear_modulus * settings.approx_min_init_elem_size * settings.sheet_thickness;
    settings.char_stretch_energy_density_scale = settings.shear_modulus * settings.sheet_thickness;
// settings.charBendEnergyDensityScale = settings.ShearModulus * settings.SheetThickness * settings.SheetThickness * settings.SheetThickness / (settings.SampleCharLength * settings.SampleCharLength);
    double totInitArea;
    std::vector<double> initAreas(settings.num_triangles);
    for (int i = 0; i < settings.num_triangles; ++i) {
        initAreas[i] = triangles[i].initArea;
    }
    totInitArea = kahanSum(initAreas);
    settings.char_stretch_energy_scale = settings.char_stretch_energy_density_scale * totInitArea;
// settings.charBendEnergyScale = settings.charBendEnergyDensityScale * totInitArea;
}


void Simulation::setup_equilibrium_dial_in_factors() {
    double tempDialInFactor = 0.0;
    if (settings.dial_in_resolution > 0 && settings.dial_in_step_time >= 0) {
        while (tempDialInFactor < 1.0) {
            DialInFactorValuesToHoldAt.push_back(tempDialInFactor);
            tempDialInFactor += settings.dial_in_resolution;
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
    gaussCurvatures = std::vector<double>(settings.num_triangles, DBL_MAX);
    meanCurvatures = std::vector<double>(settings.num_triangles, DBL_MAX);
    stretchEnergies = std::vector<double>(settings.num_triangles, DBL_MAX);
    bendEnergies = std::vector<double>(settings.num_triangles, DBL_MAX);
    stretchEnergyDensities = std::vector<double>(settings.num_triangles, DBL_MAX);
    bendEnergyDensities = std::vector<double>(settings.num_triangles, DBL_MAX);
    kineticEnergies = std::vector<double>(settings.num_nodes, DBL_MAX);
    strainMeasures = std::vector<double>(settings.num_triangles, DBL_MAX);
    cauchyStressEigenvals = std::vector<Eigen::Vector2d>(settings.num_triangles);
    cauchyStressEigenvecs = std::vector<Eigen::Matrix<double, 3, 2>>(settings.num_triangles);

    angleDeficits = std::vector<double>(settings.num_nodes, DBL_MAX);
    interiorNodeAngleDeficits = std::vector<double>(settings.num_nodes, DBL_MAX);
    boundaryNodeAngleDeficits = std::vector<double>(settings.num_nodes, DBL_MAX);
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

    if (settings.is_perturbation_of_initial_positions_enabled) {
        perturbInitialPositionsWithRandomNoise(nodes, settings);
    }
}


void Simulation::run_ansatz(int counter) {
    for (int i = 0; i < settings.num_nodes; ++i) {
        nodes[i].pos = nodeAnsatzPositions[i];
    }

    timeSinceLastEquilCheck = 0;
    simulation_status = Dialling;
    DialInFactorCounter = 0;

    if (!settings.is_dialing_from_ansatz_enabled) {

        currDialInFactor = dialInFactorToStartFrom;
        for (std::size_t i = 0; i < DialInFactorValuesToHoldAt.size(); ++i) {
            if (DialInFactorValuesToHoldAt[i] < currDialInFactor) {
                DialInFactorCounter = i;
            }
        }

        timeSinceCurrDiallingInPhaseStarted =
                ((currDialInFactor - DialInFactorValuesToHoldAt[DialInFactorCounter]) /
                 (DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                  - DialInFactorValuesToHoldAt[DialInFactorCounter])) * settings.dial_in_step_time;

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
        if (settings.for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff) {
            for (int i = 0; i < settings.num_triangles; ++i) {
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
        calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, WaitingForEquilibrium, currDialInFactor,
                                                      counter, programmed_metric_infos,
                                                      inverted_programmed_metrics, programmed_taus,
                                                      programmed_second_fundamental_forms, settings);
        updateSecondFundamentalForms(triangles, settings);


        // Temp LU decomp of triangle metric, used to check invertibility.
        Eigen::FullPivLU<Eigen::Matrix<double, 2, 2>> tempMetricDecomp;


        /* Alter inverted_programmed_metrics[progTensorSequenceCounterToStartFrom] and similar to change where
        the programmed quantities are dialling from.*/

        for (int i = 0; i < settings.num_triangles; ++i) {

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
                if (settings.is_dialing_disabled) {
                    inverted_programmed_metrics[progTensorSequenceCounterToStartFrom + 1][i] = (
                            triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                }
            }


            // Other quantities require no inverse, so are easier.
            programmed_taus[progTensorSequenceCounterToStartFrom][i] = programmed_taus[
                    progTensorSequenceCounterToStartFrom + 1][i];
            programmed_second_fundamental_forms[progTensorSequenceCounterToStartFrom][i] = triangles[i].secFF;
            if (settings.is_dialing_disabled) {
                programmed_taus[progTensorSequenceCounterToStartFrom + 1][i] = programmed_taus[
                        progTensorSequenceCounterToStartFrom + 1][i];
                programmed_second_fundamental_forms[progTensorSequenceCounterToStartFrom +
                                                    1][i] = triangles[i].secFF;
            }
        }
    }

    // If set in settings file, clamp all nodes within a certain vertical distance
    // above the node with the lowest initial z value (only works with ansazt clearly).
    if (settings.thicknesses_above_lowest_node_to_clamp_up_to > 0) {
        double minNodeZCoord = nodes[0].pos(2);
        for (int n = 0; n < settings.num_nodes; ++n) {
            if (minNodeZCoord > nodes[n].pos(2)) {
                minNodeZCoord = nodes[n].pos(2);
            }
        }

        for (int n = 0; n < settings.num_nodes; ++n) {
            if (nodes[n].pos(2) - minNodeZCoord <
                settings.thicknesses_above_lowest_node_to_clamp_up_to * settings.sheet_thickness) {
                nodes[n].isClamped = true;
            }
        }
    }
}

void Simulation::setup_imposed_seide_deformations(double &s1, int highest_node, int lowest_node,
                                                  std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    settings.lambda = 0.9;
    settings.cone_angle = asin(pow(settings.lambda, 1.5));
    log_stream.open();
    log_stream << "IMPOSING SEIDE DEFORMATIONS." << std::endl;
    log_stream.close();

    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
    // positions we will force to be those of Seide's setup (found
    // by findng the membrane stress solution for his base state, and
    // integrating the strains etc).
    double intermLengthScaleUpper = 3.0 * sqrt(settings.sheet_thickness * 0.18);
    double intermLengthScaleLower = 3.0 * sqrt(settings.sheet_thickness * 1.8);
    log_stream.open();
    log_stream
            << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
    log_stream.close();
    for (int n = 0; n < settings.num_nodes; ++n) {
        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
            nodes[n].isSeideDisplacementEnabled = true;
        }
    }

    // STORE PERFECT CONE ANSATZ
    for (int n = 0; n < settings.num_nodes; ++n) {
        nodeUnstressedConePosits[n] = nodes[n].pos;
    }
    s1 = sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
              nodes[highest_node].pos(1) * nodes[highest_node].pos(1)) / sin(settings.cone_angle);

//    Eigen::Vector3d testTriCurrCentroid;
//    testTriCurrCentroid = (nodes[triangles[settings.testTriangle].vertexLabels(0)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(1)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(2)].pos) / 3;
    //sTest = sqrt(testTriCurrCentroid(0)*testTriCurrCentroid(0) + testTriCurrCentroid(1)*testTriCurrCentroid(1)) / sin(settings.ConeAngle);

}

void Simulation::setup_glass_cones(int highest_node, int lowest_node) {
    settings.cone_angle = 1.02327019;
    settings.init_slide_z_coord_upper += -tan(settings.cone_angle) *
                                         sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
                                              nodes[highest_node].pos(1) * nodes[highest_node].pos(1));
    settings.init_slide_z_coord_lower += -tan(settings.cone_angle) *
                                         sqrt(nodes[lowest_node].pos(0) * nodes[lowest_node].pos(0) +
                                              nodes[lowest_node].pos(1) * nodes[lowest_node].pos(1));
    log_stream.open();
    log_stream << "USING TWO GLASS CONES FOR SQUASHING." << std::endl;
    log_stream.close();

    // Hijack nodes[i].isOnBoundary to instead label nodes whose
    // forces we will modify to kill any components not tangential to
    // a perfect cone base state. We do this to nodes within an intermediate
    // distance vertically from the ends of the cone.
    double intermLengthScaleUpper = 2.0 * sqrt(settings.sheet_thickness * 0.18);
    double intermLengthScaleLower = 2.0 * sqrt(settings.sheet_thickness * 1.8);
    log_stream.open();
    log_stream << "Apply normal-force-killer within the following vertical distances of the top and bottom: "
               << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
    log_stream.close();
    for (int n = 0; n < settings.num_nodes; ++n) {
        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
            nodes[n].isOnBoundary = true;
        }
    }
}

void Simulation::update_slide_properties() {
    if (!settings.is_controlled_force_enabled) {
        settings.upper_slide_displacement =
                time * settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time;
    } else { // settings.isControlledForceEnabled == true instead
        //settings.upperSlideWeight = (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) * (time * settings.slideSpeedPrefactor / bending_long_time);
        settings.slide_damping_param =
                0.4 * settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness /
                (settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time);
        if (fabs(settings.upper_tot_slide_force + settings.upper_slide_weight) /
            (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) <
            settings.total_slide_force_to_mu_t_sq_ratio_equil_threshold
            && settings.const_slide_weight_fac < 0) {
            if (timeSinceLastEquilCheck > settings.time_between_equil_checks) {
                if (equilibriumCheck(nodes, triangles, settings, log_stream) == EquilibriumReached) {
                    settings.is_slide_just_equilibrated = 1;
                    settings.upper_slide_weight +=
                            settings.slide_weight_dial_speed_fac *
                            (settings.time_step / settings.bending_long_time) *
                            (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness);
                }
                timeSinceLastEquilCheck = 0.0;
            }
        }

        // To do constant-weight slide instead
        if (settings.const_slide_weight_fac > 0) {
            settings.upper_slide_weight =
                    settings.const_slide_weight_fac * settings.shear_modulus * settings.sheet_thickness *
                    settings.sheet_thickness;
        }
    }
    settings.curr_slide_z_coord_upper = settings.init_slide_z_coord_upper - settings.upper_slide_displacement;
}


void Simulation::impose_seide_deformation(double s1, const std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    if (settings.is_seide_deformations_enabled) {
        double pInit = 3.0 * (settings.shear_modulus * settings.sheet_thickness *
                              settings.sheet_thickness); // So we don't have to start all the way from p=0, chosen based on previous sims.
        settings.p = pInit + time * settings.p_speed_prefactor * settings.shear_modulus * settings.sheet_thickness *
                             settings.sheet_thickness / settings.bending_long_time;

        for (int n = 0; n < settings.num_nodes; ++n) {
            if (nodes[n].isSeideDisplacementEnabled || step_count == 0) {

                double tCone = settings.sheet_thickness / sqrt(settings.lambda);

                double polarAng = atan2(nodeUnstressedConePosits[n](1), nodeUnstressedConePosits[n](0));
                double s = sqrt(nodeUnstressedConePosits[n](0) * nodeUnstressedConePosits[n](0) +
                                nodeUnstressedConePosits[n](1) * nodeUnstressedConePosits[n](1)) /
                           sin(settings.cone_angle);
                Eigen::Vector3d uHat;
                uHat << sin(settings.cone_angle) * cos(polarAng), sin(settings.cone_angle) * sin(polarAng), -cos(
                        settings.cone_angle);
                Eigen::Vector3d wHat;
                wHat << -cos(settings.cone_angle) * cos(polarAng), -cos(settings.cone_angle) *
                                                                   sin(polarAng), -sin(settings.cone_angle);
                double u = -settings.p * log(s / s1) /
                           (2.0 * M_PI * settings.youngs_modulus * tCone * sin(settings.cone_angle) *
                            cos(settings.cone_angle));
                double w = -settings.p * (log(s / s1) + settings.poisson_ratio) /
                           (2.0 * M_PI * settings.youngs_modulus * tCone * cos(settings.cone_angle) *
                            cos(settings.cone_angle));
                nodes[n].pos = nodeUnstressedConePosits[n] + u * uHat + w * wHat;
            }
        }
    }
}

bool is_node_lower(const Node &first, const Node &second) {
    return first.pos(2) < second.pos(2);
}

void Simulation::first_step_configuration(double &seide_quotient,
                                          std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    int lowest_node_index = std::min_element(nodes.begin(), nodes.end(), is_node_lower) - nodes.begin();
    int highest_node_index = std::max_element(nodes.begin(), nodes.end(), is_node_lower) - nodes.begin();
    settings.init_slide_z_coord_lower = nodes[0].pos(2);
    settings.init_slide_z_coord_upper = nodes[0].pos(2);


    if (settings.specify_init_slide_z_coord_upper > -100) {
        settings.init_slide_z_coord_upper = settings.specify_init_slide_z_coord_upper;
    } else {
        settings.init_slide_z_coord_upper = nodes[highest_node_index].pos(2);
    }

    if (settings.specify_init_slide_z_coord_lower > -100) {
        settings.init_slide_z_coord_lower = settings.specify_init_slide_z_coord_lower;
    } else {
        settings.init_slide_z_coord_lower = nodes[lowest_node_index].pos(2);
    }

    settings.upper_slide_displacement = 0.0;
    settings.upper_slide_vel = 0.0;
    settings.is_slide_just_equilibrated = 0;
    if (settings.is_controlled_force_enabled) {
        settings.upper_slide_weight = settings.initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq *
                                      (settings.shear_modulus * settings.sheet_thickness *
                                       settings.sheet_thickness);
    }

    // Uncomment the following line for single ridge experiment with point(ish) load to applied tip.
    //settings.initSlideZCoord_upper = settings.initSlideZCoord_lower + 2.0 * settings.SheetThickness; // to apply point(ish) load to tip

    // Readjust for the case of glass cones instead of glass slides. The
    // slide Z coords correspond to the tips of the glass cones.
    if (settings.glass_cones) {
        setup_glass_cones(highest_node_index, lowest_node_index);
    }

    // Instead set up imposed-Seide-deformations idea.
    if (settings.is_seide_deformations_enabled) {
        setup_imposed_seide_deformations(seide_quotient, highest_node_index, lowest_node_index,
                                         nodeUnstressedConePosits);
    }
}

void Simulation::begin_equilibrium_search(int counter) {
    currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter + 1];
    calcTriangleGeometries_and_DialledProgTensors(nodes, triangles, simulation_status, currDialInFactor,
                                                  counter, programmed_metric_infos,
                                                  inverted_programmed_metrics, programmed_taus,
                                                  programmed_second_fundamental_forms, settings);
    simulation_status = WaitingForEquilibrium;
    // NB an EquilCheck has not actually just occurred, but this has the
    //desired effect of ensuring that each DialInFactor value is held
    //for at least one TimeBetweenEquilChecks.
    timeSinceLastEquilCheck = 0.0;

    settings.num_damp_factor =
            settings.equilibriation_damping *
            settings.damping_scale; // Set damping factor to waiting phase value.

    log_stream.open();
    log_stream << "\tReached dial-in factor of " << DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
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
    settings.upper_tot_slide_force = upperAndLowerTotSlideForces.first;
}

void Simulation::update_dial_in_factor() {
    currDialInFactor = DialInFactorValuesToHoldAt[DialInFactorCounter] +
                       (DialInFactorValuesToHoldAt[DialInFactorCounter + 1]
                        - DialInFactorValuesToHoldAt[DialInFactorCounter]) *
                       timeSinceCurrDiallingInPhaseStarted / settings.dial_in_step_time;
}


void Simulation::save_and_print_details(int counter, double duration_us,
                                        std::pair<double, double> upperAndLowerTotSlideForces) {
    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings);
    if (settings.is_energy_densities_printed) {
        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
    }

    // Here we hack the angleDeficits to instead tell us whether a node has a Seide displacement imposed or not.
    for (int n = 0; n < settings.num_nodes; ++n) {
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
    log_stream << log_prefix()
               << "Last step's execution time " << duration_us << " us."
               << std::endl;

    log_stream.close();

    std::ofstream forceDistFile;
    forceDistFile.open(outputDirName + "/force_displacement_vals.txt", std::ofstream::app);
    forceDistFile << std::scientific << std::setprecision(14)
                  << upperAndLowerTotSlideForces.first /
                     (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) << " "
                  << upperAndLowerTotSlideForces.second /
                     (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) << " "
                  << settings.upper_slide_displacement / settings.sample_char_length << " "
                  << settings.curr_slide_z_coord_upper << " "
                  << settings.init_slide_z_coord_lower << " "
                  << settings.is_slide_just_equilibrated << " "
                  << settings.upper_slide_weight /
                     (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness)
                  << std::endl;
    forceDistFile.close();
}

std::string Simulation::log_prefix() {
    std::stringstream ss;
    ss << std::setprecision(3) << std::fixed << getRealTime() << ", step = " << step_count << ", time = " << time
       << ", dial-in factor = " <<
       currDialInFactor << ": ";
    return ss.str();
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
    if (settings.is_energy_densities_printed) {
        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings);
    }
    writeVTKDataOutput(nodes, triangles, step_count, time, currDialInFactor, counter,
                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                       cauchyStressEigenvals, cauchyStressEigenvecs, settings, outputDirName);
    throw std::runtime_error("Unexpectedly large force at step " + std::to_string(step_count) + ".");
}

void Simulation::check_for_equilibrium() {
    log_stream.open();
    log_stream << log_prefix() << "Checking for equilibrium. "
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
    simulation_status = Dialling;
    timeSinceCurrDiallingInPhaseStarted = 0.0;
    DialInFactorCounter += 1;
    settings.num_damp_factor = settings.dial_in_damping *
                               settings.damping_scale; // Set damping factor back to dialling phase value.
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
        simulation_status = Dialling;
        currDialInFactor = 0;
        DialInFactorCounter = 0;
        timeSinceCurrDiallingInPhaseStarted = 0;
    }


    /* Handle cases where settings.DialInStepTime is zero or very small. This is
    a slightly hacky way to ensure that no dialling actually occurs, and the
    simulation jumps to a status = waitingForEquilibrium state. The magic number is
    just chosen to be recognisable for debugging.*/
    if (settings.dial_in_step_time < settings.time_step) {
        settings.dial_in_step_time = 0.0;
        timeSinceCurrDiallingInPhaseStarted = 1.23456789;
    }


    // Begin dynamical evolution of node positions and velocities.
    log_stream.open();
    log_stream << "\nBeginning dynamical evolution.\n" << std::endl;
    log_stream << "\nCREATING VECTOR TO STORE UNSTRESSED CONE NODE POSITIONS.\n" << std::endl;
    log_stream.close();
    std::vector<Eigen::Vector3d> nodeUnstressedConePosits(settings.num_nodes);
    double seide_quotient = DBL_MAX;

    std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes = getCorrespondingTrianglesForNodes(
            triangles, nodes);

    while (DialInFactorCounter <= DialInFactorValuesToHoldAt.size() - 2) {
        if (step_count == 0) { first_step_configuration(seide_quotient, nodeUnstressedConePosits); }

        update_slide_properties();

        if (settings.is_seide_deformations_enabled) {
            impose_seide_deformation(seide_quotient, nodeUnstressedConePosits);
        }

        zeroForces(nodes);

        // Check if still in dialling in phase or whether it is time to wait for
        // equilibrium.
        if (timeSinceCurrDiallingInPhaseStarted >= settings.dial_in_step_time && simulation_status == Dialling) {
            begin_equilibrium_search(counter);
        }

        /* If not waiting for equilibrium, set the current value of the Dial-In
        Factor, based on linear dialling in between the previously calculated
        'checkpoint' values. */
        if (simulation_status == Dialling) { update_dial_in_factor(); }

        std::pair<double, double> upperAndLowerTotSlideForces;
        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        progress_single_step(counter, upperAndLowerTotSlideForces, correspondingTrianglesForNodes);
        auto duration = std::chrono::steady_clock::now() - begin;
        double duration_us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();


        /* Write output data regularly.
        Can be switched off with settings.PrintFrequency < 0.0.
        Doing the write-out at this point in the loop means the node positions
        and the triangle geometry data match in the output, which is desirable!*/
        if ((step_count % settings.inverse_print_rate == 0 && settings.inverse_print_rate > 0) ||
            settings.is_slide_just_equilibrated == 1) {
            save_and_print_details(counter, duration_us, upperAndLowerTotSlideForces);
            settings.is_slide_just_equilibrated = 0;
        }


        if (!settings.is_controlled_force_enabled) {
            /* If last check for equilibrium or last reaching of a new DialInFactor
             value was more than TimeBetweenEquilChecks ago, check for equilibrium.
             Also print total stretching and bending energies.*/
            if (timeSinceLastEquilCheck > settings.time_between_equil_checks &&
                simulation_status == WaitingForEquilibrium) {
                check_for_equilibrium();
            }

            /* If equilibrium reached, write output data to file, and move to next
            'dialling in' phase. */
            if (simulation_status == EquilibriumReached) {
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
        time += settings.time_step;
        timeSinceLastEquilCheck += settings.time_step;
        timeSinceCurrDiallingInPhaseStarted += settings.time_step;
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
