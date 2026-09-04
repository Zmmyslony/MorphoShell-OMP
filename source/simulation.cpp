//
// Created by Michał Zmyślony on 12/02/2024.
//

#define _USE_MATH_DEFINES

#include "simulation.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <numeric>

#include "functions/getRealTime.hpp"
#include "functions/extract_Just_Filename.hpp"
#include "initialDirAndFileHandling.hpp"
#include "readVtk.hpp"
#include "calculations/calcTrianglesIncidentOnNodes.hpp"
#include "calculations/calcTriangleAdjacencies_And_Edges.hpp"
#include "functions/kahanSum.hpp"
#include "calculations/calcNodeNeighbours.hpp"
#include "calculations/calc_nonVertexPatchNodes_and_MatForPatchDerivs.hpp"
#include "setRemainingInitCond_and_NodeMasses.hpp"
#include "functions/perturbInitialPositionsWithRandomNoise.hpp"
#include "calculations/calcEnergiesAndStresses.hpp"
#include "calculations/calcCurvatures.hpp"
#include "exportVtk.hpp"
#include "calculations/calcDeformationForces.hpp"
#include "physics/cone.h"

// #ifndef TIMING
// #define TIMING
// #endif

double getMinimumTriangleHeight(const std::vector<Triangle>& triangles) {
    double min_value = DBL_MAX;
#pragma omp parallel for reduction(min: min_value)
    for (int i = 0; i < triangles.size(); i++) {
        min_value = triangles[i].getHeight() ? (triangles[i].getHeight() < min_value) : min_value;
    }

    return min_value;
}

// Gets the height of the highest triangle

double getMaximumTriangleHeight(const std::vector<Triangle>& triangles) {
    double max_value = DBL_MIN;
#pragma omp parallel for reduction(max: max_value)
    for (int i = 0; i < triangles.size(); i++) {
        max_value = triangles[i].getHeight() ? (triangles[i].getHeight() < max_value) : max_value;
    }
    return max_value;
}
 
#if _MSC_VER
#else
std::pair<double, int> pair_min(const std::pair<double, int> &lhs, const std::pair<double, int> &rhs) {
    return lhs.first < rhs.first ? lhs : rhs;
}

std::pair<double, int> pair_max(const std::pair<double, int> &lhs, const std::pair<double, int> &rhs) {
    return lhs.first > rhs.first ? lhs : rhs;
}

#pragma omp declare reduction(pair_max: std::pair<double, int>: omp_out=pair_max(omp_out, omp_in))\
initializer(omp_priv={DBL_MIN, 0})
#endif



bool isForceThresholdExceeded(const Node &node, double force_limit) {
    return node.force.norm() >= force_limit;
}


void logForceThresholdExceeded(Node &node, std::vector<Triangle> &triangles, const Settings &settings, int step) {
    std::stringstream msg;
    msg << " ----------------------------------------" << std::endl;
    msg << " ------------CRASH REPORT----------------" << std::endl;
    msg << " ----------------------------------------" << std::endl;
    msg << "Error occured at time step " << step << std::endl;
    msg << "First offending node and its incident triangles: " << std::endl;

    msg << node.display().str();
    for (int t = 0; t < node.incidentTriLabels.size(); ++t) {
        msg << triangles[node.incidentTriLabels[t]].display().str();
    }
    std::cout << msg.str();
    throw std::runtime_error(msg.str() + "Suspiciously high force at node " + std::to_string(node.label) + std::string(" (") +
                             std::to_string(node.force.norm() / 1e5 * settings.getForceScale()) + " of limit).");
}


typedef teestream<char, std::char_traits<char>> basic_teestream;

std::pair<double, double> mean_dev(const std::vector<int>& times) {
    double mean = static_cast<double>(std::accumulate(times.begin(), times.end(), 0.0)) / static_cast<double>(times.size());
    double dev = 0;
    for (const auto& el : times) { dev += std::pow(el - mean, 2); }
    dev /= double(times.size());
    return {mean, std::sqrt(dev)};
}

//   Redirection of std::cout to both console and file following
//   https://stackoverflow.com/a/41834129/9204102
void Simulation::setupLogging() {
    log_filename = output_dir_name + "/log.txt";
    out = std::ofstream(log_filename);
    std::streambuf* out_buf = out.rdbuf(); //Get streambuf for output stream

    std::streambuf* tee_buf = new basic_teestream(out_buf, cout_buf); //create new teebuf

    std::cout.rdbuf(tee_buf); //Redirect cout
    std::cout << std::defaultfloat << std::setprecision(6);
}


void Simulation::setup_filenames(int argc, char* argv[]) {
    init_string += "Simulation run began at: " + getRealTime() + "\n";
    init_string += "The command that was run was:\n";
    for (int i = 0; i < argc; i++) {
        init_string += argv[i];
        init_string += " ";
    }
    init_string += "\n\n";

    if (argc < 2) {
        throw std::runtime_error("Insufficient number of input files. At least, core_config and "
            "simulation .vtk are required.");
    }
    std::vector<fs::path> vtk_paths;
    std::vector<fs::path> config_paths;
    for (int i = 1; i < argc; i++) {
        fs::path current_path = argv[i];
        if (current_path.extension() == ".vtk") { vtk_paths.emplace_back(current_path); } else if (current_path.
            extension() == ".cfg") { config_paths.emplace_back(current_path); } else {
            std::string message = "Unrecognised extension for input file " + std::string(argv[i]) + ".";
            throw std::runtime_error(message);
        }
    }

    if (vtk_paths.empty()) { throw std::runtime_error("Missing simulation .vtk file. "); }
    if (vtk_paths.size() > 2) {
        throw std::runtime_error("Too many input files with .vtk format. Only two are allowed -"
            " simulation and ansatz");
    }
    initialisation_filename = vtk_paths.front().string();
    if (vtk_paths.size() == 2) { ansatz_filename = vtk_paths.back().string(); }
    if (config_paths.empty()) { throw std::runtime_error("Missing configuration .cfg file."); }

    output_dir_name = directory_setup(init_string, config_paths, vtk_paths).string();

    setupLogging();
    std::cout << init_string << std::endl;
}


void Simulation::read_settings_new(int argc, char* argv[]) {
    std::vector<fs::path> paths(argc);
    for (int i = 0; i < argc; i++) { paths[i] = argv[i]; }
    settings = Settings(paths);
}


void Simulation::read_vtk_data(const CoreConfig& config) {
    readVTKData(nodes, triangles, programmed_metric_infos, inverted_programmed_metrics, programmed_taus,
                programmed_second_fundamental_forms, settings.getCore().isLceModeEnabled(),
                initialisation_filename, initial_stage, dialInFactorToStartFrom, nodeAnsatzPositions, ansatz_filename,
                config);

    stage_count = static_cast<int>(inverted_programmed_metrics.size());
    updateProgrammedValues(1);
}

void Simulation::updateProgrammedValues(int current_stage) {
    int i = current_stage;
#pragma omp parallel for
    for (int j = 0; j < triangles.size(); j++) {
        triangles[j].next_programmed_metric_info = programmed_metric_infos[i][j];
        triangles[j].next_programmed_metric_inv = inverted_programmed_metrics[i][j];
        triangles[j].next_programmed_tau = programmed_taus[i][j];
        triangles[j].next_programmed_second_fundamental_form =programmed_second_fundamental_forms[i][j];
    }
}


void Simulation::configure_nodes() const {
    std::cout << "Nodes count = " << nodes.size() << std::endl;
    std::cout << "Triangle count = " << triangles.size() << std::endl;

    double average_node_index_distance = 0;
#pragma omp parallel for reduction(+: average_node_index_distance)
    for (int i = 0; i < triangles.size(); i++) {
        const Triangle &triangle = triangles[i];
        double max_distance = std::max(
            std::max(abs(triangle.vertexLabels[2] - triangle.vertexLabels[0]),
                     abs(triangle.vertexLabels[1] - triangle.vertexLabels[0])),
            abs(triangle.vertexLabels[2] - triangle.vertexLabels[1])
        );
        average_node_index_distance += max_distance / triangles.size();
    }

    std::cout << "Average node distance: " << int(average_node_index_distance) << " (2 is ideal)." << std::endl;

    if (triangles.size() < 50) {
        std::cerr << "Your mesh has a small number of triangles. \nBeware that the code "
            "is likely to be less accurate in this case, \nand unforeseen bugs are more "
            "likely in extreme cases." << std::endl;
    }
}


void Simulation::configure_topological_properties() {
    calcTrianglesIncidentOnNodes(nodes, triangles);
    calcTriangleAdjacencies_And_Edges(nodes, triangles, edges);

    /* Calculate and print the number of boundary and non-boundary edges.
    Also label the nodes connected by boundary edges as boundary nodes.
    Also calculate the total initial perimeter of the sample, to be used as a
    characteristic sample length in estimating characteristic times, time steps
    etc.*/
    int boundary_edge_count = 0;
    std::vector<double> initBoundaryEdgeLengths(edges.size());

#pragma omp parallel for reduction(+: boundary_edge_count)
    for (int i = 0; i < edges.size(); i++) {
        if (edges[i].isBoundary()) {
            boundary_edge_count += 1;

            nodes[edges[i].nodeLabels.first].isOnBoundary = true;
            nodes[edges[i].nodeLabels.second].isOnBoundary = true;

            //If chosen in settings, clamp whole boundary in addition to clamp
            //indicators from data file.
            if (settings.getCore().isBoundaryClamped()) {
                nodes[edges[i].nodeLabels.first].clamp(settings.getCore());
                nodes[edges[i].nodeLabels.second].clamp(settings.getCore());
            }

            // Store initial length of this boundary edge
            initBoundaryEdgeLengths[i] = (nodes[edges[i].nodeLabels.first].position - nodes[edges[i].nodeLabels.second].position).norm();
        }
    }

    double initPerimeter = kahanSum(initBoundaryEdgeLengths);
    characteristic_long_length = initPerimeter;

    if (3 * triangles.size() != 2 * edges.size() - boundary_edge_count) {
        throw std::runtime_error(
            "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
    }

    std::cout << "Edges count = " << edges.size() << std::endl;
    std::cout << "Boundary edges count = " << boundary_edge_count << std::endl;

    std::cout << std::endl;
    std::cout << "Initial perimeter = " << initPerimeter << std::endl;
}


void Simulation::configure_triangles() {
    int numBoundaryTriangles = 0;
#pragma omp parallel for reduction(+ : numBoundaryTriangles)
    for (int i = 0; i < triangles.size(); i++) { if (triangles[i].isOnBoundary) { numBoundaryTriangles += 1; } }

    std::cout << std::endl;
    std::cout << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
    int hole_count = 1 + edges.size() - nodes.size() - triangles.size();
    std::cout << "Number of holes in mesh = " << hole_count
        << std::endl; // From Euler's formula for a planar graph.

    if (hole_count < 0) {
        throw std::runtime_error(
            "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
    }
}

void Simulation::orient_node_labels() {
    updateTriangleProperties(-10);
    Eigen::Vector3d tempZAxisVec;
    tempZAxisVec << 0.0, 0.0, 1.0;
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        auto sides = triangles[i].getCurrentSides();
        if (tempZAxisVec.dot(sides.col(0).cross(sides.col(1))) < 0) {
            int tempLabel = triangles[i].vertexLabels(2);
            triangles[i].vertexLabels(2) = triangles[i].vertexLabels(1);
            triangles[i].vertexLabels(1) = tempLabel;
        }
    }
    updateTriangleProperties(-10);
}


void Simulation::set_initial_conditions() {
    setRemainingInitCond_and_NodeMasses(nodes, triangles, edges,
                                        settings);
}

void Simulation::find_smallest_element() {
    characteristic_short_length = std::min_element(triangles.begin(), triangles.end(),
                                                   [](Triangle& first, Triangle& second) {
                                                       return first.getLinearSize() < second.getLinearSize();
                                                   })->getLinearSize();
    auto largest_tau = DBL_MIN;
#pragma omp parallel for reduction(max: largest_tau)
    for (int i = 0; i < programmed_taus.size(); i++) {
        for (auto& tau : programmed_taus[i]) { if (tau < largest_tau) { largest_tau = tau; } }
    }
    characteristic_length_over_tau = characteristic_short_length / sqrt(largest_tau);

    std::cout << std::endl;
    std::cout << "Sheet thickness = " << settings.getCore().getThickness() << std::endl;
    std::cout << "Approx smallest element linear size = " << characteristic_short_length << std::endl;
}


void Simulation::setup_characteristic_scales() {
    std::vector<double> initAreas(triangles.size());
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) { initAreas[i] = triangles[i].initArea; }
    double initial_area = kahanSum(initAreas);
    settings.SetupCharacteristicSizes(initial_area, characteristic_short_length);
}


void Simulation::setup_equilibrium_dial_in_factors() {
    double tempDialInFactor = 0;
    if (settings.getCore().getDialInResolution() > 0 && settings.getDurationPhase() >= 0) {
        while (tempDialInFactor < 1) {
            dial_in_phases.push_back(tempDialInFactor);
            tempDialInFactor += settings.getCore().getDialInResolution();
        }
        dial_in_phases.push_back(1);
    } else {
        std::cout << "settings.DialInResolution and settings.DialInStepTime must "
            "both be >0 and >=0 respectively, and they aren't currently. Aborting."
            << std::endl;
        throw std::runtime_error("settings.DialInResolution and settings.DialInStepTime must "
            "both be >0 and >=0 respectively, and they aren't currently. Aborting.");
    }
}

void Simulation::initialise_simulation_vectors() {
    gaussCurvatures = std::vector<double>(triangles.size(), DBL_MAX);
    meanCurvatures = std::vector<double>(triangles.size(), DBL_MAX);
    stretchEnergies = std::vector<double>(triangles.size(), DBL_MAX);
    bendEnergies = std::vector<double>(triangles.size(), DBL_MAX);
    stretchEnergyDensities = std::vector<double>(triangles.size(), DBL_MAX);
    bendEnergyDensities = std::vector<double>(triangles.size(), DBL_MAX);
    kineticEnergies = std::vector<double>(nodes.size(), DBL_MAX);
    strainMeasures = std::vector<double>(triangles.size(), DBL_MAX);
    cauchyStressEigenvals = std::vector<Eigen::Vector2d>(triangles.size());
    cauchyStressEigenvecs = std::vector<Eigen::Matrix<double, 3, 2>>(triangles.size());

    angleDeficits = std::vector<double>(nodes.size(), DBL_MAX);
    interiorNodeAngleDeficits = std::vector<double>(nodes.size(), DBL_MAX);
    boundaryNodeAngleDeficits = std::vector<double>(nodes.size(), DBL_MAX);
}


void Simulation::init(int argc, char* argv[], int threads) {
    setup_filenames(argc, argv);
    read_settings_new(argc, argv);
    settings.setThreads(threads);
    read_vtk_data(settings.getCore());
    configure_nodes();
    configure_topological_properties();
    configure_triangles();
    configureNodeAdjacency(nodes, edges);
    createNodePatches(nodes, triangles, settings.getCore().getPatchMatrixThreshold());

    node_force_proxy = std::vector<Eigen::Vector3d>(triangles.size() * 6);
    assignForceLocationsToNodes(triangles, nodes, node_force_proxy);
    orient_node_labels();

    set_initial_conditions();
    find_smallest_element();

    settings.SetupDialInTime(characteristic_long_length);
    settings.SetupStepTime(characteristic_short_length);
    settings.SetupPrintFrequency();

    setup_characteristic_scales();
    setup_equilibrium_dial_in_factors();

    initialise_simulation_vectors();

    if (settings.getCore().isInitialPositionsPerturbed()) {
        perturbInitialPositionsWithRandomNoise(nodes, characteristic_short_length);
    }
}


void Simulation::run_ansatz(int counter) {
#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i].position = nodeAnsatzPositions[i];
        nodes[i].prev_position = nodeAnsatzPositions[i];
    }

    time_equilibriation = 0;
    simulation_status = Dialling;
    phase_counter = 0;

    if (settings.getCore().isAnsatzMetricUsed()) {
        std::cout << "\nINFO Assuming ansatz state to be relaxed (initial programmed tensors calculated from geometry)."
            << std::endl;

        dial_in_factor = 0.0;
        time_phase = 0.0;

        // Calculate all necessary geometry for the ansatz state.
        updateTriangleProperties(counter);
        // updateFirstFundamentalForms(triangles, settings.getCore());
        // updateSecondFundamentalForms(triangles, settings.getCore());


        /* Alter inverted_programmed_metrics[initial_stage] and similar to change where
        the programmed quantities are dialling from.*/
#pragma omp parallel for
        for (int i = 0; i < triangles.size(); i++) {
            // Test for invertibility of metric before taking inverse.
            Eigen::Matrix<double, 3, 2> deformationGradient = triangles[i].getDeformationGradient();
            Eigen::Matrix2d metric = deformationGradient.transpose() * deformationGradient;
            Eigen::FullPivLU<Eigen::Matrix<double, 2, 2>> tempMetricDecomp(metric);
            if (!tempMetricDecomp.isInvertible()) {
                throw std::runtime_error(triangles[i].display().str() +
                    "At least one triangle [" + std::to_string(i) +
                    "] had a non-invertible metric in the ansatz state. \n"
                    "This should not occur in a reasonable mesh. Aborting.");
            }
            size_t stage = initial_stage;
            if (settings.getCore().isFirstTensorSkipped()) { stage++; }
            triangles[i].previous_metric_inverse = metric.inverse();
            // triangles[i].programmed_tau = programmed_taus[i][initial_stage];
            triangles[i].previous_second_fundamental_form = triangles[i].getSecondFundamentalForm();
        }
    } else {
        dial_in_factor = dialInFactorToStartFrom;
        for (int i = 0; i < dial_in_phases.size(); i++) {
            if (dial_in_phases[i] < dial_in_factor) { phase_counter = i; }
        }
        double dial_in_since_last_phase = dial_in_factor - dial_in_phases[phase_counter];
        double dial_in_phase_duration = dial_in_phases[phase_counter + 1] - dial_in_phases[phase_counter];

        time_phase = dial_in_since_last_phase / dial_in_phase_duration * settings.getDurationPhase();

        if (settings.getCore().isFirstTensorSkipped()) {
#pragma omp parallel for
            for (int i = 0; i < triangles.size(); i++) {
                triangles[i].programmed_metric_info = programmed_metric_infos[i][initial_stage + 1];
                triangles[i].previous_second_fundamental_form = programmed_second_fundamental_forms[i][initial_stage + 1];
            }
        }
    }
}

// Seide, 1956, Axisymmetrical Buckling of Circular Cones Under Axial Compression,
// https://doi.org/10.1115/1.4011410
void Simulation::setup_imposed_seide_deformations(double& s1, int highest_node, int lowest_node,
                                                  std::vector<Eigen::Vector3d>& nodeUnstressedConePosits) {
    lambda = 0.9;
    double cone_angle = asin(pow(lambda, 1.5));
    std::cout << "IMPOSING SEIDE DEFORMATIONS." << std::endl;

    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
    // positions we will force to be those of Seide's setup (found
    // by findng the membrane stress solution for his base state, and
    // integrating the strains etc).
    double intermLengthScaleUpper = 3.0 * sqrt(settings.getCore().getThickness() * 0.18);
    double intermLengthScaleLower = 3.0 * sqrt(settings.getCore().getThickness() * 1.8);

    std::cout
        << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
        << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;

    for (int n = 0; n < nodes.size(); ++n) {
        if (((nodes[highest_node].position(2) - nodes[n].position(2)) < intermLengthScaleUpper) ||
            ((nodes[n].position(2) - nodes[lowest_node].position(2)) < intermLengthScaleLower)) {
            nodes[n].isSeideDisplacementEnabled = true;
        }
    }

    // STORE PERFECT CONE ANSATZ
    for (int n = 0; n < nodes.size(); ++n) { nodeUnstressedConePosits[n] = nodes[n].position; }
    s1 = sqrt(nodes[highest_node].position(0) * nodes[highest_node].position(0) +
        nodes[highest_node].position(1) * nodes[highest_node].position(1)) / sin(cone_angle);
}


void Simulation::first_step_configuration() {
    slides = settings.getSlides();
    cones = settings.getCones();
    std::vector<Eigen::Vector3d> node_positions(nodes.size());
#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) { node_positions[i] = nodes[i].position; }

    for (auto& slide : slides) { slide.initialise(node_positions, settings.getDurationPhase()); }
    for (auto& cone : cones) { cone.initialise(node_positions, settings.getDurationPhase()); }
}

void Simulation::begin_equilibrium_search(int counter) {
    dial_in_factor = dial_in_phases[phase_counter + 1];
    updateTriangleProperties(counter);
    simulation_status = WaitingForEquilibrium;
    // NB an EquilCheck has not actually just occurred, but this has the
    // desired effect of ensuring that each DialInFactor value is held
    // for at least one TimeBetweenEquilChecks.
    time_equilibriation = settings.getTimeBetweenEquilibriumChecks();

    settings.useEquilibriumDamping();

    std::cout << "\tReached dial-in factor of " << dial_in_phases[phase_counter + 1]
        << ". Waiting for equilibrium." << std::endl;
}


void Simulation::add_interaction_forces() {
    if (slides.size() == 0 &&
        cones.size() == 0 &&
        settings.getMagneticField().size() == 0) {
        return;
    }

    double shared_interaction_force = 0;
    for (auto& slide : slides) {
        shared_interaction_force = 0;
#pragma omp parallel for reduction (+ : shared_interaction_force)
        for (int i = 0; i < nodes.size(); i++) {
            shared_interaction_force += slide.addInteractionForce(nodes[i].position,
                                                                  nodes[i].force,
                                                                  settings.getCore().getShearModulus(),
                                                                  settings.getCore().getThickness(),
                                                                  nodes[i].area);
        }
        slide.setTotalInteractionForce(shared_interaction_force);
    }

    for (auto& cone : cones) {
        shared_interaction_force = 0;
#pragma omp parallel for reduction (+ : shared_interaction_force)
        for (int i = 0; i < nodes.size(); i++) {
            shared_interaction_force += cone.addInteractionForce(nodes[i].position,
                                                                 nodes[i].force,
                                                                 settings.getCore().getShearModulus(),
                                                                 settings.getCore().getThickness());
        }
        cone.setTotalInteractionForce(shared_interaction_force);
    }

    for (auto& field : settings.getMagneticField()) {
        Eigen::Vector3d magnetic_field = field.magnetic_field;
#pragma omp parallel for
        for (int i = 0; i < triangles.size(); i++) { triangles[i].updateMagneticForce(magnetic_field); }
#pragma omp parallel for
        // Applies forces to each node.
        for (int i = 0; i < nodes.size(); i++) { nodes[i].addForce(); }
    }

#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) { nodes[i].apply_boundary_conditions(); }
}


void Simulation::updateTriangleProperties(int counter) {
    const double dial_in_factor_root = sqrt(dial_in_factor);
    const bool is_LCE_metric_used = settings.getCore().isLceModeEnabled();
    const bool is_stimulation_modulated = settings.getCore().isStimulationModulated();
    const double transfer_coefficient = 25 * settings.getTimeStepSize() / settings.getDurationPhase();

    // const double min_height = getMinimumTriangleHeight(triangles);
    // const double max_height = getMaximumTriangleHeight(triangles);
    const double min_height = -1;
    const double max_height = 1;

    const double stretching_pre_factor =
        0.5 * settings.getCore().getThickness() * settings.getCore().getShearModulus();
    const double bending_pre_factor =
        0.5 * pow(settings.getCore().getThickness(), 3) * settings.getCore().getShearModulus() / 12;
    const double j_pre_factor = settings.getCore().getGentFactor() / pow(settings.getCore().getThickness(), 2);
    const double poisson_ratio = settings.getCore().getPoissonRatio();

#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        if (simulation_status == Dialling) {
            triangles[i].updateProgrammedQuantities(counter, dial_in_factor, dial_in_factor_root, is_LCE_metric_used,
                                                    is_stimulation_modulated, transfer_coefficient, min_height,
                                                    max_height);
        }
        triangles[i].updateGeometricProperties(bending_pre_factor, j_pre_factor, poisson_ratio, stretching_pre_factor);
        // triangles[i].updateElasticForce(bending_pre_factor, j_pre_factor, stretching_pre_factor, poisson_ratio);
    }
}

void Simulation::add_node_forces() {
    double shared_interaction_force = 0;
    GravityConfig gravity = settings.getGravity();
    const double damping_multiplier = settings.getDampingMultiplier();
#pragma omp parallel for  reduction (+ : shared_interaction_force) shared(gravity)
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i].updateForce();
        shared_interaction_force += nodes[i].add_damping(damping_multiplier);
        nodes[i].add_gravity(gravity);
        nodes[i].apply_boundary_conditions();
    }
    damping_power_loss = shared_interaction_force;
}

long long int Simulation::progress_single_step(int counter) {
    auto begin = std::chrono::high_resolution_clock::now();
    updateTriangleProperties(counter);
    add_node_forces();
    add_interaction_forces();

    auto duration = std::chrono::high_resolution_clock::now() - begin;
    long long duration_us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    benchmarking_mechanics_duration += duration_us;
    return duration_us;
}

void Simulation::update_dial_in_factor() {
    double starting_dial_in_factor = dial_in_phases[phase_counter];
    double final_dial_in_factor = dial_in_phases[phase_counter + 1];
    double relative_time = time_phase / settings.getDurationPhase();
    dial_in_factor = starting_dial_in_factor + (final_dial_in_factor - starting_dial_in_factor) * relative_time;
}

long long int Simulation::export_vtk(int counter) {
    double stretch_prefactor = 0.5 * settings.getCore().getThickness() * settings.getCore().getShearModulus();
    for (int i = 0; i < triangles.size(); i++) {
        stretchEnergyDensities[i] = triangles[i].getStretchingEnergyDensity(stretch_prefactor);
        stretchEnergies[i] = stretchEnergyDensities[i] * triangles[i].initArea;
    }

    return ::writeVTKDataOutput(nodes, triangles, step_count, time_global, dial_in_factor, counter,
                                gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                                boundaryNodeAngleDeficits,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings,
                                output_dir_name, damping_power_loss);
}

void Simulation::save_and_print_details(int counter, long long int duration_us) {
    auto begin = std::chrono::high_resolution_clock::now();
    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings.getCore());
    if (settings.getCore().isEnergyPrinted()) {
        calcEnergiesAndStresses(nodes, triangles,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings.getCore());
    }
    auto duration = std::chrono::high_resolution_clock::now() - begin;
    long long prep_us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();

    long long export_us = export_vtk(counter);
    benchmarking_export_duration += export_us;

    std::cout << log_prefix()
        << "Last step's execution time " << duration_us << " us. Export time " << prep_us + export_us << " us ("
        << prep_us << "/" << export_us << " us prep/write)" << std::endl;
    mechanics_times.push_back(duration_us);
    pre_export_times.push_back(prep_us);
    export_times.push_back(export_us);
}

std::string Simulation::log_prefix() const {
    std::stringstream ss;
    ss << std::setprecision(3) << std::fixed << getRealTime() << ", step = " << step_count << ", time = " << time_global
        << ", dial-in factor = " <<
        dial_in_factor << ": ";
    return ss.str();
}

void Simulation::error_large_force(int counter) {
    std::cout << "At " << getRealTime() << ", stepCount = " << step_count << ", time = " << time_global
        << ", dial-in factor = " << dial_in_factor
        << " a force was suspiciously large, there is probably a problem. Writing VTK output and then aborting."
        << std::endl;

    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings.getCore());
    if (settings.getCore().isEnergyPrinted()) {
        calcEnergiesAndStresses(nodes, triangles,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings.getCore());
    }
    export_vtk(counter);
    throw std::runtime_error("Unexpectedly large force at step " + std::to_string(step_count) + ".");
}

void Simulation::check_for_equilibrium() {
    std::cout << log_prefix() << "Checking for equilibrium. " << std::endl;
    simulation_status = equilibrium_check();
    time_equilibriation = 0.0;
}

void Simulation::configure_smallest_tau() {
    double current_smallest_tau = DBL_MAX;
#pragma omp parallel for reduction(min: current_smallest_tau)
    for (int i = 0; i < triangles.size(); i++) {
        if (triangles[i].dialledProgTau < current_smallest_tau) {
            current_smallest_tau = triangles[i].next_programmed_tau;
        }
    }
    stage_stretching_wave_speed = sqrt(current_smallest_tau * settings.getCore().getShearModulus() /
                                        settings.getCore().getDensity());
}

SimulationStatus Simulation::equilibrium_check() {
#if _MSC_VER
    double max_velocity = DBL_MIN;
    double max_force = DBL_MIN;
#pragma omp parallel for reduction(max: max_velocity, max_force)
#else
    std::pair<double, int> velocity_pair = {DBL_MIN, -1};
    std::pair<double, int> force_pair = {DBL_MIN, -1};
#pragma omp parallel for reduction(pair_max: velocity_pair, force_pair)
#endif
    for (int i = 0; i < nodes.size(); i++) {
        double non_damping_force;
        if (!settings.getCore().isGradientDescentDynamics()) {
            non_damping_force = (nodes[i].force +
                                      (settings.getDampingFactor() * nodes[i].mass * nodes[i].velocity /
                                       settings.getCore().getDensity())).norm();
        } else {
            non_damping_force = nodes[i].force.norm();
        }
        std::pair<double, int> current_velocity_pair = {nodes[i].velocity.norm(), i};
        std::pair<double, int> current_force_pair = {non_damping_force, i};
#if _MSC_VER
        max_velocity = std::max(max_velocity, current_velocity_pair.first);
		max_force = std::max(max_force, current_force_pair.first);
#else
        velocity_pair = pair_max(velocity_pair, current_velocity_pair);
        force_pair = pair_max(force_pair, current_force_pair);
#endif
    }

#if _MSC_VER
    double max_non_damp_force = -1;
    int max_non_damp_force_node = max_force;

    int max_relative_speed_node = -1;
    double max_relative_speed = max_velocity / stage_stretching_wave_speed;
#else
    double max_non_damp_force = force_pair.first;
    int max_non_damp_force_node = force_pair.second;

    int max_relative_speed_node = velocity_pair.second;
    double max_relative_speed = velocity_pair.first / stage_stretching_wave_speed;
#endif
    SimulationStatus status;
    if (max_non_damp_force < settings.getForceScale()
        && max_relative_speed < settings.getCore().getEquilibriumSpeedScale()) {

        std::cout << "\tEquilibrium reached." << std::endl;
        status = EquilibriumReached;
    } else {
        std::cout << "\tEquilibrium not reached." << std::endl;
        status = WaitingForEquilibrium;
    }
    std::cout << "\tRatio of max non-damping force to characteristic force = "
              << max_non_damp_force / settings.getForceScale()
              << " (node " << max_non_damp_force_node << ")." << std::endl
              << "\tRatio of max node speed to local stretching wave speed = "
              << max_relative_speed / settings.getCore().getEquilibriumSpeedScale()
              << " (node " << max_relative_speed_node << ")." << std::endl;

    return status;
}

void Simulation::setup_reached_equilibrium() {
    simulation_status = Dialling;
    is_equilibrium_seeked = false;
    time_phase = 0;
    phase_counter += 1;
    settings.useDiallingDamping();
    if (phase_counter <= dial_in_phases.size() - 2) {
        std::cout << "New dialling in phase beginning, from value "
            << dial_in_phases[phase_counter] << ", to "
            << dial_in_phases[phase_counter + 1] << std::endl;
    }
}

void Simulation::advance_physics() {
    for (auto& slide : slides) { slide.update(settings.getTimeStepSize(), simulation_status == Dialling); }
    for (auto& cone : cones) { cone.update(settings.getTimeStepSize(), simulation_status == Dialling); }
}

void Simulation::check_if_equilibrium_search_begun(int stage_counter) {
    // Check if still in dialling in phase or whether it is time to wait for
    // equilibrium.
    if (time_phase >= settings.getDurationPhase() &&
        simulation_status == Dialling) {
        is_equilibrium_seeked = true;
        begin_equilibrium_search(stage_counter);
    }
}

bool Simulation::isDataPrinted() {
    return (step_count % settings.getStepPrintInterval() == 0 &&
        settings.getStepPrintInterval() > 0);
}

void Simulation::setup_tensor_increment(int stage_counter) {
    if (stage_counter == initial_stage &&
        ansatz_filename != "no_ansatz_file") { run_ansatz(stage_counter); } else {
        // Reset the variables that control each dialling/waiting process.
        simulation_status = Dialling;
        time_equilibriation = 0;
        time_phase = 0;
        dial_in_factor = 0;
        phase_counter = 0;
    }

    /* Handle cases where settings.DialInStepTime is zero or very small. This is
    a slightly hacky way to ensure that no dialling actually occurs, and the
    simulation jumps to a status = waitingForEquilibrium state. The magic number is
    just chosen to be recognisable for debugging.*/
    if (time_global < settings.getTimeStepSize()) {
        time_global = 0;
        //        time_phase = 0;
    }
}

/**
 * Does equilibrium check if it should be done in this step.
 * @param stage_counter
 * @param duration_us
 */
void Simulation::equilibriumTest(int stage_counter, long long duration_us) {
    /* If last check for equilibrium or last reaching of a new DialInFactor
 value was more than TimeBetweenEquilChecks ago, check for equilibrium.
 Also print total stretching and bending energies.*/
    if (simulation_status == WaitingForEquilibrium &&
        time_equilibriation >= settings.getTimeBetweenEquilibriumChecks()) { check_for_equilibrium(); }

    /* If equilibrium reached, write output data to file, and move to next
    'dialling in' phase. */
    if (simulation_status == EquilibriumReached) {
        save_and_print_details(stage_counter, duration_us);
        setup_reached_equilibrium();
    }
}

void Simulation::advance_dynamics() {
    const double dt = settings.getTimeStepSize();
    const double gradient_descent_multiplier =  settings.getCore().getDensity() / settings.getDampingFactor();
    const double force_limit = 1e5 * settings.getForceScale();
#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) {
        if (isForceThresholdExceeded(nodes[i], force_limit)) {
            logForceThresholdExceeded(nodes[i], triangles, settings, step_count);
        }

        /* Set velocities based on either Gradient Descent (over dampened) dynamics
        or Newtonian dynamics. Then advance positions accordingly. */
        if (!settings.getCore().isGradientDescentDynamics())  {
            nodes[i].advanceDynamics(dt);
        } else {
            // Gradient Descent dynamics.
            nodes[i].velocity = gradient_descent_multiplier * nodes[i].force / nodes[i].mass;
            nodes[i].position += dt * nodes[i].velocity;
        }
    }
}

/**
 * Runs a single tensor increment based on tensors provided in the input file.
 * @param stage_counter
 */
void Simulation::run_tensor_increment(int stage_counter) {
    setup_tensor_increment(stage_counter);
    updateProgrammedValues(stage_counter + 1);
    configure_smallest_tau();

    std::cout << "\nBeginning dynamical evolution.\n" << std::endl;

    std::vector<Eigen::Vector3d> nodeUnstressedConePosits(nodes.size());

    while (phase_counter < dial_in_phases.size() - 1) {
        try {
            // std::vector<std::chrono::system_clock::time_point> timestamps = {std::chrono::high_resolution_clock::now()};
            if (step_count == 0) { first_step_configuration(); }
            check_if_equilibrium_search_begun(stage_counter);
            if (simulation_status == Dialling) { update_dial_in_factor(); }
            long long duration_us = progress_single_step(stage_counter);

            // timestamps.push_back(std::chrono::high_resolution_clock::now());
            if (isDataPrinted()) { save_and_print_details(stage_counter, duration_us); }
            if (is_equilibrium_seeked) { equilibriumTest(stage_counter, duration_us); }

            // timestamps.push_back(std::chrono::high_resolution_clock::now());
            advance_dynamics();
            // timestamps.push_back(std::chrono::high_resolution_clock::now());
            advance_physics();
            advance_time();
            // timestamps.push_back(std::chrono::high_resolution_clock::now());
            // std::cout << std::endl << "Timings for step " << step_count << std::endl;
            // for (int i = 0; i < timestamps.size() - 1; i++) {
            //     std::cout << "d" << i + 1 << ": " << std::chrono::duration_cast<std::chrono::microseconds>(timestamps[i+1] - timestamps[i]).count() << "us" << std::endl;
            // }
            // if (step_count == 30) throw std::runtime_error("timing complete");
        } catch (std::runtime_error& error) {
            std::cout << "Error occurred at step " << step_count << ". Saving the output and stopping." << std::endl;
            save_and_print_details(stage_counter, 0);
            throw error;
        }
    }
}

void Simulation::advance_time() {
    time_global += settings.getTimeStepSize();
    time_equilibriation += settings.getTimeStepSize();
    time_phase += settings.getTimeStepSize();
    step_count++;
}


void Simulation::setInitialTriangleElongations() {
#pragma omp parallel for
    for (int i = 0; i < triangles.size(); i++) {
        triangles[i].setLocalElongation(triangles[i].next_programmed_metric_info(1));
    }
}

int Simulation::run_simulation() {
    std::cout << std::scientific << std::setprecision(8);
    std::ofstream forceDistFile;

    if (settings.getCore().isLceModeEnabled()) { setInitialTriangleElongations(); }
    /* Loop over the sequence of programmed tensors, dialling-in and waiting for
    equilibrium between each pair in the sequence. This loop is redundant in most use
    cases for this code, where only a single set of programmed tensors is supplied.*/
    for (int stage_counter = (int)initial_stage; stage_counter < stage_count - 1; stage_counter++) {
        run_tensor_increment(stage_counter);

        if (stage_count > 2 && stage_counter < stage_count - 2) {
            std::cout << "Moving on to next set of programmed tensors in sequence." << std::endl;
        }
    }

    std::cout << "Reached simulation time = " << time_global << " using " << step_count << " time steps" << std::endl;
    std::cout << "Simulation finished successfully." << std::endl;

#ifdef TIMING
    std::cout << std::endl;
    std::cout << "Timing results below." << std::endl;

    int timings_count = mechanics_times.size();
    std::pair<double, double> timings_mechanics = mean_dev(mechanics_times);
    std::pair<double, double> timings_export_prep = mean_dev(pre_export_times);
    std::pair<double, double> timings_export = mean_dev(export_times);

    std::cout << "Average timings for " << timings_count << " exports and " << nodes.size() << " nodes and " << triangles.size() << " triangles." << std::endl;
    std::cout << "Mechanics: " << static_cast<int>(timings_mechanics.first) << " +- " << static_cast<int>(
        timings_mechanics.second) << "us." << std::endl;
    std::cout << "Export preparation: " << static_cast<int>(timings_export_prep.first) << " +- " << static_cast<int>(
        timings_export_prep.second) << "us." << std::endl;
    std::cout << "Export writing: " << static_cast<int>(timings_export.first) << " +- " << static_cast<int>(
        timings_export.second) << "us." << std::endl;
#endif
    std::cout.flush();

    std::cout.rdbuf(cout_buf);

    auto duration = std::chrono::high_resolution_clock::now() - start_time;
    benchmarking_full_duration = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();
    return 0;
}

Simulation::Simulation(int argc, char** argv, int threads) {
    start_time = std::chrono::high_resolution_clock::now();
    init(argc, argv, threads);
}
