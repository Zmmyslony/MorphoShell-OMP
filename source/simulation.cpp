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
//#include <vtk/PolyData>

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
    log_stream.setOutputFileName(output_dir_name + "/log.txt");
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


//void Simulation::save_vtk() {
//    vtkNew<vtkPoints> nodes_vtk;
//    for (auto &node : nodes) {
//        nodes_vtk -> InsertNextPoint(node.pos(0), node.pos(1), node.pos(1));
//    }
//
//    vtkNew<vtkCellArray> triangles_vtk;
//    for (auto &triangle : triangles) {
//        vtkNew<vtkTriangle> triangle_vtk;
//        for (int i = 0; i < 3; i++) {
//            tirangle_vtk -> GetPointIds() -> SetId(i, triangle.vertexLabels(i));
//        }
//        triangles_vtk -> InsertNextCell(triangle_vtk);
//    }
//    vtkNew<vtkPolyData> polydata;
//    polydata -> SetPoints(nodes_vtk);
//    polydata -> SetPolys(triangles_vtk);
//}


void Simulation::setup_filenames(int argc, char *argv[]) {
    init_string += "Simulation run began at: " + getRealTime() + "\n";
    init_string += "The command that was run was:\n";
    for (int i = 0; i < argc; ++i) {
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
        if (current_path.extension() == ".vtk") {
            vtk_paths.emplace_back(current_path);
        } else if (current_path.extension() == ".cfg") {
            config_paths.emplace_back(current_path);
        } else {
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
    if (vtk_paths.size() == 2) {
        ansatz_filename = vtk_paths.back().string();
    }
    if (config_paths.empty()) { throw std::runtime_error("Missing configuration .cfg file."); }

    settings_filename = argv[1];

    std::string settings_file_name_str_final_piece = extract_Just_Filename(settings_filename);
    std::string initial_data_file_name_str_final_piece = extract_Just_Filename(initialisation_filename);
    std::string ansatz_data_file_name_str_final_piece("no_ansatz_file");

    if (argc == 4) {
        ansatz_filename = argv[3];
        ansatz_data_file_name_str_final_piece = extract_Just_Filename(ansatz_filename);
        init_string += "Using ansatz data file " + ansatz_data_file_name_str_final_piece + "\n";
    }

    output_dir_name = directory_setup(init_string, config_paths, vtk_paths).string();

    setup_logstream();

    log_stream.open();
    log_stream << init_string << std::endl;
    log_stream.close();
}


void Simulation::read_settings_new(int argc, char *argv[]) {
    std::vector<fs::path> paths(argc);
    for (int i = 0; i < argc; i++) {
        paths[i] = argv[i];
    }
    settings_new = SettingsNew(paths);
}


void Simulation::read_vtk_data() {
    log_stream.open();
    log_stream << "Now attempting to read data files. An error here likely \n"
                  "implies a problem with a data file, for example a mismatch between the \n"
                  "number of nodes or triangles stated and the number actually present; \n"
                  "or other similar mismatches in numbers of data values; or a format problem. \n"
                  "Remember the input files must have exactly the correct format." << std::endl;
    log_stream.close();
    std::vector<std::vector<Eigen::Vector3d>> programmed_metric_infos;
    std::vector<std::vector<Eigen::Matrix<double, 2, 2> >> inverted_programmed_metrics;
    std::vector<std::vector<double>> programmed_taus;
    std::vector<std::vector<Eigen::Matrix<double, 2, 2> >> programmed_second_fundamental_forms;

    readVTKData(nodes, triangles, programmed_metric_infos, inverted_programmed_metrics, programmed_taus,
                programmed_second_fundamental_forms, settings_new.getCore().isLceModeEnabled(),
                initialisation_filename,
                initial_stage,
                dialInFactorToStartFrom, nodeAnsatzPositions, ansatz_filename, log_stream);

    stage_count = inverted_programmed_metrics.size();

    // The first set of tensors are populated separately
    for (int i = 1; i < inverted_programmed_metrics.size(); i++) {
        for (int j = 0; j < triangles.size(); j++) {
            triangles[j].programmed_metric_infos.emplace_back(programmed_metric_infos[i][j]);
            triangles[j].programmed_metric_inv.emplace_back(inverted_programmed_metrics[i][j]);
            triangles[j].programmed_taus.emplace_back(programmed_taus[i][j]);
            triangles[j].programmed_second_fundamental_form.emplace_back(programmed_second_fundamental_forms[i][j]);
        }
    }

    num_nodes = nodes.size();
    num_triangles = triangles.size();
}


void Simulation::configure_nodes() {
    log_stream.open();
    log_stream << "Number of nodes = " << num_nodes << std::endl;
    log_stream << "Number of triangles = " << num_triangles << std::endl;

    if (num_triangles < 50) {
        std::cerr << "Your mesh has a small number of triangles. \nBeware that the code "
                     "is likely to be less accurate in this case, \nand unforeseen bugs are more "
                     "likely in extreme cases." << std::endl;
    }
    log_stream.close();


    for (int i = 0; i < num_nodes; ++i) {
        nodes[i].nodeLogStream.setOutputFileName(log_filename);
    }
}


void Simulation::configure_topological_properties() {
    calcTrianglesIncidentOnNodes(nodes, triangles);
    num_edges = calcTriangleAdjacencies_And_Edges(nodes, triangles, edges);

    for (int i = 0; i < num_edges; ++i) {
        edges[i].edgeLogStream.setOutputFileName(log_filename);
    }

    /* Calculate and print the number of boundary and non-boundary edges.
    Also label the nodes connected by boundary edges as boundary nodes.
    Also calculate the total initial perimeter of the sample, to be used as a
    characteristic sample length in estimating characteristic times, time steps
    etc.*/
    int numBoundaryEdges = 0;
    std::vector<double> initBoundaryEdgeLengths(num_edges);

    for (int i = 0; i < num_edges; ++i) {
        if (edges[i].isOnBoundary) {

            numBoundaryEdges += 1;

            nodes[edges[i].nodeLabels(0)].isOnBoundary = true;
            nodes[edges[i].nodeLabels(1)].isOnBoundary = true;

            //If chosen in settings, clamp whole boundary in addition to clamp
            //indicators from data file.
            if (settings_new.getCore().isBoundaryClamped()) {
                nodes[edges[i].nodeLabels(0)].isClamped = true;
                nodes[edges[i].nodeLabels(1)].isClamped = true;
            }

            // Store initial length of this boundary edge
            initBoundaryEdgeLengths[i] = (nodes[edges[i].nodeLabels(0)].pos - nodes[edges[i].nodeLabels(1)].pos).norm();
        }
    }

    double initPerimeter = kahanSum(initBoundaryEdgeLengths);
    characteristic_long_length = initPerimeter;
    log_stream.open();
    log_stream << "Initial perimeter = " << initPerimeter << std::endl;
    log_stream.close();


    if (3 * num_triangles != 2 * num_edges - numBoundaryEdges) {
        throw std::runtime_error(
                "Something has gone wrong in calculating triangle adjacencies and/or edges: the current edge and triangle counts violate a topological identity.");
    }

    log_stream.open();
    log_stream << "Number of edges = " << num_edges << std::endl;
    log_stream << "Number of boundary edges = " << numBoundaryEdges << std::endl;
    log_stream << "Number of non-boundary edges = " << num_edges - numBoundaryEdges << std::endl;
    log_stream.close();
}


void Simulation::configure_triangles() {
    for (int i = 0; i < num_triangles; ++i) {
        triangles[i].triLogStream.setOutputFileName(log_filename);
    }

    int numBoundaryTriangles = 0;
    for (int i = 0; i < num_triangles; ++i) {
        if (triangles[i].isOnBoundary) {
            numBoundaryTriangles += 1;
        }
    }

    log_stream.open();
    log_stream << "Number of boundary triangles = " << numBoundaryTriangles << std::endl;
    log_stream << "Number of holes in mesh = " << 1 + num_edges - num_nodes - num_triangles
               << std::endl; // From Euler's formula for a planar graph.
    log_stream.close();

    if (1 + num_edges - num_nodes - num_triangles < 0) {
        throw std::runtime_error(
                "Something is very wrong with the mesh, because the code thinks it has a negative number of holes! A first thing to check is that all nodes touch at least one tri.");
    }
}

void Simulation::set_node_patches() {
    log_stream.open();
    log_stream << "\n" << "Beginning patch selection and related pre-calculations." << std::endl;
    log_stream.close();

    calc_nonVertexPatchNodes_and_MatForPatchDerivs(nodes, triangles, log_stream,
                                                   settings_new.getCore().getPatchMatrixThreshold());

    log_stream.open();
    log_stream << "Successfully completed patch setup." << "\n" << std::endl;
    log_stream.close();
}


void Simulation::orient_node_labels() {
    updateTriangleValues(nodes, triangles, simulation_status, -12345, 98765, settings_new);
    Eigen::Vector3d tempZAxisVec;
    tempZAxisVec << 0.0, 0.0, 1.0;
    for (int i = 0; i < num_triangles; ++i) {
        if (tempZAxisVec.dot(triangles[i].currSides.col(0).cross(triangles[i].currSides.col(1))) < 0) {
            int tempLabel = triangles[i].vertexLabels(2);
            triangles[i].vertexLabels(2) = triangles[i].vertexLabels(1);
            triangles[i].vertexLabels(1) = tempLabel;
        }
    }
    std::cout << "CHECK THIS - should shuffle normals before patch matrix calc I think?" << std::endl;
    updateTriangleValues(nodes, triangles, simulation_status, -12345, 98765, settings_new);
}


void Simulation::set_initial_conditions() {
    setRemainingInitCond_and_NodeMasses(nodes, triangles, edges,
                                        settings_new);
}

void Simulation::find_smallest_element() {
    characteristic_short_length = std::min_element(triangles.begin(), triangles.end(),
                                                   [](Triangle &first, Triangle &second) {
                                                       return first.getLinearSize() < second.getLinearSize();
                                                   })->getLinearSize();
    // This is slightly incorrect as it takes the smallest element and largest tau, instead of taking
    // smallest ratio of size to tau, though it is just erring on the side of caution.
    std::vector<double> largest_tau_vector(stage_count);
    double largest_tau = DBL_MIN;
    for (int i = 0; i < triangles.size(); i++) {
        for (auto &tau: triangles[i].programmed_taus) {
            if (tau < largest_tau) { largest_tau = tau; }
        }
    }
    characteristic_length_over_tau = characteristic_short_length / sqrt(largest_tau);

    log_stream.open();
    log_stream << "Sheet thickness = " << settings_new.getCore().getThickness() << std::endl;
    log_stream << "Approx smallest element linear size = " << characteristic_short_length << std::endl;
    log_stream.close();

}


void Simulation::setup_characteristic_scales() {
    std::vector<double> initAreas(triangles.size());
    for (int i = 0; i < triangles.size(); ++i) {
        initAreas[i] = triangles[i].initArea;
    }
    double initial_area = kahanSum(initAreas);
    settings_new.SetupCharacteristicSizes(initial_area, characteristic_short_length);
}


void Simulation::setup_equilibrium_dial_in_factors() {
    double tempDialInFactor = 0;
    if (settings_new.getCore().getDialInResolution() > 0 && settings_new.getDurationPhase() >= 0) {
        while (tempDialInFactor < 1) {
            dial_in_phases.push_back(tempDialInFactor);
            tempDialInFactor += settings_new.getCore().getDialInResolution();
        }
        dial_in_phases.push_back(1);
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
    gaussCurvatures = std::vector<double>(num_triangles, DBL_MAX);
    meanCurvatures = std::vector<double>(num_triangles, DBL_MAX);
    stretchEnergies = std::vector<double>(num_triangles, DBL_MAX);
    bendEnergies = std::vector<double>(num_triangles, DBL_MAX);
    stretchEnergyDensities = std::vector<double>(num_triangles, DBL_MAX);
    bendEnergyDensities = std::vector<double>(num_triangles, DBL_MAX);
    kineticEnergies = std::vector<double>(num_nodes, DBL_MAX);
    strainMeasures = std::vector<double>(num_triangles, DBL_MAX);
    cauchyStressEigenvals = std::vector<Eigen::Vector2d>(num_triangles);
    cauchyStressEigenvecs = std::vector<Eigen::Matrix<double, 3, 2>>(num_triangles);

    angleDeficits = std::vector<double>(num_nodes, DBL_MAX);
    interiorNodeAngleDeficits = std::vector<double>(num_nodes, DBL_MAX);
    boundaryNodeAngleDeficits = std::vector<double>(num_nodes, DBL_MAX);
}


void Simulation::init(int argc, char *argv[]) {
    setup_filenames(argc, argv);
    read_settings_new(argc, argv);
    read_vtk_data();

    configure_nodes();
    configure_topological_properties();
    configure_triangles();
    configureNodeAdjacency(nodes, edges);
    set_node_patches();
    orient_node_labels();

    set_initial_conditions();
    find_smallest_element();

    log_stream.open();
    log_stream << settings_new.SetupDialInTime(characteristic_long_length);
    log_stream << settings_new.SetupStepTime(characteristic_short_length);
    log_stream << settings_new.SetupPrintFrequency();
    log_stream.close();

    setup_characteristic_scales();
    setup_equilibrium_dial_in_factors();

    initialise_simulation_vectors();

    if (settings_new.getCore().isInitialPositionsPerturbed()) {
        perturbInitialPositionsWithRandomNoise(nodes, characteristic_short_length);
    }
}


void Simulation::run_ansatz(int counter) {
    for (int i = 0; i < num_nodes; ++i) {
        nodes[i].pos = nodeAnsatzPositions[i];
    }

    time_equilibriation = 0;
    simulation_status = Dialling;
    phase_counter = 0;

    if (settings_new.getCore().isAnsatzMetricUsed()) {
        log_stream.open();
        log_stream
                << "\nINFO Assuming ansatz state to be relaxed (initial programmed tensors calculated from geometry)."
                << std::endl;
        log_stream.close();

        dial_in_factor = 0.0;
        time_phase = 0.0;

        // Calculate all necessary geometry for the ansatz state.
        updateTriangleValues(nodes, triangles, WaitingForEquilibrium, dial_in_factor,
                             counter, settings_new);
        updateSecondFundamentalForms(triangles, settings_new.getCore());


        // Temp LU decomp of triangle metric, used to check invertibility.
        Eigen::FullPivLU<Eigen::Matrix<double, 2, 2>> tempMetricDecomp;


        /* Alter inverted_programmed_metrics[initial_stage] and similar to change where
        the programmed quantities are dialling from.*/

        for (int i = 0; i < num_triangles; ++i) {
            // Test for invertibility of metric before taking inverse.
            tempMetricDecomp.compute((triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse());
            if (!tempMetricDecomp.isInvertible()) {
                throw std::runtime_error(
                        "At least one triangle had a non-invertible metric in the ansatz state. \n"
                        "This should not occur in a reasonable mesh. Aborting.");
            }

            if (settings_new.getCore().isFirstTensorSkipped()) {
                triangles[i].programmed_metric_inv[initial_stage + 1] = (
                        triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                triangles[i].programmed_taus[initial_stage + 1] =
                        triangles[i].programmed_taus[initial_stage + 1];
                triangles[i].programmed_second_fundamental_form[initial_stage +
                                                                1] = triangles[i].secFF;
            } else {
                triangles[i].programmed_metric_inv[initial_stage] = (
                        triangles[i].defGradient.transpose() * triangles[i].defGradient).inverse();
                triangles[i].programmed_taus[initial_stage] =
                        triangles[i].programmed_taus[initial_stage + 1];
                triangles[i].programmed_second_fundamental_form[initial_stage] = triangles[i].secFF;
            }

        }
    } else {
        dial_in_factor = dialInFactorToStartFrom;
        for (int i = 0; i < dial_in_phases.size(); ++i) {
            if (dial_in_phases[i] < dial_in_factor) {
                phase_counter = i;
            }
        }
        double dial_in_since_last_phase = dial_in_factor - dial_in_phases[phase_counter];
        double dial_in_phase_duration = dial_in_phases[phase_counter + 1] - dial_in_phases[phase_counter];

        time_phase = dial_in_since_last_phase / dial_in_phase_duration * settings_new.getDurationPhase();

        if (settings_new.getCore().isFirstTensorSkipped()) {
            for (int i = 0; i < num_triangles; ++i) {
                triangles[i].programmed_metric_infos[initial_stage] =
                        triangles[i].programmed_metric_infos[initial_stage + 1];
                triangles[i].programmed_second_fundamental_form[initial_stage] =
                        triangles[i].programmed_second_fundamental_form[initial_stage + 1];
//                programmed_second_fundamental_forms[initial_stage][i] =
//                        programmed_second_fundamental_forms[initial_stage + 1][i];
//                programmed_metric_infos[initial_stage][i] =
//                        programmed_metric_infos[initial_stage + 1][i];
            }
        }
    }
}

// Seide, 1956, Axisymmetrical Buckling of Circular Cones Under Axial Compression,
// https://doi.org/10.1115/1.4011410
void Simulation::setup_imposed_seide_deformations(double &s1, int highest_node, int lowest_node,
                                                  std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    lambda = 0.9;
    double cone_angle = asin(pow(lambda, 1.5));
    log_stream.open();
    log_stream << "IMPOSING SEIDE DEFORMATIONS." << std::endl;
    log_stream.close();

    // Use nodes[i].isSeideDisplacementEnabled to label nodes whose
    // positions we will force to be those of Seide's setup (found
    // by findng the membrane stress solution for his base state, and
    // integrating the strains etc).
    double intermLengthScaleUpper = 3.0 * sqrt(settings_new.getCore().getThickness() * 0.18);
    double intermLengthScaleLower = 3.0 * sqrt(settings_new.getCore().getThickness() * 1.8);
    log_stream.open();
    log_stream
            << "Imposing Seide displacements within the following (initial) vertical distances of the top and bottom: "
            << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
    log_stream.close();
    for (int n = 0; n < num_nodes; ++n) {
        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
            nodes[n].isSeideDisplacementEnabled = true;
        }
    }

    // STORE PERFECT CONE ANSATZ
    for (int n = 0; n < num_nodes; ++n) {
        nodeUnstressedConePosits[n] = nodes[n].pos;
    }
    s1 = sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
              nodes[highest_node].pos(1) * nodes[highest_node].pos(1)) / sin(cone_angle);

//    Eigen::Vector3d testTriCurrCentroid;
//    testTriCurrCentroid = (nodes[triangles[settings.testTriangle].vertexLabels(0)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(1)].pos +
//                           nodes[triangles[settings.testTriangle].vertexLabels(2)].pos) / 3;
    //sTest = sqrt(testTriCurrCentroid(0)*testTriCurrCentroid(0) + testTriCurrCentroid(1)*testTriCurrCentroid(1)) / sin(settings.ConeAngle);

}

// TODO Add cone squashing -> like slide but using a cone.
//void Simulation::setup_glass_cones(int highest_node, int lowest_node) {
//    double cone_angle = 1.02327019;
//    settings.init_slide_z_coord_upper += -tan(cone_angle) *
//                                         sqrt(nodes[highest_node].pos(0) * nodes[highest_node].pos(0) +
//                                              nodes[highest_node].pos(1) * nodes[highest_node].pos(1));
//    settings.init_slide_z_coord_lower += -tan(cone_angle) *
//                                         sqrt(nodes[lowest_node].pos(0) * nodes[lowest_node].pos(0) +
//                                              nodes[lowest_node].pos(1) * nodes[lowest_node].pos(1));
//    log_stream.open();
//    log_stream << "USING TWO GLASS CONES FOR SQUASHING." << std::endl;
//    log_stream.close();
//
//    // Hijack nodes[i].isOnBoundary to instead label nodes whose
//    // forces we will modify to kill any components not tangential to
//    // a perfect cone base state. We do this to nodes within an intermediate
//    // distance vertically from the ends of the cone.
//    double intermLengthScaleUpper = 2.0 * sqrt(settings_new.getCore().getThickness() * 0.18);
//    double intermLengthScaleLower = 2.0 * sqrt(settings_new.getCore().getThickness() * 1.8);
//    log_stream.open();
//    log_stream << "Apply normal-force-killer within the following vertical distances of the top and bottom: "
//               << intermLengthScaleUpper << ", " << intermLengthScaleLower << std::endl;
//    log_stream.close();
//    for (int n = 0; n < num_nodes; ++n) {
//        if (((nodes[highest_node].pos(2) - nodes[n].pos(2)) < intermLengthScaleUpper) ||
//            ((nodes[n].pos(2) - nodes[lowest_node].pos(2)) < intermLengthScaleLower)) {
//            nodes[n].isOnBoundary = true;
//        }
//    }
//}

//void Simulation::update_slide_properties() {
//    if (!settings.is_controlled_force_enabled) {
//        settings.upper_slide_displacement =
//                time * settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time;
//    } else { // settings.isControlledForceEnabled == true instead
//        //settings.upperSlideWeight = (settings.ShearModulus * settings.SheetThickness * settings.SheetThickness) * (time * settings.slideSpeedPrefactor / bending_long_time);
//        settings.slide_damping_param =
//                0.4 * settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness /
//                (settings.slide_speed_prefactor * settings.sample_char_length / settings.bending_long_time);
//        if (fabs(settings.upper_tot_slide_force + settings.upper_slide_weight) /
//            (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness) <
//            settings.total_slide_force_to_mu_t_sq_ratio_equil_threshold
//            && settings.const_slide_weight_fac < 0) {
//            if (time_equilibriation > settings.time_between_equil_checks) {
//                if (equilibriumCheck(nodes, triangles, settings, log_stream) == EquilibriumReached) {
//                    settings.is_slide_just_equilibrated = 1;
//                    settings.upper_slide_weight +=
//                            settings.slide_weight_dial_speed_fac *
//                            (settings.time_step_size / settings.bending_long_time) *
//                            (settings.shear_modulus * settings.sheet_thickness * settings.sheet_thickness);
//                }
//                time_equilibriation = 0.0;
//            }
//        }
//
//        // To do constant-weight slide instead
//        if (settings.const_slide_weight_fac > 0) {
//            settings.upper_slide_weight =
//                    settings.const_slide_weight_fac * settings.shear_modulus * settings.sheet_thickness *
//                    settings.sheet_thickness;
//        }
//    }
//    settings.curr_slide_z_coord_upper = settings.init_slide_z_coord_upper - settings.upper_slide_displacement;
//}


//void Simulation::impose_seide_deformation(double s1, const std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
//    if (settings.is_seide_deformations_enabled) {
//        double pInit = 3.0 * (settings.shear_modulus * settings.sheet_thickness *
//                              settings.sheet_thickness); // So we don't have to start all the way from p=0, chosen based on previous sims.
//        settings.p = pInit + time * settings.p_speed_prefactor * settings.shear_modulus * settings.sheet_thickness *
//                             settings.sheet_thickness / settings.bending_long_time;
//
//        for (int n = 0; n < num_nodes; ++n) {
//            if (nodes[n].isSeideDisplacementEnabled || step_count == 0) {
//
//                double tCone = settings.sheet_thickness / sqrt(settings.lambda);
//
//                double polarAng = atan2(nodeUnstressedConePosits[n](1), nodeUnstressedConePosits[n](0));
//                double s = sqrt(nodeUnstressedConePosits[n](0) * nodeUnstressedConePosits[n](0) +
//                                nodeUnstressedConePosits[n](1) * nodeUnstressedConePosits[n](1)) /
//                           sin(settings.cone_angle);
//                Eigen::Vector3d uHat;
//                uHat << sin(settings.cone_angle) * cos(polarAng), sin(settings.cone_angle) * sin(polarAng), -cos(
//                        settings.cone_angle);
//                Eigen::Vector3d wHat;
//                wHat << -cos(settings.cone_angle) * cos(polarAng), -cos(settings.cone_angle) *
//                                                                   sin(polarAng), -sin(settings.cone_angle);
//                double u = -settings.p * log(s / s1) /
//                           (2.0 * M_PI * settings.youngs_modulus * tCone * sin(settings.cone_angle) *
//                            cos(settings.cone_angle));
//                double w = -settings.p * (log(s / s1) + settings.poisson_ratio) /
//                           (2.0 * M_PI * settings.youngs_modulus * tCone * cos(settings.cone_angle) *
//                            cos(settings.cone_angle));
//                nodes[n].pos = nodeUnstressedConePosits[n] + u * uHat + w * wHat;
//            }
//        }
//    }
//}


void Simulation::first_step_configuration(double &seide_quotient,
                                          std::vector<Eigen::Vector3d> &nodeUnstressedConePosits) {
    slides = settings_new.getSlides();
    std::vector<Eigen::Vector3d> node_positions(nodes.size());
    for (int i = 0; i < nodes.size(); i++) {
        node_positions[i] = nodes[i].pos;
    }

    for (auto &slide: slides) {
        slide.initialise(node_positions, settings_new.getDurationPhase());
    }
}

void Simulation::begin_equilibrium_search(int counter) {
    dial_in_factor = dial_in_phases[phase_counter + 1];
    updateTriangleValues(nodes, triangles, simulation_status, dial_in_factor,
                         counter, settings_new);
    simulation_status = WaitingForEquilibrium;
    // NB an EquilCheck has not actually just occurred, but this has the
    //desired effect of ensuring that each DialInFactor value is held
    //for at least one TimeBetweenEquilChecks.
    time_equilibriation = 0.0;

    settings_new.useEquilibriumDamping();

    log_stream.open();
    log_stream << "\tReached dial-in factor of " << dial_in_phases[phase_counter + 1]
               << ". Waiting for equilibrium." << std::endl;
    log_stream.close();
}

//void Simulation::update_second_fundamental_form() {
//
//}

void
Simulation::add_elastic_forces(const std::vector<std::vector<std::pair<int, int>>> &correspondingTrianglesForNodes) {
    double stretchingPreFac = 0.5 * settings_new.getCore().getThickness() * settings_new.getCore().getShearModulus();
    double bendingPreFac =
            0.5 * pow(settings_new.getCore().getThickness(), 3) * settings_new.getCore().getShearModulus() / 12;
    double JPreFactor = settings_new.getCore().getGentFactor() / pow(settings_new.getCore().getThickness(), 2);
    std::vector<Eigen::Vector3d> forcesForEachTriangle(6 * triangles.size());

#pragma omp parallel
    {
#pragma omp for
        // Calculates forces experienced by nodes coming from each triangle.
        for (int i = 0; i < triangles.size(); i++) {
            triangles[i].updateSecondFundamentalForm(bendingPreFac, JPreFactor,
                                                     settings_new.getCore().getPoissonRatio());
            triangles[i].updateHalfPK1Stress(stretchingPreFac);
            Eigen::Matrix<double, 3, 3> stretchForces = triangles[i].getStretchingForces();

            Eigen::Matrix<double, 3, 3> triangleEdgeNormals = triangles[i].getTriangleEdgeNormals();
            Eigen::Matrix<double, 3, 3> normalDerivPiece =
                    0.5 * triangles[i].currAreaInv * (triangles[i].patchSecDerivs.transpose() * triangleEdgeNormals);

            for (int n = 0; n < 3; ++n) {
                forcesForEachTriangle[6 * i + n] =
                        triangles[i].getBendingForce(normalDerivPiece, n) + stretchForces.col(n);
            }
            for (int n = 3; n < 6; ++n) {
                forcesForEachTriangle[6 * i + n] = triangles[i].getBendingForce(normalDerivPiece, n);
            }
        }

#pragma omp for
        // Applies forces to each node.
        for (int i = 0; i < nodes.size(); i++) {
            for (auto &trianglesForNode: correspondingTrianglesForNodes[i]) {
                int index = 6 * trianglesForNode.first + trianglesForNode.second;
                nodes[i].force += forcesForEachTriangle[index];
            }
        }
    }
};


void Simulation::add_non_elastic_forces() {
#pragma omp parallel for
    for (int i = 0; i < nodes.size(); i++) {
        nodes[i].add_damping(settings_new);
        nodes[i].add_gravity(settings_new.getGravity());
//        nodes[i].add_prod_force(settings_new);
//        nodes[i].add_load_force();

        for (auto &slide: slides) {
            slide.addInteractionForce(nodes[i].pos,
                                      nodes[i].force,
                                      settings_new.getCore().getShearModulus(),
                                      settings_new.getCore().getThickness());
        }

        nodes[i].apply_boundary_conditions();
    }
}

void Simulation::progress_single_step(int counter,
                                      const std::vector<std::vector<std::pair<int, int>>> &correspondingTrianglesForNodes) {

    updateTriangleValues(nodes, triangles, simulation_status, dial_in_factor,
                         counter, settings_new);



    /* Calculate secFF estimates for triangles, and related quantities such
    as the derivative of the bending energy wrt the secFF components.*/
//    updateSecondFundamentalForms(triangles, settings_new.getCore());

    // Calculate current strain and bending force on each node.
//    update_elastic_forces(nodes, triangles, settings_new.getCore(), correspondingTrianglesForNodes);

    add_elastic_forces(correspondingTrianglesForNodes);
    /* Add force contributions from e.g. damping, loads, 'prod' perturbation, and
     account for BCs e.g. clamping. */

    add_non_elastic_forces();

}

void Simulation::update_dial_in_factor() {
    double starting_dial_in_factor = dial_in_phases[phase_counter];
    double final_dial_in_factor = dial_in_phases[phase_counter + 1];
    double relative_time = time_phase / settings_new.getDurationPhase();
    dial_in_factor = starting_dial_in_factor + (final_dial_in_factor - starting_dial_in_factor) * relative_time;
}


void Simulation::save_and_print_details(int counter, double duration_us) {
    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings_new.getCore());
    if (settings_new.getCore().isEnergyPrinted()) {
        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings_new.getCore());
    }

    // Here we hack the angleDeficits to instead tell us whether a node has a Seide displacement imposed or not.
    for (int n = 0; n < num_nodes; ++n) {
        if (nodes[n].isSeideDisplacementEnabled) {
            angleDeficits[n] = 1;
        } else {
            angleDeficits[n] = 0;
        }
    }
    writeVTKDataOutput(nodes, triangles, step_count, time_global, dial_in_factor, counter,
                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                       cauchyStressEigenvals, cauchyStressEigenvecs, settings_new, output_dir_name);

    log_stream.open();
    log_stream << log_prefix()
               << "Last step's execution time " << duration_us << " us."
               << std::endl;

    log_stream.close();
}

std::string Simulation::log_prefix() const {
    std::stringstream ss;
    ss << std::setprecision(3) << std::fixed << getRealTime() << ", step = " << step_count << ", time = " << time_global
       << ", dial-in factor = " <<
       dial_in_factor << ": ";
    return ss.str();
}

void Simulation::error_large_force(int counter) {
    log_stream.open();
    log_stream << "At " << getRealTime() << ", stepCount = " << step_count << ", time = " << time_global
               << ", dial-in factor = " << dial_in_factor
               << " a force was suspiciously large, there is probably a problem. Writing VTK output and then aborting."
               << std::endl;
    log_stream.close();

    calcCurvatures(nodes, triangles, gaussCurvatures, meanCurvatures, angleDeficits,
                   interiorNodeAngleDeficits, boundaryNodeAngleDeficits, settings_new.getCore());
    if (settings_new.getCore().isEnergyPrinted()) {
        calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                                stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                                cauchyStressEigenvals, cauchyStressEigenvecs, settings_new.getCore());
    }
    writeVTKDataOutput(nodes, triangles, step_count, time_global, dial_in_factor, counter,
                       gaussCurvatures, meanCurvatures, angleDeficits, interiorNodeAngleDeficits,
                       boundaryNodeAngleDeficits, stretchEnergyDensities, bendEnergyDensities,
                       stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                       cauchyStressEigenvals, cauchyStressEigenvecs, settings_new, output_dir_name);
    throw std::runtime_error("Unexpectedly large force at step " + std::to_string(step_count) + ".");
}

void Simulation::check_for_equilibrium() {
    log_stream.open();
    log_stream << log_prefix() << "Checking for equilibrium. "
               << std::endl;
    log_stream.close();

    simulation_status = equilibriumCheck(nodes, settings_new, triangles, log_stream);
    time_equilibriation = 0.0;

    calcEnergiesAndStresses(nodes, triangles, stretchEnergyDensities, bendEnergyDensities,
                            stretchEnergies, bendEnergies, kineticEnergies, strainMeasures,
                            cauchyStressEigenvals, cauchyStressEigenvecs, settings_new.getCore());
}

void Simulation::setup_reached_equilibrium() {
    simulation_status = Dialling;
    time_phase = 0;
    phase_counter += 1;
    settings_new.useDiallingDamping();
    if (phase_counter <= dial_in_phases.size() - 2) {
        log_stream.open();
        log_stream << "New dialling in phase beginning, from value "
                   << dial_in_phases[phase_counter] << ", to "
                   << dial_in_phases[phase_counter + 1] << std::endl;
        log_stream.close();
    }
}

void Simulation::run_tensor_increment(int stage_counter) {
    /* A third input file specifying an ansatz for the node positions may have
        been read in if it was given as a command line argument. This could for
        example correspond to an output file from this code, to carry on where some
        previous simulation left off. In this case, the nodes are moved to their
        ansatz positions here, and other relevant variables are set up. */
    if (ansatz_filename != "no_ansatz_file" &&
        stage_counter == initial_stage) {
        run_ansatz(stage_counter);
    } else {
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
    if (time_global < settings_new.getTimeStepSize()) {
        time_global = 0;
        time_phase = 0;
    }


    // Begin dynamical evolution of node positions and velocities.
    log_stream.open();
    log_stream << "\nBeginning dynamical evolution.\n" << std::endl;
//    log_stream << "\nCREATING VECTOR TO STORE UNSTRESSED CONE NODE POSITIONS.\n" << std::endl;
    log_stream.close();
    std::vector<Eigen::Vector3d> nodeUnstressedConePosits(num_nodes);
    double seide_quotient = DBL_MAX;

    std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes = getCorrespondingTrianglesForNodes(
            triangles, nodes);

    while (phase_counter < dial_in_phases.size() - 1) {
        if (step_count == 0) { first_step_configuration(seide_quotient, nodeUnstressedConePosits); }
        for (auto &slide: slides) {
            slide.update(settings_new.getTimeStepSize(), simulation_status == Dialling);
        }
//        update_slide_properties();

//        if (settings_new.getCore().isSeideDeformations()) {
//            impose_seide_deformation(seide_quotient, nodeUnstressedConePosits);
//        }

        zeroForces(nodes);

        // Check if still in dialling in phase or whether it is time to wait for
        // equilibrium.
        if (time_phase >= settings_new.getDurationPhase() &&
            simulation_status == Dialling) {
            begin_equilibrium_search(stage_counter);
        }

        /* If not waiting for equilibrium, set the current value of the Dial-In
        Factor, based on linear dialling in between the previously calculated
        'checkpoint' values. */
        if (simulation_status == Dialling) { update_dial_in_factor(); }

        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
        progress_single_step(stage_counter, correspondingTrianglesForNodes);
        auto duration = std::chrono::steady_clock::now() - begin;
        double duration_us = std::chrono::duration_cast<std::chrono::microseconds>(duration).count();


        /* Write output data regularly.
        Can be switched off with settings.PrintFrequency < 0.0.
        Doing the write-out at this point in the loop means the node positions
        and the triangle geometry data match in the output, which is desirable!*/
        if ((step_count % settings_new.getStepPrintInterval() == 0 &&
             settings_new.getStepPrintInterval() > 0)) {
            save_and_print_details(stage_counter, duration_us);
        }


        if (is_equilibrium_seeked) {
            /* If last check for equilibrium or last reaching of a new DialInFactor
             value was more than TimeBetweenEquilChecks ago, check for equilibrium.
             Also print total stretching and bending energies.*/
            if (time_equilibriation > settings_new.getTimeBetweenEquilibriumChecks() &&
                simulation_status == WaitingForEquilibrium) {
                check_for_equilibrium();
            }

            /* If equilibrium reached, write output data to file, and move to next
            'dialling in' phase. */
            if (simulation_status == EquilibriumReached) {
                save_and_print_details(stage_counter, duration_us);
                setup_reached_equilibrium();
            }
        }

        /* Advance node positions and velocities using dynamical solver. Error
        caught here if any node force is suspiciously high, indicating probable
        'blowing, up', and the code aborts in that case. */
        try { advanceDynamics(nodes, triangles, settings_new, log_stream); }
        catch (const std::runtime_error &error) {
            error_large_force(stage_counter);
        }

        advance_time();
    }
}

void Simulation::advance_time() {
    time_global += settings_new.getTimeStepSize();
    time_equilibriation += settings_new.getTimeStepSize();
    time_phase += settings_new.getTimeStepSize();
    step_count++;
}


void Simulation::run_simulation() {
    log_stream << std::scientific << std::setprecision(8);
    std::ofstream forceDistFile;

    /* Loop over the sequence of programmed tensors, dialling-in and waiting for
    equilibrium between each pair in the sequence. This loop is redundant in most use
    cases for this code, where only a single set of programmed tensors is supplied.*/
    for (std::size_t stage_counter = initial_stage; stage_counter < stage_count - 1; stage_counter++) {
        run_tensor_increment(stage_counter);

        if (stage_count > 2 && stage_counter < stage_count - 2) {
            log_stream.open();
            log_stream << "Moving on to next set of programmed tensors in sequence." << std::endl;
            log_stream.close();
        }
    }

// Print some helpful final things.
    log_stream.open();
    log_stream << "Reached simulation time = " << time_global << " using " << step_count << " time steps" << std::endl;
    log_stream.close();
}

Simulation::Simulation(int argc, char **argv) {
    init(argc, argv);
}

