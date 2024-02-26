//
// Created by Michał Zmyślony on 22/03/2022.
//
#define _USE_MATH_DEFINES

#include "Settings.hpp"

#include <map>
#include <iostream>
#include <cmath>
#include <math.h>

double *Settings::getParameterAddressDouble(const std::string &parameterName) {
    std::map<std::string, double *> doubleMap{
            {"init_slide_z_coord_lower",                     &init_slide_z_coord_lower},
            {"init_slide_z_coord_upper",                     &init_slide_z_coord_upper},
            {"curr_slide_z_coord_upper",                     &curr_slide_z_coord_upper},
            {"slide_stiffness_prefactor",                    &slide_stiffness_prefactor},
            {"slide_speed_prefactor",                        &slide_speed_prefactor},
            {"upper_slide_displacement",                     &upper_slide_displacement},
            {"slide_friction_coefficient",                   &slide_friction_coefficient},
            {"thicknesses_above_lowest_node_to_clamp_up_to", &thicknesses_above_lowest_node_to_clamp_up_to},
            {"bending_long_time",                            &bending_long_time},
            {"cone_angle",                                   &cone_angle},
            {"lambda",                                       &lambda},
            {"slide_weight_dial_speed_fac",                  &slide_weight_dial_speed_fac},
            {"total_slide_force_to_mu_t_sq_ratio_equil_threshold",
                                                             &total_slide_force_to_mu_t_sq_ratio_equil_threshold},
            {"initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq",
                                                             &initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq},
            {"slide_damping_param",                          &slide_damping_param},
            {"upper_slide_weight",                           &upper_slide_weight},
            {"upper_slide_vel",                              &upper_slide_vel},
            {"upper_tot_slide_force",                        &upper_tot_slide_force},
            {"const_slide_weight_fac",                       &const_slide_weight_fac},
            {"spacer_height",                                &spacer_height},
            {"p",                                            &p},
            {"p_speed_prefactor",                            &p_speed_prefactor},
            {"specify_init_slide_z_coord_upper",             &specify_init_slide_z_coord_upper},
            {"specify_init_slide_z_coord_lower",             &specify_init_slide_z_coord_lower},
            {"print_frequency",                              &print_frequency},
            {"num_damp_factor",                              &num_damp_factor},
            {"dial_in_damping",                              &dial_in_damping},
            {"equilibriation_damping",                       &equilibriation_damping},
            {"gent_factor",                                  &gent_factor},
            {"time_step_prefactor",                          &time_step_prefactor},
            {"time_step_size",                                    &time_step},
            {"char_force",                                   &char_force},
            {"char_speed",                                   &char_speed},
            {"dial_in_resolution",                           &dial_in_resolution},
            {"dial_in_step_time",                            &dial_in_step_time},
            {"time_between_equil_checks_prefactor",          &time_between_equil_checks_prefactor},
            {"time_between_equil_checks",                    &time_between_equil_checks},
            {"dial_in_step_time_prefactor",                  &dial_in_step_time_prefactor},
            {"prod_force_time",                              &prod_force_time},
            {"prod_strength",                                &prod_strength},
            {"load_force_time",                              &load_force_time},
            {"load_strength",                                &load_strength},
            {"poisson_ratio",                                &poisson_ratio},
            {"sample_char_length",                           &sample_char_length},
            {"char_stretch_energy_density_scale",            &char_stretch_energy_density_scale},
            {"char_stretch_energy_scale",                    &char_stretch_energy_scale},
            {"sheet_thickness",                              &sheet_thickness},
            {"shear_modulus",                                &shear_modulus},
            {"youngs_modulus",                               &youngs_modulus},
            {"init_density",                                 &init_density},
            {"approx_min_init_elem_size",                    &approx_min_init_elem_size},
            {"char_force_scale",                             &char_force_scale},
            {"patch_matrix_dimensionless_conditioning_threshold",
                                                             &patch_matrix_dimensionless_conditioning_threshold},
            {"gravity_sign",                                 &gravity_sign}
    };
    double *parameterAddress = doubleMap[parameterName];
    return parameterAddress;
}


int *Settings::getParameterAddressInt(const std::string &parameterName) {
    std::map<std::string, int *> intMap{
            {"test_triangle",              &test_triangle},
            {"is_slide_just_equilibrated", &is_slide_just_equilibrated},
            {"num_nodes",                  &num_nodes},
            {"num_triangles",              &num_triangles},
            {"num_edges",                  &num_edges},
            {"num_boundary_nodes",         &num_boundary_nodes},
            {"number_of_cores",            &number_of_cores}
    };
    int *parameterAddress = intMap[parameterName];
    return parameterAddress;
}


bool *Settings::getParameterAddressBool(const std::string &parameterName) {
    std::map<std::string, bool *> boolMap{
            {"is_dialing_disabled",                          &is_dialing_disabled},
            {"glass_cones",                                  &glass_cones},
            {"is_controlled_force_enabled",                  &is_controlled_force_enabled},
            {"is_seide_deformations_enabled",                &is_seide_deformations_enabled},
            {"is_lce_mode_enabled",                          &is_lce_mode_enabled},
            {"is_energy_densities_printed",                  &is_energy_densities_printed},
            {"is_triangle_areas_printed",                    &is_triangle_areas_printed},
            {"is_angle_deficits_printed",                    &is_angle_deficits_printed},
            {"is_gradient_descent_dynamics_enabled",         &is_gradient_descent_dynamics_enabled},
            {"for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff",
                                                             &for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff},
            {"is_dialing_from_ansatz_enabled",               &is_dialing_from_ansatz_enabled},
            {"is_perturbation_of_initial_positions_enabled", &is_perturbation_of_initial_positions_enabled},
            {"is_boundary_clamped",                          &is_boundary_clamped}
    };
    bool *parameterAddress = boolMap[parameterName];
    return parameterAddress;
}


void Settings::SetupDialIn(CustomOutStreamClass &logStream) {
    /* Calculate 'dialling in' time and damping coefficient based on toy model
stretching and bending analyses. The 'Long times' are approximate characteristic
times for the longest-wavelength modes in the system for stretching and bending.
The damping coefficient is chosen to approximately critically damp the longest-
wavelength bending mode.
*/
    const double PI = M_PI;

    logStream.open();
    double minWavevector = 2 * PI / sample_char_length;
    damping_scale =
            2 * sqrt(init_density * shear_modulus / (6 * (1 - poisson_ratio))) * sheet_thickness * minWavevector *
            minWavevector;
    num_damp_factor = dial_in_damping * damping_scale;
    logStream << "Numerical Damping Coefficient = " << num_damp_factor << std::endl;

    double stretchingLongTime = sqrt(init_density / shear_modulus) / minWavevector;
    bending_long_time = sqrt(init_density * 6 * (1 - poisson_ratio) / shear_modulus) /
                        (sheet_thickness * minWavevector * minWavevector);

    if (stretchingLongTime > bending_long_time) {
        dial_in_step_time = dial_in_step_time_prefactor * stretchingLongTime;
        logStream << "Dial-in step time = " << dial_in_step_time
                  << ", set based on stretching rather than bending." << std::endl;
    } else {
        dial_in_step_time = dial_in_step_time_prefactor * bending_long_time;
        logStream << "Dial-in step time = " << dial_in_step_time
                  << ", set based on bending rather than stretching." << std::endl;
    }
    logStream.close();
}


void Settings::SetupStepTime(CustomOutStreamClass &logStream) {
    /* Calculate time step based on toy model stretching and bending analyses (take
 whichever gives shortest characteristic time), and print. */
    double stretchingTimeStep;
    double bendingTimeStep;

    if (!is_gradient_descent_dynamics_enabled) {
        stretchingTimeStep = time_step_prefactor * sqrt(init_density / shear_modulus) *
                             smallest_size_over_root_tau;
        bendingTimeStep = time_step_prefactor * sqrt(init_density * 6 * (1 - poisson_ratio) / shear_modulus) *
                          (approx_min_init_elem_size / sheet_thickness) * approx_min_init_elem_size / (2 * M_PI);
    } else {
        /* Use second damping factor, since we should only be using gradient descent
        for already dialled-in states. Note the damping factor cancels out in the
        dynamics and has no effect - but only if you are consistent with which one
        you use!*/
        double numDampFac2 = equilibriation_damping * damping_scale;
        stretchingTimeStep =
                time_step_prefactor * (numDampFac2 / shear_modulus) * pow(smallest_size_over_root_tau / (2 * M_PI), 2);
        bendingTimeStep =
                time_step_prefactor *
                (6 * (1 - poisson_ratio) * numDampFac2 / (pow(sheet_thickness, 2) * shear_modulus))
                * pow(smallest_size_over_root_tau / (2 * M_PI), 4);

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
        time_step = stretchingTimeStep;
        logStream << "Time step = " << time_step << ", set based on stretching rather than bending."
                  << std::endl;
    } else {
        time_step = bendingTimeStep;
        logStream << "Time step = " << time_step << ", set based on bending rather than stretching."
                  << std::endl;
    }
    logStream.close();

    /* Fix up DialInStepTime if gradient descent is being used, just to give
reasonable printout and equilibrium check frequencies without extreme settings.*/
    if (is_gradient_descent_dynamics_enabled) {
        dial_in_step_time = 1000 * time_step;
        time_between_equil_checks = time_between_equil_checks_prefactor * dial_in_step_time;
    }
}


void Settings::SetupPrintFrequency(CustomOutStreamClass &logStream) {
    /* Calculate print frequency based on DialInStepTime / TimeStep, tuned by a
dimensionless parameter in the setting file and rounded to the nearest integer.
This rounding should work fine as long as the PrintFrequency isn't ridiculously
huge. To avoid this we require that PrintFrequency is not in [0,1). If
PrintFrequency is very large (more likely), so that InversePrintRate rounds to 0,
we set InversePrintRate to 1 to get a printout after every time step. This can
be a useful thing to do sometimes. If PrintFrequency is set to a negative value
we instead set InversePrintRate to -1, which means no output files are regularly
written at all.
*/
    if (print_frequency < 0) {
        inverse_print_rate = -1;
    } else {
        if (print_frequency < 1) {
            logStream
                    << "Error: To avoid potential divide-by-zero problems we do NOT allow settings.PrintFrequency to lie in [0, 1). This is overkill somewhat, and could be relaxed if need be."
                    << std::endl;
            return;
        } else {
            inverse_print_rate = lround(dial_in_step_time / (time_step * print_frequency));
            if (inverse_print_rate == 0) {
                inverse_print_rate = 1;
            }
        }
    }
}

//void Settings::SetupSmallestElements(CustomOutStreamClass &logStream, std::vector<Triangle> &triangles,
//                                     std::vector<std::vector<double>> &programmed_taus) {
//
//}