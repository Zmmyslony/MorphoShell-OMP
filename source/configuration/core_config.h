//
// Created by Michał Zmyślony on 19/02/2024.
//

#ifndef MORPHOSHELL_CORE_CONFIG_H
#define MORPHOSHELL_CORE_CONFIG_H

#include "config_base.h"
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;


class CoreConfig {
    double thickness{};
    double shear_modulus{};

    // Number of cores used for parallelization.
    int core_number = 0;
    // Ratio of force to characteristic force required for the state to be considered in equilibrium.
    double equilibrium_force_scale = 1e-9;

    // Ratio of speed to characteristic speed required for the state to be considered in equilibrium.
    double equilibrium_speed_scale = 1e-4;
    double poisson_ratio = 0.5;
    double density = 1000;

    double damping_prefactor = 2;
    double dial_in_damping = 2;
    double equilibriation_damping = 10;

    // Number of vtks saved and details displayed per each dial-in.
    double print_frequency = 5;

    // Non-
    double dial_in_resolution = 1.1;
    // Ratio of dial-in time to the time between equilibrium checks.
    double time_between_equilibrium_checks = 0.25;
    // Multiplies
    double dial_in_time_prefactor = 0.5;
    // Multiplies the time in each time step, i.e. larger value means more time per step.
    double time_step_prefactor = 0.5;

    // LCE mode - using director, default - using programmed metric and second fundamental form.
    bool is_lce_mode_enabled = false;
    bool is_simple_sec_ff_used = true;
    bool is_boundary_clamped = false;

    bool is_energy_printed = true;
    bool is_triangle_area_printed = true;
    bool is_angle_deficit_printed = false;
    bool is_gradient_descent_dynamics = false;

    // When enabled, current geometry is used as the preferred one in the beginning.
    bool is_ansatz_metric_used = false;

    // Used for finding second fundamental form coefficients. Worse meshes require larger values.
    double patch_matrix_threshold = 10;

public:
    explicit CoreConfig(const ConfigBase &config_base);
    explicit CoreConfig(const fs::path& path);



};


#endif //MORPHOSHELL_CORE_CONFIG_H
