//
// Created by Michał Zmyślony on 19/02/2024.
//

#include "core_config.h"
#include <thread>


CoreConfig::CoreConfig(const ConfigBase &config_base) {
    // Required fields.
    if (!config_base.get("thickness", thickness)) {
        throw std::runtime_error("Undefined \"thickness\" in core config.");
    }
    if (!config_base.get("shear_modulus", shear_modulus)) {
        throw std::runtime_error("Undefined \"shear_modulus\" in core config.");
    }

    config_base.get("core_number", core_number);
    if (core_number == 0) {
        int processor_count = (int) std::thread::hardware_concurrency();
        if (processor_count == 0) {
            std::cerr << "Core number set to zero and unable to detect number of physical cores. Defaulting to 1."
                      << std::endl;
        } else {
            std::cerr << "Core number set to zero. Defaulting to physical core count: " << processor_count << "."
                      << std::endl;
        }
    }

    // Optional fields that have defined defaults.
    config_base.get("equilibrium_force_scale", equilibrium_force_scale);
    config_base.get("equilibrium_speed_scale", equilibrium_speed_scale);
    config_base.get("poisson_ratio", poisson_ratio);
    config_base.get("density", density);
    config_base.get("damping_prefactor", damping_prefactor);
    config_base.get("dial_in_damping", dial_in_damping);
    config_base.get("equilibriation_damping", equilibriation_damping);
    config_base.get("print_frequency", print_frequency);
    config_base.get("dial_in_resolution", dial_in_resolution);
    config_base.get("time_between_equilibrium_checks", time_between_equilibrium_checks);
    config_base.get("dial_in_time_prefactor", dial_in_time_prefactor);
    config_base.get("time_step_prefactor", time_step_prefactor);
    config_base.get("is_lce_mode_enabled", is_lce_mode_enabled);
    config_base.get("is_simple_sec_ff_used", is_simple_sec_ff_used);
    config_base.get("is_boundary_clamped", is_boundary_clamped);
    config_base.get("is_energy_printed", is_energy_printed);
    config_base.get("is_triangle_area_printed", is_triangle_area_printed);
    config_base.get("is_angle_deficit_printed", is_angle_deficit_printed);
    config_base.get("is_gradient_descent_dynamics", is_gradient_descent_dynamics);
    config_base.get("is_ansatz_metric_used", is_ansatz_metric_used);
    config_base.get("patch_matrix_threshold", patch_matrix_threshold);

}

CoreConfig::CoreConfig(const fs::path &path) :
        CoreConfig(ConfigBase(path)) {}
