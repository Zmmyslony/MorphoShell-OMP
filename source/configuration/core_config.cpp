//
// Created by Michał Zmyślony on 19/02/2024.
//

#define _USE_MATH_DEFINES

#include "core_config.h"
#include <thread>
#include <cmath>
#include <math.h>
#include <omp.h>

CoreConfig::CoreConfig() = default;

CoreConfig::CoreConfig(const ConfigBase &config_base) {
    // Required fields.
    if (!config_base.get("thickness", thickness)) {
        throw std::runtime_error("Undefined \"thickness\" in core config.");
    }
    if (!config_base.get("shear_modulus", shear_modulus)) {
        throw std::runtime_error("Undefined \"shear_modulus\" in core config.");
    } else {
        // Correction from SI to a system where mm are natural units.
        shear_modulus *= 1e-3;
    }

    config_base.get("core_number", core_number);
    if (core_number == 0) {
        int processor_count = (int) omp_get_num_procs();
        if (processor_count == 0) {
            core_number = 1;
            std::cerr << "Core number set to zero and unable to detect number of threads. Defaulting to 1."
                      << std::endl;
        } else {
            core_number = processor_count;
            std::cerr << "Core number set to zero. Defaulting to threads count: " << processor_count << "."
                      << std::endl;
        }
    }
    omp_set_num_threads(core_number);

    // Optional fields that have defined defaults.
    config_base.get("equilibrium_force_scale", equilibrium_force_scale);
    config_base.get("equilibrium_speed_scale", equilibrium_speed_scale);
    config_base.get("poisson_ratio", poisson_ratio);
    config_base.get("density", density);
    // SI correction
    density *= 1e-9;

    config_base.get("damping_prefactor", damping_prefactor);
    config_base.get("dial_in_damping", dial_in_damping);
    config_base.get("equilibriation_damping", equilibriation_damping);
    config_base.get("print_frequency", print_frequency);
    config_base.get("dial_in_resolution", dial_in_resolution);
    config_base.get("interval_equilibrium_check", time_between_equilibrium_checks);
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
    config_base.get("is_x_fixed_bc", is_x_fixed_bc);
    config_base.get("is_y_fixed_bc", is_y_fixed_bc);
    config_base.get("is_z_fixed_bc", is_z_fixed_bc);
}

CoreConfig::CoreConfig(const fs::path &path) :
        CoreConfig(ConfigBase(path)) {}

double CoreConfig::getThickness() const {
    return thickness;
}

double CoreConfig::getShearModulus() const {
    return shear_modulus;
}

int CoreConfig::getCoreNumber() const {
    return core_number;
}

double CoreConfig::getEquilibriumForceScale() const {
    return equilibrium_force_scale;
}

double CoreConfig::getEquilibriumSpeedScale() const {
    return equilibrium_speed_scale;
}

double CoreConfig::getPoissonRatio() const {
    return poisson_ratio;
}

double CoreConfig::getDensity() const {
    return density;
}

double CoreConfig::getDampingPrefactor() const {
    return damping_prefactor;
}

double CoreConfig::getDialInDamping() const {
    return dial_in_damping;
}

double CoreConfig::getEquilibriationDamping() const {
    return equilibriation_damping;
}

double CoreConfig::getPrintFrequency() const {
    return print_frequency;
}

double CoreConfig::getDialInResolution() const {
    return dial_in_resolution;
}

double CoreConfig::getTimeBetweenEquilibriumChecks() const {
    return time_between_equilibrium_checks;
}

double CoreConfig::getDialInTimePrefactor() const {
    return dial_in_time_prefactor;
}

double CoreConfig::getTimeStepPrefactor() const {
    return time_step_prefactor;
}

bool CoreConfig::isLceModeEnabled() const {
    return is_lce_mode_enabled;
}

bool CoreConfig::isSimpleSecFfUsed() const {
    return is_simple_sec_ff_used;
}

bool CoreConfig::isBoundaryClamped() const {
    return is_boundary_clamped;
}

bool CoreConfig::isEnergyPrinted() const {
    return is_energy_printed;
}

bool CoreConfig::isTriangleAreaPrinted() const {
    return is_triangle_area_printed;
}

bool CoreConfig::isAngleDeficitPrinted() const {
    return is_angle_deficit_printed;
}

bool CoreConfig::isGradientDescentDynamics() const {
    return is_gradient_descent_dynamics;
}

bool CoreConfig::isAnsatzMetricUsed() const {
    return is_ansatz_metric_used;
}

double CoreConfig::getPatchMatrixThreshold() const {
    return patch_matrix_threshold;
}

double CoreConfig::getDampingScale(double size_factor) const {
    double min_wave_vector = 2 * M_PI / size_factor;
    return 2 * sqrt(density * shear_modulus / (6 * (1 - poisson_ratio))) * thickness * pow(min_wave_vector, 2);

}

double CoreConfig::getStretchingTimeStep(double size_factor) const {
    return time_step_prefactor * sqrt(density / shear_modulus) * size_factor;
}

double CoreConfig::getBendingTimeStep(double size_factor) const {
    return time_step_prefactor * sqrt(density * 6 * (1 - poisson_ratio) / shear_modulus) / (2 * M_PI * thickness) *
           pow(size_factor, 2);
}

double CoreConfig::getStretchingTimeStepGradientDescent(double size_factor) const {
    /* Use equilibrium damping factor, since we should only be using gradient descent
    for already dialled-in states. Note the damping factor cancels out in the
    dynamics and has no effect - but only if you are consistent with which one
    you use!*/
    double numerical_damping = equilibriation_damping * damping_prefactor;
    return time_step_prefactor * (numerical_damping / shear_modulus) * pow(size_factor / (2 * M_PI), 2);

}

double CoreConfig::getBendingTimeStepGradientDescent(double size_factor) const {
    /* Use equilibrium damping factor, since we should only be using gradient descent
    for already dialled-in states. Note the damping factor cancels out in the
    dynamics and has no effect - but only if you are consistent with which one
    you use!*/
    double numerical_damping = equilibriation_damping * damping_prefactor;
    return time_step_prefactor * numerical_damping * (6 * (1 - poisson_ratio) / (pow(thickness, 2) * shear_modulus))
           * pow(size_factor / (2 * M_PI), 4);
}

double CoreConfig::getStretchingTimeScale(double size_factor) const {
    double min_wave_vector = 2 * M_PI / size_factor;
    return sqrt(density / shear_modulus) / min_wave_vector * dial_in_time_prefactor;
}

double CoreConfig::getBendingTimeScale(double size_factor) const {
    double min_wave_vector = 2 * M_PI / size_factor;
    return sqrt(density * 6 * (1 - poisson_ratio) / shear_modulus) / (thickness * pow(min_wave_vector, 2)) *
    dial_in_time_prefactor;
}

bool CoreConfig::isInitialPositionsPerturbed() const {
    return is_initial_positions_perturbed;
}

bool CoreConfig::isFirstTensorSkipped() const {
    return is_first_tensor_skipped;
}

double CoreConfig::getGentFactor() const {
    return gentFactor;
}

bool CoreConfig::isSeideDeformations() const {
    return is_seide_deformations;
}

bool CoreConfig::isXFixedBc() const {
    return is_x_fixed_bc;
}

bool CoreConfig::isYFixedBc() const {
    return is_y_fixed_bc;
}

bool CoreConfig::isZFixedBc() const {
    return is_z_fixed_bc;
}

double CoreConfig::getUnits() const {
    return units;
}

