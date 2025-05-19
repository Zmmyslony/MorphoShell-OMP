//
// Created by Michał Zmyślony on 19/02/2024.
//

#ifndef MORPHOSHELL_CORE_CONFIG_H
#define MORPHOSHELL_CORE_CONFIG_H

#define _USE_MATH_DEFINES

#include <boost/filesystem/path.hpp>
#include <cfloat>
#include "config_base.h"

namespace fs = boost::filesystem;

/**
 * Core settings controlling the simulation. Need to include thickness and shear_modulus.
 */
class CoreConfig {
    double thickness = DBL_MIN;
    double shear_modulus = DBL_MIN;

    // Number of cores used for parallelization.
    int core_number = 0;
    // Ratio of force to characteristic force required for the state to be considered in equilibrium.
    double equilibrium_force_scale = 1e-9;

    // Ratio of speed to characteristic speed required for the state to be considered in equilibrium.
    double equilibrium_speed_scale = 1e-4;
    double poisson_ratio = 0.5;
    double density = 1000;

    double damping_prefactor = 2;
    double dial_in_damping = 1;
    double equilibration_damping = 1;

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
    // If enabled together with the LCE mode, the elongation and curvature will dynamically change based on height.
    bool is_stimulation_modulated = false;
    bool is_simple_sec_ff_used = true;
    bool is_boundary_clamped = false;
    bool is_initial_positions_perturbed = true;

    bool is_energy_printed = true;
    bool is_stress_printed = false;
    bool is_triangle_area_printed = true;
    bool is_angle_deficit_printed = false;
    bool is_gradient_descent_dynamics = false;
    bool is_seide_deformations = false;
    bool is_force_printed = false;

    // When enabled, current geometry is assumed to be the relaxed one, i.e. uses calculated metric and bend tensors as programmed.
    bool is_ansatz_metric_used = false;

    //
    bool is_first_tensor_skipped = false;

    // Used for finding second fundamental form coefficients. Worse meshes require larger values.
    double patch_matrix_threshold = 10;

    double gentFactor = 1;

    // Boundary conditions - which force should remain zero throughout the simulation
    bool is_x_fixed_bc = true;
    bool is_y_fixed_bc = true;
    bool is_z_fixed_bc = true;

    // Natural units in the simulation (usually mm) - this setting allows you to use other physics in SI units.
    double units = 1e-3;

public:
    CoreConfig();

    explicit CoreConfig(const ConfigBase &config_base);

    explicit CoreConfig(const fs::path &path);

    double getThickness() const;

    double getShearModulus() const;

    int getCoreNumber() const;

    double getEquilibriumForceScale() const;

    double getEquilibriumSpeedScale() const;

    double getPoissonRatio() const;

    double getDensity() const;

    double getDampingPrefactor() const;

    double getDialInDamping() const;

    double getEquilibrationDamping() const;

    double getPrintFrequency() const;

    double getDialInResolution() const;

    double getTimeBetweenEquilibriumChecks() const;

    double getDialInTimePrefactor() const;

    double getTimeStepPrefactor() const;

    bool isLceModeEnabled() const;

    bool isSimpleSecFfUsed() const;

    bool isBoundaryClamped() const;

    bool isEnergyPrinted() const;

    bool isStressPrinted() const;

    bool isForcePrinted() const;

    bool isTriangleAreaPrinted() const;

    bool isAngleDeficitPrinted() const;

    bool isGradientDescentDynamics() const;

    bool isAnsatzMetricUsed() const;

    double getPatchMatrixThreshold() const;

    double getStretchingTimeStep(double size_factor) const;

    double getBendingTimeStep(double size_factor) const;

    double getStretchingTimeStepGradientDescent(double size_factor) const;

    double getBendingTimeStepGradientDescent(double size_factor) const;

    double getStretchingTimeScale(double size_factor) const;

    double getBendingTimeScale(double size_factor) const;

    double getDampingScale(double size_factor) const;

    bool isInitialPositionsPerturbed() const;

    bool isFirstTensorSkipped() const;

    double getGentFactor() const;

    bool isSeideDeformations() const;

    bool isXFixedBc() const;

    bool isYFixedBc() const;

    bool isZFixedBc() const;

    double getUnits() const;

    bool isStimulationModulated() const;

};


#endif //MORPHOSHELL_CORE_CONFIG_H
