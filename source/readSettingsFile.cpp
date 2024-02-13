/*
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
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

Function to read in settings file (using libconfig++ library) and put these
these settings in a Settings which is then returned.
See Settings.hpp for details of settings.

Libconfig++ is distributed under the Lesser GPL license (2.1 or later), and copyright is held by Mark A Lindner (at least in large part).
*/

#include <libconfig.h++>  //External library, may need downloading
#include <stdexcept>
#include <vector>

#include "readSettingsFile.hpp"
#include "Settings.hpp"


// TODO  write a method for Settings which returns the value of parameter from its string name
//       this will allow for clearing up the code and removing a lot of lines

bool tryFillingDoubleParameter(Settings &settings, const libconfig::Config &config,
                               const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressDouble(valueName));
}

bool tryFillingBoolParameter(Settings &settings, const libconfig::Config &config,
                             const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressBool(valueName));
}

bool tryFillingIntParameter(Settings &settings, const libconfig::Config &config,
                            const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressInt(valueName));
}


void checkIfParameterIsDefinedCorrectlyGeneral(Settings &settings, const libconfig::Config &config,
                                               const std::string &valueName, const std::string &suffix,
                                               bool (*fillingFunction)(Settings &, const libconfig::Config &,
                                                                       const std::string &)) {
    if (!fillingFunction(settings, config, valueName)) {
        std::string errorMessage =
                "Error reading " + valueName + " setting - it may be missing, or have the wrong format." + suffix;
        throw std::runtime_error(errorMessage);
    }
}

void checkIfDoubleParameterIsDefinedCorrectly(Settings &settings, const libconfig::Config &config,
                                              const std::string &valueName) {
    std::string suffix = "\n Remember that floating point numbers must have a decimal point in the settings file.";
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, suffix, tryFillingDoubleParameter);
}


void checkIfIntParameterIsDefinedCorrectly(Settings &settings, const libconfig::Config &config,
                                           const std::string &valueName) {
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, "", tryFillingIntParameter);
}


void checkIfBoolParameterIsDefinedCorrectly(Settings &settings, const libconfig::Config &config,
                                            const std::string &valueName) {
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, "", tryFillingBoolParameter);
}


void readInAllParametersWithAnIntegrityCheck(Settings &settings, const libconfig::Config &config) {
    std::vector<std::string>
            doubleConfigurationParameters = {"slide_stiffness_prefactor", "slide_speed_prefactor", "slide_friction_coefficient",
                                             "thicknesses_above_lowest_node_to_clamp_up_to", "p_speed_prefactor",
                                             "const_slide_weight_fac", "slide_weight_dial_speed_fac", "spacer_height",
                                             "specify_init_slide_z_coord_upper", "specify_init_slide_z_coord_lower",
                                             "total_slide_force_to_mu_t_sq_ratio_equil_threshold",
                                             "initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq", "print_frequency",
                                             "char_force",
                                             "char_speed", "dial_in_resolution",
                                             "dial_in_step_time_prefactor", "time_between_equil_checks_prefactor",
                                             "prod_force_time", "prod_strength", "load_force_time", "load_strength",
                                             "dial_in_damping", "equilibriation_damping", "gent_factor",
                                             "time_step_prefactor", "sheet_thickness", "shear_modulus", "poisson_ratio",
                                             "init_density", "patch_matrix_dimensionless_conditioning_threshold", "gravity_sign"};
    std::vector<std::string>
            intConfigurationParameters = {"test_triangle"};

    std::vector<std::string>
            boolConfigurationParameters = {"is_dialing_disabled", "glass_cones", "is_controlled_force_enabled",
                                           "is_seide_deformations_enabled", "is_lce_mode_enabled", "is_energy_densities_printed",
                                           "is_triangle_areas_printed",
                                           "is_angle_deficits_printed", "is_gradient_descent_dynamics_enabled",
                                           "for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff",
                                           "is_dialing_from_ansatz_enabled", "is_perturbation_of_initial_positions_enabled",
                                           "is_boundary_clamped"};

    for (auto &floatParameterName: doubleConfigurationParameters) {
        checkIfDoubleParameterIsDefinedCorrectly(settings, config, floatParameterName);
    }

    for (auto &intParameterName: intConfigurationParameters) {
        checkIfIntParameterIsDefinedCorrectly(settings, config, intParameterName);
    }

    for (auto &boolParameterName: boolConfigurationParameters) {
        checkIfBoolParameterIsDefinedCorrectly(settings, config, boolParameterName);
    }

    if (settings.for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff &&
        settings.is_dialing_from_ansatz_enabled) {
        throw std::runtime_error("Error: The ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF\n \
        and isDialingFromAnsatzEnabled settings give incompatible behaviours, so you're not allowed to have both == true simultaneously. Aborting.");
    }
}


void readSettingsFile(Settings &settings, const char *settings_file_name) {
    libconfig::Config config;
    config.readFile(settings_file_name);
    readInAllParametersWithAnIntegrityCheck(settings, config);
    settings.youngs_modulus = 2.0 * settings.shear_modulus * (1.0 + settings.poisson_ratio);
}

void readSettingsFile(Settings &settings, const std::string& settings_filename) {
    readSettingsFile(settings, settings_filename.c_str());
}