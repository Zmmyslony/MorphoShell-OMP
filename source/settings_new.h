//
// Created by Michał Zmyślony on 19/02/2024.
//

#ifndef MORPHOSHELL_SETTINGS_NEW_H
#define MORPHOSHELL_SETTINGS_NEW_H

#include <boost/filesystem/path.hpp>

#include "configuration/core_config.h"
#include "configuration/gravity_config.h"
#include "physics/slide.h"

namespace fs = boost::filesystem;

class SettingsNew {
    CoreConfig core;
    GravityConfig gravity;
    std::vector<Slide> slides;

    double time_step_size = 0;
    double duration_phase = 0;
    double interval_equilibrium_check = 0;
    // Interval in steps between
    int interval_print_step = 0;
    // Constant prefactor calculated from settings and element size
    double num_damping_prefactor = 0;
    double damping_factor = 0;
    double force_scale = 0;
    double stretch_energy_density_scale = 0;
    double stretch_energy_scale = 0;

public:
    SettingsNew();

    explicit SettingsNew(const std::vector<fs::path> &config_paths);

    [[nodiscard]] const CoreConfig &getCore() const;

    [[nodiscard]] const GravityConfig &getGravity() const;

    [[nodiscard]] const std::vector<Slide> &getSlides() const;

    std::string SetupPrintFrequency();

    std::string SetupStepTime(double size_factor);

    std::string SetupDialInTime(double size_factor);

    void SetupCharacteristicSizes(double total_area, double size_factor);
    void useEquilibriumDamping();
    void useDiallingDamping();

    double getTimeStepSize() const;

    double getDurationPhase() const;

    double getDampingFactor() const;

    double getForceScale() const;

    int getStepPrintInterval() const;

    double getStretchEnergyDensityScale() const;

    double getStretchEnergyScale() const;

    double getTimeBetweenEquilibriumChecks() const;
};


#endif //MORPHOSHELL_SETTINGS_NEW_H
