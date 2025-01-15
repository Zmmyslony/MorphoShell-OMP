//
// Created by Michał Zmyślony on 19/02/2024.
//

#define _USE_MATH_DEFINES

#include "settings_new.h"
#include <map>
#include <cstdlib>

#include <cmath>
#include <math.h>


std::map<std::string, int> config_map{
        {"NULL_PLACEHOLDER", -1},
        {"core",             0},
        {"gravity",          1},
        {"slide",            2},
        {"cone",             3}

};


SettingsNew::SettingsNew(const std::vector<fs::path> &config_paths) {
    bool is_core_read = false;
    bool is_gravity_read = false;

    for (auto &path: config_paths) {
        if (path.extension() != ".cfg") { continue; }

        ConfigBase config(path);
        std::string type;
        config.get("type", type);

        switch (config_map[type]) {
            case 0:
                if (is_core_read) { throw std::runtime_error("Repetition of core config."); }
                core = CoreConfig(config);
                is_core_read = true;
                break;
            case 1:
                if (is_gravity_read) { throw std::runtime_error("Repetition of gravity config."); }
                gravity = GravityConfig(config);
                is_gravity_read = true;
                break;
            case 2:
                slides.emplace_back(config);
                break;
            case 3:
                cones.emplace_back(config);
                break;
            default:
                throw std::runtime_error("Unknown configuration file provided.");
        }
    }

    if (!is_core_read) { throw std::runtime_error("Core config missing from input arguments."); }

}

const CoreConfig &SettingsNew::getCore() const {
    return core;
}

const GravityConfig &SettingsNew::getGravity() const {
    return gravity;
}

const std::vector<Slide> &SettingsNew::getSlides() const {
    return slides;
}

std::string SettingsNew::SetupDialInTime(double size_factor) {
    /* Calculate 'dialling in' time and damping coefficient based on toy model
    stretching and bending analyses. The 'Long times' are approximate characteristic
    times for the longest-wavelength modes in the system for stretching and bending.
    The damping coefficient is chosen to approximately critically damp the longest-
    wavelength bending mode.
    */
    std::stringstream log_stream;

    num_damping_prefactor = getCore().getDampingPrefactor() * getCore().getDampingScale(size_factor);
    log_stream << "Numerical damping Coefficient = " << num_damping_prefactor << std::endl;

    double stretching_time_scale = getCore().getStretchingTimeScale(size_factor);
    double bending_time_scale = getCore().getBendingTimeScale(size_factor);

    if (stretching_time_scale > bending_time_scale) {
        duration_phase = stretching_time_scale;
        log_stream << "Phase duration = " << duration_phase
                   << ", set based on stretching rather than bending." << std::endl;
    } else {
        duration_phase = bending_time_scale;
        log_stream << "Phase duration = " << duration_phase
                   << ", set based on bending rather than stretching." << std::endl;
    }
    return log_stream.str();
}


std::string SettingsNew::SetupStepTime(double size_factor) {
    /* Calculate time step based on toy model stretching and bending analyses (take
    whichever gives shortest characteristic time), and print. */
    std::stringstream log_stream;

    double stretchingTimeStep;
    double bendingTimeStep;

    if (!core.isGradientDescentDynamics()) {
        stretchingTimeStep = core.getStretchingTimeStep(size_factor);
        bendingTimeStep = core.getBendingTimeStep(size_factor);
    } else {

        stretchingTimeStep = core.getStretchingTimeStepGradientDescent(size_factor);
        bendingTimeStep = core.getBendingTimeStepGradientDescent(size_factor);

        log_stream
                << "\nUSING GRADIENT DESCENT DYNAMICS." << std::endl
                << "BE WARNED - this feature is only intended for use on nearly converged states supplied as fully dialled-in ansatzes."
                << std::endl
                << "The number of timesteps required for a full simulation with gradient descent is expected to be huge (~10^8)."
                << std::endl
                << "As such, the gradient descent long bending timescale "
                << "is not even used for DialInStepTime,\nas you really shouldn't be doing any dialling in!\nDialInStepTime is instead set such that usual non-gradient-descent settings"
                << " should still give reasonable printout and and equilibrium check frequencies." << std::endl;

    }


    log_stream << "Short stretching and bending timescales: " << stretchingTimeStep << ", " << bendingTimeStep
               << std::endl;

    if (stretchingTimeStep < bendingTimeStep) {
        time_step_size = stretchingTimeStep;
        log_stream << "Time step = " << time_step_size << ", set based on stretching rather than bending."
                   << std::endl;
    } else {
        time_step_size = bendingTimeStep;
        log_stream << "Time step = " << time_step_size << ", set based on bending rather than stretching."
                   << std::endl;
    }

    /* Fix up DialInStepTime if gradient descent is being used, just to give
    reasonable printout and equilibrium check frequencies without extreme settings.*/
    if (core.isGradientDescentDynamics()) {
        duration_phase = 1000 * time_step_size;
        interval_equilibrium_check = core.getDialInTimePrefactor() * duration_phase;
    }

    return log_stream.str();
}


std::string SettingsNew::SetupPrintFrequency() {
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
    std::stringstream log_stream;
    if (core.getPrintFrequency() < 0) {
        interval_print_step = -1;
        return log_stream.str();
    }

    if (core.getPrintFrequency() < 1) {
        log_stream
                << "Error: To avoid potential divide-by-zero problems we do NOT allow settings.PrintFrequency to lie in [0, 1). This is overkill somewhat, and could be relaxed if need be."
                << std::endl;
        return log_stream.str();
    }

    interval_print_step = lround(duration_phase / (time_step_size * core.getPrintFrequency()));
    if (interval_print_step == 0) {
        interval_print_step = 1;
    }
    return log_stream.str();
}

void SettingsNew::SetupCharacteristicSizes(double total_area, double size_factor) {
    interval_equilibrium_check = core.getTimeBetweenEquilibriumChecks() * duration_phase;
    force_scale = core.getShearModulus() * size_factor * core.getThickness();
    stretch_energy_density_scale = core.getShearModulus() * core.getThickness();
    stretch_energy_scale = stretch_energy_density_scale * total_area;
// settings.charBendEnergyDensityScale = settings.ShearModulus * settings.SheetThickness * settings.SheetThickness * settings.SheetThickness / (settings.SampleCharLength * settings.SampleCharLength);
// settings.charBendEnergyScale = settings.charBendEnergyDensityScale * totInitArea;
}

void SettingsNew::useEquilibriumDamping() {
    damping_factor = num_damping_prefactor * core.getEquilibrationDamping();
}

void SettingsNew::useDiallingDamping() {
    damping_factor = num_damping_prefactor * core.getDialInDamping();
}

double SettingsNew::getTimeStepSize() const {
    return time_step_size;
}

double SettingsNew::getDurationPhase() const {
    return duration_phase;
}

double SettingsNew::getDampingFactor() const {
    return damping_factor;
}

double SettingsNew::getForceScale() const {
    return force_scale;
}

int SettingsNew::getStepPrintInterval() const {
    return interval_print_step;
}

double SettingsNew::getStretchEnergyDensityScale() const {
    return stretch_energy_density_scale;
}

double SettingsNew::getStretchEnergyScale() const {
    return stretch_energy_scale;
}

double SettingsNew::getTimeBetweenEquilibriumChecks() const {
    return interval_equilibrium_check;
}

const std::vector<Cone> &SettingsNew::getCones() const {
    return cones;
}

SettingsNew::SettingsNew() = default;



