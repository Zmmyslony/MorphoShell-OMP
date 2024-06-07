//
// Created by Michał Zmyślony on 14/02/2024.
//

#include "gravity_config.h"

#include <cmath>


GravityConfig::GravityConfig() = default;


GravityConfig::GravityConfig(const fs::path &config_path) :
        GravityConfig(ConfigBase(config_path)) {}


GravityConfig::GravityConfig(const ConfigBase &config_base) {
    config_base.get("is_gravity_enabled", is_gravity_enabled);
    if (is_gravity_enabled) {
        config_base.get("x_normal", x_normal);
        config_base.get("y_normal", y_normal);
        config_base.get("z_normal", z_normal);
    }
    double norm = sqrt(pow(x_normal, 2) +
                       pow(y_normal, 2) +
                       pow(z_normal, 2));
    config_base.get("gravity_magnitude", gravity_magnitude);
    if (norm != 0) {
        x_normal /= norm;
        y_normal /= norm;
        z_normal /= norm;
    } else {
        is_gravity_enabled = false;
        throw std::runtime_error("When gravity is enabled, at least one of the components need to non-zero.");
    }

}

bool GravityConfig::isGravityEnabled() const {
    return is_gravity_enabled;
}

double GravityConfig::getXGravityComponent() const {
    return x_normal;
}

double GravityConfig::getYGravityComponent() const {
    return y_normal;
}

double GravityConfig::getZGravityComponent() const {
    return z_normal;
}



