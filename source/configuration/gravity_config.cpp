//
// Created by Michał Zmyślony on 14/02/2024.
//

#include "gravity_config.h"


GravityConfig::GravityConfig() = default;


GravityConfig::GravityConfig(const fs::path &config_path) :
        GravityConfig(ConfigBase(config_path)) {}


GravityConfig::GravityConfig(const ConfigBase &config_base) {
    config_base.get("is_gravity_enabled", is_gravity_enabled);
    if (is_gravity_enabled) {
        config_base.get("x_gravity_component", x_gravity_component);
        config_base.get("y_gravity_component", y_gravity_component);
        config_base.get("z_gravity_component", z_gravity_component);
    }
    double norm = sqrt(pow(x_gravity_component, 2) +
                       pow(y_gravity_component, 2) +
                       pow(z_gravity_component, 2));

    if (norm != 0) {
        x_gravity_component /= norm;
        y_gravity_component /= norm;
        z_gravity_component /= norm;
    } else {
        is_gravity_enabled = false;
        throw std::runtime_error("When gravity is enabled, at least one of the components need to non-zero.");
    }

}

bool GravityConfig::isGravityEnabled() const {
    return is_gravity_enabled;
}

double GravityConfig::getXGravityComponent() const {
    return x_gravity_component;
}

double GravityConfig::getYGravityComponent() const {
    return y_gravity_component;
}

double GravityConfig::getZGravityComponent() const {
    return z_gravity_component;
}



