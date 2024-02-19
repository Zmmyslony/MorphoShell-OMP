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
}



