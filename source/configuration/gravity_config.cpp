//
// Created by Michał Zmyślony on 14/02/2024.
//

#include "gravity_config.h"


GravityConfig::GravityConfig() = default;


void GravityConfig::read_fields_from_config() {
    get("is_gravity_enabled", is_gravity_enabled);
    if (is_gravity_enabled) {
        get("x_gravity_component", x_gravity_component);
        get("y_gravity_component", y_gravity_component);
        get("z_gravity_component", z_gravity_component);
    }
}

GravityConfig::GravityConfig(const fs::path &config_path) : ConfigBase(config_path) {
    read_fields_from_config();
}

GravityConfig::GravityConfig(const ConfigBase &config_base) : ConfigBase(config_base) {
    read_fields_from_config();
}



