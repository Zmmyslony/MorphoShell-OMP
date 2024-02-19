//
// Created by Michał Zmyślony on 14/02/2024.
//

#ifndef SHELLOMORPH_GRAVITY_CONFIG_H
#define SHELLOMORPH_GRAVITY_CONFIG_H


#include "config_base.h"

class GravityConfig {
    bool is_gravity_enabled = false;
    double x_gravity_component = 0;
    double y_gravity_component = 0;
    double z_gravity_component = 0;

public:
    GravityConfig();

    explicit GravityConfig(const fs::path &config_path);

    explicit GravityConfig(const ConfigBase &config_base);
};


#endif //SHELLOMORPH_GRAVITY_CONFIG_H
