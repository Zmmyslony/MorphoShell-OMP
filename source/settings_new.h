//
// Created by Michał Zmyślony on 19/02/2024.
//

#ifndef MORPHOSHELL_SETTINGS_NEW_H
#define MORPHOSHELL_SETTINGS_NEW_H

#include <boost/filesystem/path.hpp>

#include "configuration/core_config.h"
#include "configuration/gravity_config.h"

namespace fs = boost::filesystem;

class SettingsNew {
    CoreConfig core;
    GravityConfig gravity;

public:
    explicit SettingsNew(const std::vector<fs::path> &config_paths);
};


#endif //MORPHOSHELL_SETTINGS_NEW_H
