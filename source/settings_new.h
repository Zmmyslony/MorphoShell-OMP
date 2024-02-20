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

public:
    explicit SettingsNew(const std::vector<fs::path> &config_paths);

    [[nodiscard]] const CoreConfig &getCore() const;

    [[nodiscard]] const GravityConfig &getGravity() const;

    [[nodiscard]] const std::vector<Slide> &getSlides() const;
};


#endif //MORPHOSHELL_SETTINGS_NEW_H
