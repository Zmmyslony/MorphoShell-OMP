//
// Created by Michał Zmyślony on 19/02/2024.
//

#ifndef MORPHOSHELL_SLIDE_CONFIGURATION_H
#define MORPHOSHELL_SLIDE_CONFIGURATION_H

#include "config_base.h"

#include <map>

const std::map<std::string, int> slide_type {
        {"fixed", 0},
        {"load_controlled", 1},
        {"displacement_controlled", 2}
};

/**
 * Configuration of the slide squashing. Slides are characterised by their initial position, movement direction and
 * type.
 */
class SlideConfiguration {
    std::vector<int> slide_types;
    std::vector<double> x_slide_normal;
    std::vector<double> y_slide_normal;
    std::vector<double> z_slide_normal;

    std::vector<bool> is_origin_provided;
    std::vector<double> x_slide_initial;
    std::vector<double> y_slide_initial;
    std::vector<double> z_slide_initial;

public:
    SlideConfiguration();
    explicit SlideConfiguration(const ConfigBase &config);
    explicit SlideConfiguration(const fs::path &path);

};


#endif //MORPHOSHELL_SLIDE_CONFIGURATION_H
