//
// Created by Michał Zmyślony on 19/02/2024.
//

#include "slide_configuration.h"

SlideConfiguration::SlideConfiguration() = default;

SlideConfiguration::SlideConfiguration(const ConfigBase &config) {
    std::vector<std::string> slide_types_s;
    config.get("slide_types", slide_types_s);

    for (const auto &slide_type_s: slide_types_s) {
        if (slide_type_s == clean_line("fixed")) {
            slide_types.push_back(slide_type["fixed"]);
        } else if (slide_type_s == clean_line("load_controlled")) {
            slide_types.push_back(slide_type["load_controlled"]);
        } else if (slide_type_s == clean_line("displacement_controlled")) {
            slide_types.push_back(slide_type["displacement_controlled"]);
        } else {
            throw std::runtime_error(
                    "Unknown type of slide. Possible options are 'fixed', 'load_controlled' or "
                    "'displacement_controlled'");
        }
    }

    config.get("x_slide_normal", x_slide_normal);
    config.get("y_slide_normal", y_slide_normal);
    config.get("z_slide_normal", z_slide_normal);

    config.get("is_origin_provided", is_origin_provided);
    config.get("x_slide_initial", x_slide_initial);
    config.get("y_slide_initial", y_slide_initial);
    config.get("z_slide_initial", z_slide_initial);

    int slides_count = slide_types.size();
    if (x_slide_normal.size() != slides_count ||
        y_slide_normal.size() != slides_count ||
        z_slide_normal.size() != slides_count ||
        is_origin_provided.size() != slides_count ||
        x_slide_initial.size() != slides_count ||
        y_slide_initial.size() != slides_count ||
        z_slide_initial.size() != slides_count) {

        throw std::runtime_error("Uneven number of slide properties. Ensure that normals and initial "
                                 "coordinates are provided.");
    }
}

SlideConfiguration::SlideConfiguration(const fs::path &path) :
        SlideConfiguration(ConfigBase(path)) {}
