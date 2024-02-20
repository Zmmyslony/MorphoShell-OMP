//
// Created by Michał Zmyślony on 19/02/2024.
//

#include "settings_new.h"
#include <map>


std::map<std::string, int> config_map{
        {"core", 0},
        {"gravity", 1}
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
                if (is_core_read) { throw std::runtime_error("Repetition of gravity config."); }
                core = CoreConfig(config);
                is_core_read = true;
                break;
            case 1:
                if (is_gravity_read) { throw std::runtime_error("Repetition of gravity config."); }
                gravity = GravityConfig(config);
                is_gravity_read = true;
                break;
            default:
                throw std::runtime_error("Unknown configuration file provided.");
        }
    }

    if (!is_core_read) { throw std::runtime_error("Core config missing from input arguments."); }

}
