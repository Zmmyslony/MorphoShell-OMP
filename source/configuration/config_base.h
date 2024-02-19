//
// Created by Michał Zmyślony on 14/02/2024.
//

#ifndef SHELLOMORPH_CONFIG_BASE_H
#define SHELLOMORPH_CONFIG_BASE_H

#include <iostream>
#include <unordered_set>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

class ConfigBase {
    std::vector<std::pair<std::string, std::string>> config_options;

public:
    ConfigBase();
    explicit ConfigBase(const fs::path &config_path);

    /**
     * Reads option into target address. If option is not found, the value is unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::string &target) const;
    bool get(const std::string &name, int &target) const;
    bool get(const std::string &name, double &target) const;
    bool get(const std::string &name, bool &target) const;
    bool get(const std::string &name, std::vector<std::string> &target) const;
    bool get(const std::string &name, std::vector<double> &target) const;
    bool get(const std::string &name, std::vector<int> &target) const;
    bool get(const std::string &name, std::vector<bool> &target) const;

};


#endif //SHELLOMORPH_CONFIG_BASE_H
