//
// Created by Michał Zmyślony on 14/02/2024.
//

#ifndef SHELLOMORPH_CONFIG_BASE_H
#define SHELLOMORPH_CONFIG_BASE_H

#include <iostream>
#include <unordered_map>
#include <boost/filesystem/path.hpp>

namespace fs = boost::filesystem;

/**
 * Abstract class containing pairs of keys and values read from the file. Each pair is read from a single line, ignoring
 * spaces, underscores and capitalisation. Key and value are separated by finding '='.
 */
class ConfigBase {
    std::unordered_map<std::string, std::string> config_options;

public:
    ConfigBase();
    explicit ConfigBase(const fs::path &config_path);

    void print() const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::string &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, int &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, double &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, bool &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::vector<std::string> &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::vector<double> &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::vector<int> &target) const;

    /**
     * Reads option into target address. If option is not found, the value remains unchanged.
     * @param name String name of the option
     * @param target Address into which the option is to be read
     * @return true = success, false = failure
     */
    bool get(const std::string &name, std::vector<bool> &target) const;

    /**
     * Tests whether an option is equal to string
     * @param name
     * @param test_string
     * @return
     */
    bool isEqual(const std::string &name, const std::string &test_string) const;
};

std::string clean_line(std::string);


#endif //SHELLOMORPH_CONFIG_BASE_H
