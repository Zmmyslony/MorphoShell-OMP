//
// Created by Michał Zmyślony on 14/02/2024.
//

#include "config_base.h"
#include <algorithm>
#include <cctype>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <istream>


ConfigBase::ConfigBase() = default;

std::string clean_line(std::string line) {
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
    std::transform(line.begin(), line.end(), line.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return line;
}


bool is_comment(const std::string &option) {
    if (option.empty() ||
        option[0] == '#' ||
        option.size() > 2 && option[0] == '\\' && option[1] == '\\'
            ) {
        return true;
    }
    return false;
}

ConfigBase::ConfigBase(const fs::path &config_path) {
    std::string line;
    std::ifstream file(config_path.string());

    while (std::getline(file, line)) {
        std::string cleaned_line = clean_line(line);
        std::string element;
        std::stringstream line_stream(line);

        std::vector<std::string> row;

        while (std::getline(line_stream, element, '=')) {
            row.push_back(element);
        }

        if (is_comment(row[0])) {
            continue;
        }
        if (row.size() > 2) {
            std::cerr << "Too many entries for " << row[0] << " in " << config_path.string() <<
                                       ".\nUsing only the first read value." << std::endl;
        } else if (row.size() < 2) {
            std::cerr << "Too few entries for " << row[0] << " in " << config_path.string() << "." << std::endl;
        }

        std::pair<std::string, std::string> config_pair = {row[0], row[1]};
        config_options.push_back(config_pair);
    }
}

bool ConfigBase::get(const std::string &name, std::string &target) const {
    std::string normalised_name = clean_line(name);
    for (auto &pair : config_options) {
        if (pair.first == normalised_name) {
            target = pair.second;
            return true;
        }
    }
    target = "NULL_PLACEHOLDER";
    return false;
}

bool ConfigBase::get(const std::string &name, double &target) const {
    std::string read_value;
    bool is_in_config = get(name, read_value);
    if (is_in_config) {
        target = std::stod(read_value);
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, int &target) const {
    std::string read_value;
    bool is_in_config = get(name, read_value);
    if (is_in_config) {
        target = std::stoi(read_value);
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, bool &target) const {
    std::string read_value;
    bool is_in_config = get(name, read_value);
    if (is_in_config) {
        std::istringstream(read_value) >> std::boolalpha >> target;
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, std::vector<std::string> &target) const {
    std::string read_value;
    bool is_in_config = get(name, read_value);

    std::vector<std::string> string_vector;
    if (is_in_config) {
        std::string element;
        std::stringstream stream (read_value);
        while (std::getline(stream, element, ',')) {
            string_vector.push_back(element);
        }
        target = string_vector;
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, std::vector<double> &target) const {
    std::vector<std::string> read_vector;
    bool is_in_config = get(name, read_vector);

    std::vector<double> result;
    if (is_in_config) {
        for (auto &element : read_vector) {
            result.push_back(std::stod(element));
        }

        target = result;
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, std::vector<int> &target) const {
    std::vector<std::string> read_vector;
    bool is_in_config = get(name, read_vector);

    std::vector<int> result;
    if (is_in_config) {
        for (auto &element : read_vector) {
            result.push_back(std::stoi(element));
        }

        target = result;
        return true;
    }
    return false;
}

bool ConfigBase::get(const std::string &name, std::vector<bool> &target) const {
    std::vector<std::string> read_vector;
    bool is_in_config = get(name, read_vector);

    std::vector<bool> result;
    if (is_in_config) {
        for (auto &element : read_vector) {
            bool element_b;
            std::istringstream(element) >> std::boolalpha >> element_b;
            result.push_back(element_b);
        }

        target = result;
        return true;
    }
    return false;
}


bool ConfigBase::isEqual(const std::string &name, const std::string &test_string) const {
    std::string target;
    get(name, target);
    std::string cleaned_string = clean_line(test_string);
    return target == cleaned_string;
}



