//
// Created by Michał Zmyślony on 20/02/2024.
//

#include "slide.h"

double Slide::distance(const Eigen::Vector3d &pos) const {
    return (pos - position).dot(normal);
}


void Slide::loadNormals(const ConfigBase &config) {
    bool is_normals_provided = config.get("x_normal", normal[0]) &&
                               config.get("y_normal", normal[1]) &&
                               config.get("z_normal", normal[2]);
    if (!is_normals_provided) { throw std::runtime_error("Slide normals are missing from config file."); }
    normal.normalize();
}

void Slide::loadPositions(const ConfigBase &config) {
    is_origin_provided = config.get("x_position", position[0]) ||
                         config.get("y_position", position[1]) ||
                         config.get("z_position", position[2]);
}

void Slide::loadType(const ConfigBase &config) {
    std::string control_type;
    config.get("control_type", control_type);

    if (control_type == clean_line("fixed")) {
        slide_type = slide_types["fixed"];
    } else if (control_type == clean_line("load_controlled")) {
        slide_type = slide_types["load_controlled"];
    } else if (control_type == clean_line("displacement_controlled")) {
        slide_type = slide_types["displacement_controlled"];
        bool is_displacement_loaded = config.get("x_displacement", displacement[0]) ||
                                      config.get("y_displacement", displacement[1]) ||
                                      config.get("z_displacement", displacement[2]);
        if (!is_displacement_loaded) {
            throw std::runtime_error("Slide is displacement controlled, yet no displacement was specified.");
        }
    } else {
        throw std::runtime_error(
                "Unknown type of slide. Possible options are 'fixed', 'load_controlled' or "
                "'displacement_controlled'.");
    }
}


void Slide::loadInteractionScales(const ConfigBase &config) {
    config.get("load", load);
    config.get("weight", weight);
    config.get("force_prefactor", force_prefactor);
    config.get("friction_coefficient", friction_coefficient);

}

void Slide::validate() const {
    if (force_prefactor <= 0) { throw std::runtime_error("'force_prefactor' needs to be a positive number."); }
    if (friction_coefficient < 0) {
        throw std::runtime_error("'friction_coefficient' needs to be a non-negative number.");
    }
    if (slide_type == slide_types["load_controlled"] && weight <= 0) {
        throw std::runtime_error("'weight' needs to be positive for a 'load_controlled' slide.");
    }
}

Slide::Slide(const ConfigBase &config) {
    loadNormals(config);
    loadPositions(config);
    loadType(config);
    loadInteractionScales(config);
    validate();
}


Slide::Slide(const fs::path &config_path) : Slide(ConfigBase(config_path)) {}

void Slide::initialise(const std::vector<Eigen::Vector3d> &node_pos, double dial_in_time) {
    if (!is_origin_provided) {
        position = {0, 0, 0};
        double furthest_distance = DBL_MAX;
        for (auto &pos: node_pos) {
            double current_distance = distance(pos);
            if (current_distance < furthest_distance) { furthest_distance = current_distance; }
        }
        position += furthest_distance * normal;
    }
    if (slide_type == slide_types["displacement_controlled"]) {
        velocity = displacement / dial_in_time;
    }
}

void Slide::update(double time_step_size, bool is_dialling_in) {
    switch (slide_type) {
        case 0:
            return;
        case 1:
            velocity += time_step_size * (load - interaction_load) * normal / weight;
            position += time_step_size * velocity;
            break;
        case 2:
            if (is_dialling_in) {
                position += time_step_size * velocity;
            }
            break;
        default:
            throw std::runtime_error("Unknown type of slide control.");
    }
    load = 0;
}

void Slide::addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force,
                                double shear_modulus, double thickness) {
    Eigen::Vector3d interaction = {0, 0, 0};
    double distance_v = distance(pos);
    if (distance_v >= 0) { return; }
    double normal_force = distance_v * force_prefactor * shear_modulus * thickness;
    double friction_force = normal_force * friction_coefficient;

    Eigen::Vector3d tangent_force = node_force - node_force.dot(normal) * normal;
    double tangent_force_v = tangent_force.norm();
    if (tangent_force_v <= friction_force) {
        node_force -= tangent_force;
    } else {
        node_force -=  tangent_force.normalized() * friction_force;
    }
    node_force -= normal_force * normal;
    interaction_load += normal_force;
}

const Eigen::Vector3d &Slide::getPosition() const {
    return position;
}

const Eigen::Vector3d &Slide::getVelocity() const {
    return velocity;
}

double Slide::getLoad() const {
    return load;
}

