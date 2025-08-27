//
// Created by Michał Zmyślony on 27/02/2024.
//

#include "rigid_body.h"
#include <cfloat>


void RigidBody::loadNormals(const ConfigBase &config) {
    bool is_normals_provided = config.get("x_normal", normal[0]) &&
                               config.get("y_normal", normal[1]) &&
                               config.get("z_normal", normal[2]);
    if (!is_normals_provided) { throw std::runtime_error("RigidBody normals are missing from config file."); }
    if (normal.norm() == 0) {throw std::runtime_error("Norm of RigidBody normals is equal to zero.");}
    normal.normalize();
}

void RigidBody::loadPositions(const ConfigBase &config) {
    is_origin_provided = config.get("x_position", position[0]) ||
                         config.get("y_position", position[1]) ||
                         config.get("z_position", position[2]);
}

void RigidBody::loadType(const ConfigBase &config) {
    std::string control_type;
    config.get("control_type", control_type);

    if (control_type == clean_line("fixed")) {
        type = FIXED;
    } else if (control_type == clean_line("load_controlled")) {
        type = LOAD_CONTROLLED;
    } else if (control_type == clean_line("displacement_controlled")) {
        type = DISPLACEMENT_CONTROLLED;
        bool is_displacement_loaded = config.get("x_displacement", displacement[0]) ||
                                      config.get("y_displacement", displacement[1]) ||
                                      config.get("z_displacement", displacement[2]);
        if (!is_displacement_loaded) {
            throw std::runtime_error("RigidBody is displacement controlled, yet no displacement was specified.");
        }
    } else {
        throw std::runtime_error(
                "Unknown type of RigidBody. Possible options are 'fixed', 'load_controlled' or "
                "'displacement_controlled'.");
    }
}


void RigidBody::loadInteractionScales(const ConfigBase &config) {
    config.get("load", load);
    config.get("weight", weight);
    config.get("force_prefactor", force_prefactor);
    config.get("friction_coefficient", friction_coefficient);

}

void RigidBody::validate() const {
    if (force_prefactor <= 0) { throw std::runtime_error("'force_prefactor' needs to be a positive number."); }
    if (friction_coefficient < 0) {
        throw std::runtime_error("'friction_coefficient' needs to be a non-negative number.");
    }
    if (type == LOAD_CONTROLLED && weight <= 0) {
        throw std::runtime_error("'weight' needs to be positive for a 'load_controlled' RigidBody.");
    }
}

RigidBody::RigidBody(const ConfigBase &config) {
    loadNormals(config);
    loadPositions(config);
    loadType(config);
    loadInteractionScales(config);
    validate();
}


RigidBody::RigidBody(const fs::path &config_path) : RigidBody(ConfigBase(config_path)) {}

void RigidBody::update(double time_step_size, bool is_dialling_in) {
    switch (type) {
        case FIXED:
            return;
        case LOAD_CONTROLLED:
            velocity += time_step_size * (load - interaction_load) * normal / weight;
            position += time_step_size * velocity;
            break;
        case DISPLACEMENT_CONTROLLED:
            if (is_dialling_in) {
                position += time_step_size * velocity;
            }
            break;
        default:
            throw std::runtime_error("Unknown type of RigidBody control.");
    }
    load = 0;
}

const Eigen::Vector3d &RigidBody::getPosition() const {
    return position;
}

const Eigen::Vector3d &RigidBody::getVelocity() const {
    return velocity;
}

double RigidBody::getLoad() const {
    return load;
}

void RigidBody::setInteractionLoad(double interaction_load) {
    RigidBody::interaction_load = interaction_load;
}
