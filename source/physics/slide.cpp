//
// Created by Michał Zmyślony on 20/02/2024.
//

#include "slide.h"
#include <cfloat>

Slide::Slide(const ConfigBase &config) : RigidBody(config) {}

Slide::Slide(const fs::path &config_path) : RigidBody(config_path) {}

double Slide::distance(const Eigen::Vector3d &pos) const {
    return (pos - position).dot(normal);
}

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
    if (type == DISPLACEMENT_CONTROLLED) {
        velocity = displacement / dial_in_time;
    }
}

double Slide::addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force, double shear_modulus,
                                  double thickness, double node_area) const {
    double distance_v = distance(pos);
    double normal_force_value = 0;
    if (distance_v < 0) {
        normal_force_value -= distance_v * force_prefactor * shear_modulus * thickness;

        Eigen::Vector3d normal_force = node_force.dot(normal) * normal;
        Eigen::Vector3d tangent_force = node_force - normal_force;
        double friction_force = normal_force.norm() * friction_coefficient;
        double tangent_force_v = tangent_force.norm();
        if (tangent_force_v <= friction_force) {
            node_force -= tangent_force;
        } else {
            node_force -= tangent_force.normalized() * friction_force;
        }
    }
    if (adhesion_constant > 0) {
        normal_force_value -= node_area * adhesion_constant * exp(-pow(distance_v / adhesion_distance, 2));
    }

    node_force += normal_force_value * normal;
    return normal_force_value;
}

