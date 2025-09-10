//
// Created by Michał Zmyślony on 26/02/2024.
//
#define _USE_MATH_DEFINES

#include "cone.h"
#include <cfloat>
#include <math.h>

Cone::Cone(const ConfigBase &config) : RigidBody(config) {
    config.get("cone_angle", cone_angle);
}

Cone::Cone(const fs::path &config_path) : Cone(ConfigBase(config_path)) {}


Eigen::Vector3d Cone::coneNormal(const Eigen::Vector3d &pos) const {
    Eigen::Vector3d relative_position = pos - position;
    double height = relative_position.dot(normal);
    Eigen::Vector3d height_vector = normal * height;
    Eigen::Vector3d radial_vector = relative_position - height_vector;

    Eigen::Vector3d cone_normal = normal * sin(cone_angle) + radial_vector.normalized() * cos(cone_angle);
    return cone_normal;
}

double Cone::distance(const Eigen::Vector3d &pos) const {
    Eigen::Vector3d relative_position = pos - position;
    Eigen::Vector3d cone_normal = coneNormal(pos);
    double distance = relative_position.dot(cone_normal);
    return distance;
}

void Cone::initialise(const std::vector<Eigen::Vector3d> &node_pos, double dial_in_time) {
    if (!is_origin_provided) {
        position = {0, 0, 0};
        double furthest_distance = DBL_MAX;
        for (auto &pos: node_pos) {
            double current_distance = distance(pos);
            if (current_distance < furthest_distance) { furthest_distance = current_distance; }
        }
        position += furthest_distance * normal / sin(cone_angle);
    }
    if (type == DISPLACEMENT_CONTROLLED) {
        velocity = displacement / dial_in_time;
    }
    if (cone_angle == M_PI_2) {
        std::cerr << "Warning: cone angle set to π / 2 - consider using a slide instead." << std::endl;
        char test = 'n';
        do {
            std::cout << "Continue? [y/N]" << std::endl;
            std::cin >> test;
        } while (!std::cin.fail() && test != 'y' && test != 'n');
        if (test == 'n') {
            std::terminate();
        }
    }
}

double Cone::addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force, double shear_modulus,
                                 double thickness) const {
    Eigen::Vector3d relative_position = pos - position;
    Eigen::Vector3d cone_normal = coneNormal(pos);
    double distance_v = relative_position.dot(cone_normal);
    if (distance_v >= 0) { return 0; }
    double normal_response_force = -distance_v * force_prefactor * shear_modulus * thickness;

    Eigen::Vector3d normal_force = node_force.dot(cone_normal) * cone_normal;
    Eigen::Vector3d tangent_force = node_force - normal_force;

    double friction_force = normal_force.norm() * friction_coefficient;
    double tangent_force_v = tangent_force.norm();

    if (tangent_force_v <= friction_force) {
        node_force -= tangent_force;
    } else {
        node_force -= tangent_force.normalized() * friction_force;
    }
    node_force += normal_response_force * cone_normal;
    return normal_response_force;
}