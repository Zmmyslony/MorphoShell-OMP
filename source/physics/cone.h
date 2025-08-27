//
// Created by Michał Zmyślony on 26/02/2024.
//

#ifndef MORPHOSHELL_CONE_H
#define MORPHOSHELL_CONE_H

#define _USE_MATH_DEFINES

#include <Eigen/Dense>

#include "../configuration/config_base.h"
#include "rigid_body.h"
#include <cmath>
// For M_PI_2
#include <math.h>

#define FIXED 0
#define LOAD_CONTROLLED 1
#define DISPLACEMENT_CONTROLLED 2

class Cone : public RigidBody {
    double cone_angle = M_PI_2;

    /**
     * Calculates distance from the cone surface.
     * @param pos
     * @return
     */
    [[nodiscard]] double distance(const Eigen::Vector3d &pos) const;

    /**
     * Calculated distance vector from the cone_surface.
     * @param pos
     * @return
     */
    Eigen::Vector3d coneNormal(const Eigen::Vector3d &pos) const;

public:

    explicit Cone(const ConfigBase &config);

    explicit Cone(const fs::path &config_path);

    /**
     * If position is unspecified, sets them depending on the nodes and sets the velocity for displacement controlled.
     * experiment.
     * @param node_pos positions of the nodes
     * @param dial_in_time total duration of the dial-in
     */
    void initialise(const std::vector<Eigen::Vector3d> &node_pos, double dial_in_time);

    /**
     * Applies interaction force to a given node.
     * @param pos position of the node.
     * @param node_force forces corresponding to a given node that will be modified.
     * @param shear_modulus
     * @param thickness
     */
    double addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force, double shear_modulus,
                               double thickness) const;
};


#endif //MORPHOSHELL_CONE_H
