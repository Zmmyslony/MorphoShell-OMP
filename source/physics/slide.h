//
// Created by Michał Zmyślony on 20/02/2024.
//

#ifndef MORPHOSHELL_SLIDE_H
#define MORPHOSHELL_SLIDE_H

#include <Eigen/Dense>

#include "../configuration/config_base.h"
#include "rigid_body.h"

class Slide : public RigidBody {
    // Pull-off constant divided by the area.
    double adhesion_constant = 0;
    double adhesion_distance = 0.01;
public:
    explicit Slide(const ConfigBase &config);

    explicit Slide(const fs::path &config_path);

    /**
     * Calculates distance from the slide.
     * @param pos
     * @return
     */
    [[nodiscard]] double distance(const Eigen::Vector3d &pos) const;

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
    double
    addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force, double shear_modulus, double thickness,
                        double node_area) const;

};


#endif //MORPHOSHELL_SLIDE_H
