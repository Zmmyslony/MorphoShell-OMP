//
// Created by Michał Zmyślony on 20/02/2024.
//

#ifndef MORPHOSHELL_SLIDE_H
#define MORPHOSHELL_SLIDE_H

#include <Eigen/Dense>
#include <map>

#include "../configuration/config_base.h"


std::map<std::string, int> slide_types{
        {"fixed",                   0},
        {"load_controlled",         1},
        {"displacement_controlled", 2}
};

class Slide {
    int slide_type = slide_types["fixed"];
    // If it is not provided, the position will be found as the position of the last node in direction opposite to normal.
    bool is_origin_provided = false;
    Eigen::Vector3d normal = {0, 0, 1};
    Eigen::Vector3d position = {0, 0, 0};
    Eigen::Vector3d velocity = {0, 0, 0};
    Eigen::Vector3d displacement = {0, 0, 0};

    double interaction_load = 0;
    double load = 0;
    double weight = 1;
    double force_prefactor = 1;
    double friction_coefficient = 0;

    void loadNormals(const ConfigBase &config);

    void loadPositions(const ConfigBase &config);

    void loadType(const ConfigBase &config);

    void loadInteractionScales(const ConfigBase &config);

    void validate() const;

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
     * Updates position and velocity of the slide depending on operation mode.
     * @param time_step_size
     * @param is_dialling_in
     */
    void update(double time_step_size, bool is_dialling_in);


    /**
     * Applies interaction force to a given node.
     * @param pos position of the node.
     * @param node_force forces corresponding to a given node that will be modified.
     * @param shear_modulus
     * @param thickness
     */
    void
    addInteractionForce(const Eigen::Vector3d &pos, Eigen::Vector3d &node_force, double shear_modulus,
                        double thickness);

    const Eigen::Vector3d &getPosition() const;

    const Eigen::Vector3d &getVelocity() const;

    double getLoad() const;
};


#endif //MORPHOSHELL_SLIDE_H
