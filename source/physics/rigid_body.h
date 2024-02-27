//
// Created by Michał Zmyślony on 27/02/2024.
//

#ifndef MORPHOSHELL_RIGID_BODY_H
#define MORPHOSHELL_RIGID_BODY_H

#include <Eigen/Dense>

#include "../configuration/config_base.h"

#define FIXED 0
#define LOAD_CONTROLLED 1
#define DISPLACEMENT_CONTROLLED 2

class RigidBody {
protected:
    int type = FIXED;
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

    explicit RigidBody(const fs::path &config_path);

    explicit RigidBody(const ConfigBase &config);

    [[nodiscard]] const Eigen::Vector3d &getPosition() const;

    [[nodiscard]] const Eigen::Vector3d &getVelocity() const;

    [[nodiscard]] double getLoad() const;

    void update(double time_step_size, bool is_dialling_in);

    void setInteractionLoad(double interaction_load);
};


#endif //MORPHOSHELL_RIGID_BODY_H
