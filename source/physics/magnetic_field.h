//
// Created by Michał Zmyślony on 29/08/2025.
//

#ifndef MORPHOSHELL_MAGNETIC_FIELD_H
#define MORPHOSHELL_MAGNETIC_FIELD_H

#include <Eigen/Dense>
#include "../configuration/config_base.h"

class MagneticField {

public:
    Eigen::Vector3d magnetic_field;
    explicit MagneticField(const ConfigBase &config);
};


#endif //MORPHOSHELL_MAGNETIC_FIELD_H