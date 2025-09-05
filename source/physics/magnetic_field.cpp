//
// Created by Michał Zmyślony on 29/08/2025.
//

#include "magnetic_field.h"

MagneticField::MagneticField(const ConfigBase &config) {
    std::vector<double> magnetic_field_vec;
     config.get("magnetic_field", magnetic_field_vec);
    magnetic_field[0] = magnetic_field_vec[0];
    magnetic_field[1] = magnetic_field_vec[1];
    magnetic_field[2] = magnetic_field_vec[2];
}