/*
/////////////////////////////////////////////////////
Copyright (C) 2020, Daniel Duffy, dld34@cam.ac.uk. All rights reserved.
Please cite Daniel Duffy and Dr John Biggins if you use any part of this
code in work that you publish or distribute.

This file is part of Shellmorph.

Shellmorph is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Shellmorph is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////

Function to read in settings file (using libconfig++ library) and put these
these settings in a SettingsStruct which is then returned.
See SettingsStruct.hpp for details of settings.

Libconfig++ is distributed under the Lesser GPL license (2.1 or later), and copyright is held by Mark A Lindner (at least in large part).
*/

#include <libconfig.h++>  //External library, may need downloading
#include <stdexcept>
#include <vector>

#include "readSettingsFile.hpp"
#include "SettingsStruct.hpp"


// TODO  write a method for SettingsStruct which returns the value of parameter from its string name
//       this will allow for clearing up the code and removing a lot of lines

bool tryFillingDoubleParameter(SettingsStruct &settings, const libconfig::Config &config,
                               const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressDouble(valueName));
}

bool tryFillingBoolParameter(SettingsStruct &settings, const libconfig::Config &config,
                             const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressBool(valueName));
}

bool tryFillingIntParameter(SettingsStruct &settings, const libconfig::Config &config,
                            const std::string &valueName) {
    return config.lookupValue(valueName, *settings.getParameterAddressInt(valueName));
}


void checkIfParameterIsDefinedCorrectlyGeneral(SettingsStruct &settings, const libconfig::Config &config,
                                               const std::string &valueName, const std::string &suffix,
                                               bool (*fillingFunction)(SettingsStruct &, const libconfig::Config &,
                                                                       const std::string &)) {
    if (!fillingFunction(settings, config, valueName)) {
        std::string errorMessage =
                "Error reading " + valueName + " setting - it may be missing, or have the wrong format." + suffix;
        throw std::runtime_error(errorMessage);
    }
}

void checkIfDoubleParameterIsDefinedCorrectly(SettingsStruct &settings, const libconfig::Config &config,
                                              const std::string &valueName) {
    std::string suffix = "\n Remember that floating point numbers must have a decimal point in the settings file.";
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, suffix, tryFillingDoubleParameter);
}


void checkIfIntParameterIsDefinedCorrectly(SettingsStruct &settings, const libconfig::Config &config,
                                           const std::string &valueName) {
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, "", tryFillingIntParameter);
}


void checkIfBoolParameterIsDefinedCorrectly(SettingsStruct &settings, const libconfig::Config &config,
                                            const std::string &valueName) {
    checkIfParameterIsDefinedCorrectlyGeneral(settings, config, valueName, "", tryFillingBoolParameter);
}


void readInAllParametersWithAnIntegrityCheck(SettingsStruct &settings, const libconfig::Config &config) {
    std::vector<std::string>
            doubleConfigurationParameters = {"slideStiffnessPrefactor", "slideSpeedPrefactor", "slideFrictionCoeff",
                                             "ThicknessesAboveLowestNodeToClampUpTo", "pSpeedPrefactor",
                                             "constSlideWeightFac", "slideWeightDialSpeedFac", "SpacerHeight",
                                             "SpecifyInitSlideZCoord_upper", "SpecifyInitSlideZCoord_lower",
                                             "totalSlideForceToMuTSqRatioEquilThreshold",
                                             "InitialSlideWeightForCtrldForceInUnitsOfMuTsq", "PrintFrequency",
                                             "Equil_Force_To_CharForce_Ratio_Threshold",
                                             "Equil_Speed_To_SoundSpeed_Ratio_Threshold", "DialInResolution",
                                             "DialInStepTimePrefactor", "TimeBetweenEquilChecksPrefactor",
                                             "ProdForceTime", "ProdStrength", "LoadForceTime", "LoadStrength",
                                             "DampingPrefactor1", "DampingPrefactor2", "GentFactor",
                                             "TimeStepPrefactor", "SheetThickness", "ShearModulus", "PoissonRatio",
                                             "InitDensity", "PatchMatrixDimensionlessConditioningThreshold"};
    std::vector<std::string>
            intConfigurationParameters = {"testTriangle"};
    std::vector<std::string>
            boolConfigurationParameters = {"isDialingDisabled", "GlassCones", "isSeideDeformationsEnabled",
                                           "isControlledForceEnabled", "isLCEModeEnabled", "isEnergyDensitiesPrinted",
                                           "isTriangleAreasPrinted",
                                           "isAngleDeficitsPrinted", "isGradientDescentDynamicsEnabled",
                                           "ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF",
                                           "isDialingFromAnsatzEnabled", "isPerturbationOfInitialPositionsEnabled",
                                           "isBoundaryClamped"};

    for (auto &floatParameterName: doubleConfigurationParameters) {
        checkIfDoubleParameterIsDefinedCorrectly(settings, config, floatParameterName);
    }

    for (auto &intParameterName: intConfigurationParameters) {
        checkIfIntParameterIsDefinedCorrectly(settings, config, intParameterName);
    }

    for (auto &boolParameterName: boolConfigurationParameters) {
        checkIfBoolParameterIsDefinedCorrectly(settings, config, boolParameterName);
    }

    if (settings.ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF &&
        settings.isDialingFromAnsatzEnabled) {
        throw std::runtime_error("Error: The ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF\n \
        and isDialingFromAnsatzEnabled settings give incompatible behaviours, so you're not allowed to have both == true simultaneously. Aborting.");
    }
}


void readSettingsFile(SettingsStruct &settings, const char *settings_file_name) {
    libconfig::Config config;
    config.readFile(settings_file_name);
    readInAllParametersWithAnIntegrityCheck(settings, config);
}