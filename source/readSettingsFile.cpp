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
//
//void readSettingsFile(SettingsStruct &settings, const char *settings_file_name) {
//
//    libconfig::Config config;
//    config.readFile(settings_file_name);
//
//
//    // FOR CONE SQUASHING/BUCKLING BETWEEN TWO SLIDES.
//    if (config.lookupValue("slideStiffnessPrefactor", settings.slideStiffnessPrefactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading slideStiffnessPrefactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("slideSpeedPrefactor", settings.slideSpeedPrefactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading slideSpeedPrefactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//    if (config.lookupValue("slideFrictionCoeff", settings.slideFrictionCoeff)) { ; }
//    else {
//        throw std::runtime_error("Error reading slideFrictionCoeff setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//    if (config.lookupValue("ThicknessesAboveLowestNodeToClampUpTo",
//                           settings.ThicknessesAboveLowestNodeToClampUpTo)) { ; }
//    else {
//        throw std::runtime_error("Error reading ThicknessesAboveLowestNodeToClampUpTo setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//    if (config.lookupValue("isDialingDisabled", settings.isDialingDisabled)) { ; }
//    else {
//        throw std::runtime_error("Error reading isDialingDisabled setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("GlassCones", settings.GlassCones)) { ; }
//    else {
//        throw std::runtime_error("Error reading GlassCones setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("isSeideDeformationsEnabled", settings.isSeideDeformationsEnabled)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isSeideDeformationsEnabled setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("pSpeedPrefactor", settings.pSpeedPrefactor)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading pSpeedPrefactor setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("testTriangle", settings.testTriangle)) { ; }
//    else {
//        throw std::runtime_error("Error reading testTriangle setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("isControlledForceEnabled", settings.isControlledForceEnabled)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isControlledForceEnabled setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("constSlideWeightFac", settings.constSlideWeightFac)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading constSlideWeightFac setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("slideWeightDialSpeedFac", settings.slideWeightDialSpeedFac)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading slideWeightDialSpeedFac setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("SpacerHeight", settings.SpacerHeight)) { ; }
//    else {
//        throw std::runtime_error("Error reading SpacerHeight setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("SpecifyInitSlideZCoord_upper", settings.SpecifyInitSlideZCoord_upper)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading SpecifyInitSlideZCoord_upper setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("SpecifyInitSlideZCoord_lower", settings.SpecifyInitSlideZCoord_lower)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading SpecifyInitSlideZCoord_lower setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("totalSlideForceToMuTSqRatioEquilThreshold",
//                           settings.totalSlideForceToMuTSqRatioEquilThreshold)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading totalSlideForceToMuTSqRatioEquilThreshold setting - it may be missing, or have the wrong format.");
//    }
//    if (config.lookupValue("InitialSlideWeightForCtrldForceInUnitsOfMuTsq",
//                           settings.InitialSlideWeightForCtrldForceInUnitsOfMuTsq)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading InitialSlideWeightForCtrldForceInUnitsOfMuTsq setting - it may be missing, or have the wrong format.");
//    }
//
//
//    if (config.lookupValue("PrintFrequency", settings.PrintFrequency)) { ; }
//    else {
//        throw std::runtime_error("Error reading PrintFrequency setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("isLCEModeEnabled", settings.isLCEModeEnabled)) { ; }
//    else {
//        throw std::runtime_error("Error reading isLCEModeEnabled setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("isEnergyDensitiesPrinted", settings.isEnergyDensitiesPrinted)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isEnergyDensitiesPrinted setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("isTriangleAreasPrinted", settings.isTriangleAreasPrinted)) { ; }
//    else {
//        throw std::runtime_error("Error reading isTriangleAreasPrinted setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("isAngleDeficitsPrinted", settings.isAngleDeficitsPrinted)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isAngleDeficitsPrinted setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("Equil_Force_To_CharForce_Ratio_Threshold",
//                           settings.Equil_Force_To_CharForce_Ratio_Threshold)) { ; }
//    else {
//        throw std::runtime_error("Error reading Equil_Force_To_CharForce_Ratio_Threshold setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("Equil_Speed_To_SoundSpeed_Ratio_Threshold",
//                           settings.Equil_Speed_To_SoundSpeed_Ratio_Threshold)) { ; }
//    else {
//        throw std::runtime_error("Error reading Equil_Speed_To_SoundSpeed_Ratio_Threshold setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("isGradientDescentDynamicsEnabled", settings.isGradientDescentDynamicsEnabled)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isGradientDescentDynamicsEnabled setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF",
//                           settings.ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("isDialingFromAnsatzEnabled", settings.isDialingFromAnsatzEnabled)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isDialingFromAnsatzEnabled setting - it may be missing, or have the wrong format.");
//    }
//
//    if (settings.ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF == true &&
//        settings.isDialingFromAnsatzEnabled == true) {
//        throw std::runtime_error("Error: The ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF\n \
//        and isDialingFromAnsatzEnabled settings give incompatible behaviours, so you're not allowed to have both == true simultaneously. Aborting.");
//    }
//
//    if (config.lookupValue("DialInResolution", settings.DialInResolution)) { ; }
//    else {
//        throw std::runtime_error("Error reading DialInResolution setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("DialInStepTimePrefactor", settings.DialInStepTimePrefactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading DialInStepTimePrefactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("TimeBetweenEquilChecksPrefactor", settings.TimeBetweenEquilChecksPrefactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading TimeBetweenEquilChecksPrefactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("ProdForceTime", settings.ProdForceTime)) { ; }
//    else {
//        throw std::runtime_error("Error reading ProdForceTime setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("ProdStrength", settings.ProdStrength)) { ; }
//    else {
//        throw std::runtime_error("Error reading ProdStrength setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("LoadForceTime", settings.LoadForceTime)) { ; }
//    else {
//        throw std::runtime_error("Error reading LoadForceTime setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("LoadStrength", settings.LoadStrength)) { ; }
//    else {
//        throw std::runtime_error("Error reading LoadStrength setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("isPerturbationOfInitialPositionsEnabled",
//                           settings.isPerturbationOfInitialPositionsEnabled)) { ; }
//    else {
//        throw std::runtime_error(
//                "Error reading isPerturbationOfInitialPositionsEnabled setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("isBoundaryClamped", settings.isBoundaryClamped)) { ; }
//    else {
//        throw std::runtime_error("Error reading isBoundaryClamped setting - it may be missing, or have the wrong format.");
//    }
//
//    if (config.lookupValue("DampingPrefactor1", settings.DampingPrefactor1)) { ; }
//    else {
//        throw std::runtime_error("Error reading DampingPrefactor1 setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("DampingPrefactor2", settings.DampingPrefactor2)) { ; }
//    else {
//        throw std::runtime_error("Error reading DampingPrefactor2 setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("GentFactor", settings.GentFactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading GentFactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("TimeStepPrefactor", settings.TimeStepPrefactor)) { ; }
//    else {
//        throw std::runtime_error("Error reading TimeStepPrefactor setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("SheetThickness", settings.SheetThickness)) { ; }
//    else {
//        throw std::runtime_error("Error reading SheetThickness setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("ShearModulus", settings.ShearModulus)) { ; }
//    else {
//        throw std::runtime_error("Error reading ShearModulus setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("PoissonRatio", settings.PoissonRatio)) { ; }
//    else {
//        throw std::runtime_error("Error reading PoissonRatio setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("InitDensity", settings.InitDensity)) { ; }
//    else {
//        throw std::runtime_error("Error reading InitDensity setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//
//    if (config.lookupValue("PatchMatrixDimensionlessConditioningThreshold",
//                           settings.PatchMatrixDimensionlessConditioningThreshold)) { ; }
//    else {
//        throw std::runtime_error("Error reading PatchMatrixDimensionlessConditioningThreshold setting - it may be missing, or have the wrong format.\n \
//        Remember floating point numbers must have a decimal point in the settings file.");
//    }
//}
