//
// Created by Michał Zmyślony on 22/03/2022.
//

#include "Settings.hpp"

#include <map>

double *Settings::getParameterAddressDouble(const std::string &parameterName) {
    std::map<std::string, double *> doubleMap{
            {"initSlideZCoord_lower",                         &initSlideZCoord_lower},
            {"initSlideZCoord_upper",                         &initSlideZCoord_upper},
            {"currSlideZCoord_upper",                         &currSlideZCoord_upper},
            {"slideStiffnessPrefactor",                       &slideStiffnessPrefactor},
            {"slideSpeedPrefactor",                           &slideSpeedPrefactor},
            {"upperSlideDisplacement",                        &upperSlideDisplacement},
            {"slideFrictionCoeff",                            &slideFrictionCoeff},
            {"ThicknessesAboveLowestNodeToClampUpTo",         &ThicknessesAboveLowestNodeToClampUpTo},
            {"BendingLongTime",                               &BendingLongTime},
            {"ConeAngle",                                     &ConeAngle},
            {"lambda",                                        &lambda},
            {"slideWeightDialSpeedFac",                       &slideWeightDialSpeedFac},
            {"totalSlideForceToMuTSqRatioEquilThreshold",     &totalSlideForceToMuTSqRatioEquilThreshold},
            {"InitialSlideWeightForCtrldForceInUnitsOfMuTsq", &InitialSlideWeightForCtrldForceInUnitsOfMuTsq},
            {"slideDampingParam",                             &slideDampingParam},
            {"upperSlideWeight",                              &upperSlideWeight},
            {"upperSlideVel",                                 &upperSlideVel},
            {"upperTotSlideForce",                            &upperTotSlideForce},
            {"constSlideWeightFac",                           &constSlideWeightFac},
            {"SpacerHeight",                                  &SpacerHeight},
            {"p",                                             &p},
            {"pSpeedPrefactor",                               &pSpeedPrefactor},
            {"SpecifyInitSlideZCoord_upper",                  &SpecifyInitSlideZCoord_upper},
            {"SpecifyInitSlideZCoord_lower",                  &SpecifyInitSlideZCoord_lower},
            {"PrintFrequency",                                &PrintFrequency},
            {"NumDampFactor",                                 &NumDampFactor},
            {"DampingPrefactor1",                             &DampingPrefactor1},
            {"DampingPrefactor2",                             &DampingPrefactor2},
            {"GentFactor",                                    &GentFactor},
            {"TimeStepPrefactor",                             &TimeStepPrefactor},
            {"TimeStep",                                      &TimeStep},
            {"Equil_Force_To_CharForce_Ratio_Threshold",      &Equil_Force_To_CharForce_Ratio_Threshold},
            {"Equil_Speed_To_SoundSpeed_Ratio_Threshold",     &Equil_Speed_To_SoundSpeed_Ratio_Threshold},
            {"DialInResolution",                              &DialInResolution},
            {"DialInStepTime",                                &DialInStepTime},
            {"TimeBetweenEquilChecksPrefactor",               &TimeBetweenEquilChecksPrefactor},
            {"TimeBetweenEquilChecks",                        &TimeBetweenEquilChecks},
            {"DialInStepTimePrefactor",                       &DialInStepTimePrefactor},
            {"ProdForceTime",                                 &ProdForceTime},
            {"ProdStrength",                                  &ProdStrength},
            {"LoadForceTime",                                 &LoadForceTime},
            {"LoadStrength",                                  &LoadStrength},
            {"PoissonRatio",                                  &PoissonRatio},
            {"SampleCharLength",                              &SampleCharLength},
            {"charStretchEnergyDensityScale",                 &charStretchEnergyDensityScale},
            {"charStretchEnergyScale",                        &charStretchEnergyScale},
            {"SheetThickness",                                &SheetThickness},
            {"ShearModulus",                                  &ShearModulus},
            {"YoungsModulus",                                 &YoungsModulus},
            {"InitDensity",                                   &InitDensity},
            {"ApproxMinInitElemSize",                         &ApproxMinInitElemSize},
            {"charForceScale",                                &charForceScale},
            {"PatchMatrixDimensionlessConditioningThreshold", &PatchMatrixDimensionlessConditioningThreshold}
    };
    double *parameterAddress = doubleMap[parameterName];
    return parameterAddress;
}


int *Settings::getParameterAddressInt(const std::string &parameterName) {
    std::map<std::string, int *> intMap{
            {"testTriangle",          &testTriangle},
            {"slideJustReachedEquil", &slideJustReachedEquil},
            {"NumNodes",              &NumNodes},
            {"NumTriangles",          &NumTriangles},
            {"NumEdges",              &NumEdges},
            {"numBoundaryNodes",      &numBoundaryNodes},
            {"numberOfCores",         &numberOfCores}
    };
    int *parameterAddress = intMap[parameterName];
    return parameterAddress;
}


bool *Settings::getParameterAddressBool(const std::string &parameterName) {
    std::map<std::string, bool *> boolMap{
            {"isDialingDisabled",                                                               &isDialingDisabled},
            {"GlassCones",                                                                      &GlassCones},
            {"isControlledForceEnabled",                                                        &isControlledForceEnabled},
            {"isSeideDeformationsEnabled",                                                      &isSeideDeformationsEnabled},
            {"isLCEModeEnabled",                                                                &isLCEModeEnabled},
            {"isEnergyDensitiesPrinted",                                                        &isEnergyDensitiesPrinted},
            {"isTriangleAreasPrinted",                                                          &isTriangleAreasPrinted},
            {"isAngleDeficitsPrinted",                                                          &isAngleDeficitsPrinted},
            {"isGradientDescentDynamicsEnabled",                                                &isGradientDescentDynamicsEnabled},
            {"ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF", &ForInitialPortionOfProgTensorsSequence_DialProgTauButJumpProgMetricAndProgSecFF},
            {"isDialingFromAnsatzEnabled",                                                      &isDialingFromAnsatzEnabled},
            {"isPerturbationOfInitialPositionsEnabled",                                         &isPerturbationOfInitialPositionsEnabled},
            {"isBoundaryClamped",                                                               &isBoundaryClamped}
    };
    bool *parameterAddress = boolMap[parameterName];
    return parameterAddress;
}


//long int &Settings::getParameterAddressLongInt(const std::string &parameterName) {
//    std::map<std::string, long int &> longIntMap{
//            {"InversePrintRate", InversePrintRate}
//    };
//    long int &parameterAddress = longIntMap[parameterName];
//    return parameterAddress;
//}