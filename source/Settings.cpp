//
// Created by Michał Zmyślony on 22/03/2022.
//
#define _USE_MATH_DEFINES

#include "Settings.hpp"

#include <map>
#include <iostream>
#include <cmath>
#include <math.h>

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
            {"BendingLongTime",                               &bending_long_time},
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
            {"PatchMatrixDimensionlessConditioningThreshold", &PatchMatrixDimensionlessConditioningThreshold},
            {"GravitySign",                                   &gravity_sign}
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


void Settings::SetupDialIn(CustomOutStreamClass &logStream) {
    /* Calculate 'dialling in' time and damping coefficient based on toy model
stretching and bending analyses. The 'Long times' are approximate characteristic
times for the longest-wavelength modes in the system for stretching and bending.
The damping coefficient is chosen to approximately critically damp the longest-
wavelength bending mode.
*/
    const double PI = M_PI;

    logStream.open();
    double minWavevector = 2 * PI / SampleCharLength;
    dampingScale =
            2 * sqrt(InitDensity * ShearModulus / (6 * (1 - PoissonRatio))) * SheetThickness * minWavevector *
            minWavevector;
    NumDampFactor = DampingPrefactor1 * dampingScale;
    logStream << "Numerical Damping Coefficient = " << NumDampFactor << std::endl;

    double stretchingLongTime = sqrt(InitDensity / ShearModulus) / minWavevector;
    bending_long_time = sqrt(InitDensity * 6 * (1 - PoissonRatio) / ShearModulus) /
                        (SheetThickness * minWavevector * minWavevector);

    if (stretchingLongTime > bending_long_time) {
        DialInStepTime = DialInStepTimePrefactor * stretchingLongTime;
        logStream << "Dial-in step time = " << DialInStepTime
                  << ", set based on stretching rather than bending." << std::endl;
    } else {
        DialInStepTime = DialInStepTimePrefactor * bending_long_time;
        logStream << "Dial-in step time = " << DialInStepTime
                  << ", set based on bending rather than stretching." << std::endl;
    }
    logStream.close();
}


void Settings::SetupStepTime(CustomOutStreamClass &logStream) {
    /* Calculate time step based on toy model stretching and bending analyses (take
 whichever gives shortest characteristic time), and print. */
    double stretchingTimeStep;
    double bendingTimeStep;

    if (!isGradientDescentDynamicsEnabled) {
        stretchingTimeStep = TimeStepPrefactor * sqrt(InitDensity / ShearModulus) *
                             smallestSizeOverRootTau;
        bendingTimeStep = TimeStepPrefactor * sqrt(InitDensity * 6 * (1 - PoissonRatio) / ShearModulus) *
                          (ApproxMinInitElemSize / SheetThickness) * ApproxMinInitElemSize / (2 * M_PI);
    } else {
        /* Use second damping factor, since we should only be using gradient descent
        for already dialled-in states. Note the damping factor cancels out in the
        dynamics and has no effect - but only if you are consistent with which one
        you use!*/
        double numDampFac2 = DampingPrefactor2 * dampingScale;
        stretchingTimeStep =
                TimeStepPrefactor * (numDampFac2 / ShearModulus) * pow(smallestSizeOverRootTau / (2 * M_PI), 2);
        bendingTimeStep =
                TimeStepPrefactor * (6 * (1 - PoissonRatio) * numDampFac2 / (pow(SheetThickness, 2) * ShearModulus))
                * pow(smallestSizeOverRootTau / (2 * M_PI), 4);

        logStream.open();
        logStream
                << "\nUSING GRADIENT DESCENT DYNAMICS.\nBE WARNED - this feature is only intended for use on nearly converged states supplied as fully dialled-in ansatzes."
                << "\nThe number of timesteps required for a full simulation with gradient descent is expected to be huge (~10^8).\nAs such, the gradient descent long bending timescale "
                << "is not even used for DialInStepTime,\nas you really shouldn't be doing any dialling in!\nDialInStepTime is instead set such that usual non-gradient-descent settings"
                << " should still give reasonable printout and and equilibrium check frequencies." << std::endl;
        logStream.close();
    }

    logStream.open();
    logStream << "Short stretching and bending timescales: " << stretchingTimeStep << ", " << bendingTimeStep
              << std::endl;

    if (stretchingTimeStep < bendingTimeStep) {
        TimeStep = stretchingTimeStep;
        logStream << "Time step = " << TimeStep << ", set based on stretching rather than bending."
                  << std::endl;
    } else {
        TimeStep = bendingTimeStep;
        logStream << "Time step = " << TimeStep << ", set based on bending rather than stretching."
                  << std::endl;
    }
    logStream.close();

    /* Fix up DialInStepTime if gradient descent is being used, just to give
reasonable printout and equilibrium check frequencies without extreme settings.*/
    if (isGradientDescentDynamicsEnabled) {
        DialInStepTime = 1000 * TimeStep;
        TimeBetweenEquilChecks = TimeBetweenEquilChecksPrefactor * DialInStepTime;
    }
}


void Settings::SetupPrintFrequency(CustomOutStreamClass &logStream) {
    /* Calculate print frequency based on DialInStepTime / TimeStep, tuned by a
dimensionless parameter in the setting file and rounded to the nearest integer.
This rounding should work fine as long as the PrintFrequency isn't ridiculously
huge. To avoid this we require that PrintFrequency is not in [0,1). If
PrintFrequency is very large (more likely), so that InversePrintRate rounds to 0,
we set InversePrintRate to 1 to get a printout after every time step. This can
be a useful thing to do sometimes. If PrintFrequency is set to a negative value
we instead set InversePrintRate to -1, which means no output files are regularly
written at all.
*/
    if (PrintFrequency < 0) {
        InversePrintRate = -1;
    } else {
        if (PrintFrequency < 1) {
            logStream
                    << "Error: To avoid potential divide-by-zero problems we do NOT allow settings.PrintFrequency to lie in [0, 1). This is overkill somewhat, and could be relaxed if need be."
                    << std::endl;
            return;
        } else {
            InversePrintRate = lround(DialInStepTime / (TimeStep * PrintFrequency));
            if (InversePrintRate == 0) {
                InversePrintRate = 1;
            }
        }
    }
}

//void Settings::SetupSmallestElements(CustomOutStreamClass &logStream, std::vector<Triangle> &triangles,
//                                     std::vector<std::vector<double>> &programmed_taus) {
//
//}