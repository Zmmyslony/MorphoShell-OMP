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

This is the header file for the struct that that will be hold settings
and user-set parameters, and passed around between functions etc. It also holds
some values which are not settings as such, like the number of nodes (determined
by initial data). It should always be passed by reference, and almost always
const (though occasionally one could want to change 'settings' within a function
- e.g. time step may change.*/

#ifndef _SETTINGS_TAG_
#define _SETTINGS_TAG_

#include <string>
#include <vector>

#include "CustomOutStreamClass.hpp"
//#include "Triangle.hpp"

struct Settings {

    // FOR CONE SQUASHING/BUCKLING BETWEEN TWO SLIDES.
    // Z-coordinates of bottom and top slides.
    int number_of_cores = 1;
    double gravity_sign;
    double init_slide_z_coord_lower;
    double init_slide_z_coord_upper;
    double curr_slide_z_coord_upper;
    // Spring stiffness of 'glass' slides (modelled as linear springs) is this
    // prefactor * mu * t
    double slide_stiffness_prefactor;
    // Slide speed is this prefactor * char longest-wavelength mode speed.
    double slide_speed_prefactor;
    // Upper slide's downwards displacement from its initial position.
    double upper_slide_displacement;
    // Coefficient of slide friction
    double slide_friction_coefficient;
    // If using an ansatz, all nodes less than this many thicknesses vertically
    // above the node with the lowest initial (ansatz) z coord will be clamped.
    // Set to negative value to turn off.
    double thicknesses_above_lowest_node_to_clamp_up_to;
    double bending_long_time;
    double damping_scale;
    double smallest_size_over_root_tau;
    bool is_dialing_disabled;
    bool glass_cones;
    double cone_angle;
    double lambda; //LCE lambda.
    int test_triangle; // which triangle to print out stress for.

    bool is_controlled_force_enabled; // If true do controlled-force experiment rather than controlled-displacement.
    double slide_weight_dial_speed_fac; // Control how big the jumps in weight are every time slide equil is reached.
    double total_slide_force_to_mu_t_sq_ratio_equil_threshold; // Control how fussy the slide is about being in equilibrium (material equilibrium checked for separately).
    double initial_slide_weight_for_ctrld_force_in_units_of_mu_tsq; // Initial weight of slide in units of mu t^2, when doing a controlled force experiment.
    int is_slide_just_equilibrated; // 0 unless equil just reached, then set to 1 to trigger printout, then immediately returned to 0.
    double slide_damping_param; // Ratio between total force on upper slide, and it's velocity (we do viscous dynamics for the slide).
    double upper_slide_weight; // Slide weight as multiple of charForceScale. If doing controlled-force experiment, slideSpeedPrefactor set's dial speed of this weight rather than slide displacement.
    // Upper slide's downwards velocity, in a controlled-force experiment.
    double upper_slide_vel;
    // Total upwards force on upper slide, only needs to be stored here for controlled-force experiment.
    double upper_tot_slide_force;
    // To instead hold slide weight constant while you change lambda, set this to be
    // positive. The weight will be this times mu t^2. If this is negative instead slide
    // weight is varied according to the other settings.
    double const_slide_weight_fac;
    double spacer_height; // Height of spacer in units of thickness.

    bool is_seide_deformations_enabled; // In regions near the top and bottom rim, hopefully better at getting Seide buckling result than glass squashing.
    double p; // Seide's load force.
    double p_speed_prefactor; // Control how quickly p increases from 0.


    // To allow specifying initial slide positions directly in settings file.
    double specify_init_slide_z_coord_upper;
    double specify_init_slide_z_coord_lower;

    /* Inverse print rate - results data will be saved in .vtk file every
    InversePrintRate steps. This is based on the dial-in time and the time step,
    and should be tuned via the PrintFrequency setting. */
    long int inverse_print_rate;

    /* Dimensionless tunable setting, that sets roughly how many output files
    will be generated over the course of one dial-in step. Output will be printed
    every DialInStepTime / Timestep, rounded to the nearest integer > 0. Must be
    specified with a decimal point as it is a double. Switch off these prints
    entirely by setting PrintFrequency to a negative number. Currently,
    PrintFrequency values in [0.0, 1.0) are not permitted. */
    double print_frequency;

    /* If true, the LCE parameter is dialled directly, rather than the components
    of the programmed metric (which are non-linear functions of lambda). This is
    useful for making physically correct evolution videos, but less so for efficiently
    reaching final states because it does not play well with isDialingFromAnsatzEnabled.*/
    bool is_lce_mode_enabled;

    /* Bool to determine whether the output files include stretching and
    bending energy densities for each triangle in addition to curvatures.*/
    bool is_energy_densities_printed;

    /* Bool to determine whether the output files include triangle areas.
    */
    bool is_triangle_areas_printed;

    /* Bool to determine whether the output files include the 'angular
    deficit' for each node: 2*Pi - Sum(incident triangle angles).*/
    bool is_angle_deficits_printed;

    // Total numbers of nodes, triangles, and edges.
    int num_nodes, num_triangles, num_edges;

    /* Numerical damping factor, such that damping force on a node is equal to
     -mass*velocity*NumDampFactor/InitDensity
    The value is chosen based on a toy stretching + bending model, to
    approximately critically damp the longest wavelength bending mode. The
    value can be tuned from the settings file using two dimensionless prefactors.
    The first is used during dialling in, the second is used instead when waiting
    for equilibrium.
    */
    double num_damp_factor;
    double dial_in_damping;
    double equilibriation_damping;

    /* Dimensionless number which modifies the bending energy density in the
    *spirit* of the Gent 3D elastic energy. This parameter is *analogous* to the
    original Gent parameter - not the same! Let the usual bend energy density
    that is ~ curvature^2 be called Wb. The energy density modification is
    such that for Wb below some scale (i.e. for suitably low curvatures), the
    bend energy density used is Wb to a very good approx. If Wb is bigger than
    some scale, the bend energy density is taken to be a function of Wb, such
    that it grows much more sharply with curvature. This is intended to reflect
    the physical fact that curvatures on the scale of the thickness should be
    much more strongly penalised than they are by Wb. Thus, the energy
    density scale where the behaviour transitions should be something like Wb
    evaluated for a curvature of 1/thickness. We take it to be that scale,
    multiplied by the tunable GentFactor below, which should therefore be rougly
    O(1). We generally use values between 1 and 10.*/
    double gent_factor;

    /* Prefactor in time step calculation, chosen
empirically to avoid `blowing up'. */
    double time_step_prefactor;

    // Time step
    double time_step;

    /* Force (largest) to characteristic force and speed (largest) to sound
    speed ratios below which system is considered to be in equilbirum. */
    double char_force, char_speed;

    /* Bool determining whether Gradient Descent dynamics are used as opposed to
    the default Newtonian Dynamics. Gradient Descent is equivalent to fully
    overdamped Newtonian Dynamics, i.e. with no inertial term in the equation
    of motion, giving (for linear damping):
    velocity = non-damping force / damping coefficient,
    where damping coefficient = mass*NumDampFactor/InitDensity in this code.*/
    bool is_gradient_descent_dynamics_enabled;

    /* When this bool is set to true, the progMetric and
    progSecFF that we would normally be dialling FROM at the start of the run,
    (the flat sheet trivial programmed tensors, or *possibly* those from a later
    point in the progTensors sequence if an ansatz was used) are
    just set to the values they would be dialling TO, so that no dialling
    of these quantities happens between the initial (starting) entry in
    the progTensor sequence and the next. progTau however IS dialled in as
    normal. Upon reaching the next entry in the sequence, normal behaviour
    resumes so there will be a phase of waiting for equilibrium,
    followed potentially by normal dialling in of further progTensors in
    the sequence, if there are anymore in the sequence.
    To use it, you probably want to supply an ansatz file with its
    currDialInFactor CHANGED TO = 0. See main() for more details, including the
    intended use case. */
    bool for_initial_portion_of_prog_tensors_sequence_dial_prog_tau_but_jump_prog_metric_and_prog_sec_ff;

    /* If settings.isDialingFromAnsatzEnabled == true AND an ansatz data file
    is being used, the state of dial-in given by the header line in the ansatz
    file is ignored, and instead a new dialling in phase is started from
    whatever the *calculated* metric and secFF are in the ansatz state.
    Tau is not a calculable geometric quantity so that we still dial in from
    1.0 as usual. I can see an argument for just setting it to the programmed
    value immediately too, but it should make almost no difference anyway.
    This setting is very handy if you can only make an
    approximate guess at a solution, and want to start from there and see
    where things go without explosive dynamics.*/
    bool is_dialing_from_ansatz_enabled;


    /* Variables used to control the programmed properties' time variation.
    This is not necessarily directly physical, especially for a fast optical
    activation instead of slow heating; the idea here is to gradually `dial in'
    the programmed (energetically favoured) metric and second fundamental form.
    This prevents anything too explosive happening the code, and can be used to
    make the simulation quasistatic (which can help avoid undesired isometries
    for example), though it may not be necessary.
    The range between 0, and 1 is split into a sequence of discrete values based
    on the chosen settings.DialInResolution. The dial-in factor will then
    evolve linearly from one such value to the next, over
    settings.DialInStepTime. After each such step, the value is held constant
    until equilibrium is reached (to within the chosen tolerance). In the
    'waiting' stage, the check for equilibrium occurs at a rate determined by
    the tunable (and dimensionless)
    settings.TimeBetweenEquilChecksPrefactor in the settings file.
    TimeBetweenEquilChecks will be set to this value * settings.DialInStepTime.
    settings. DialInStepTime is determined by the code to approximately
    equal the characteristic decay time for the longest wavelength bending mode
    when it is critically damped. The value can be tuned in the settings file
    with a dimensionless prefactor. */
    double dial_in_resolution;
    double dial_in_step_time;
    double time_between_equil_checks_prefactor;
    double time_between_equil_checks;
    double dial_in_step_time_prefactor;

    // Amount of time (from start) that 'prod' force should be applied for
    double prod_force_time;

    /* Factor determining strength of prod force. Values less than 0.01 are
    probably recommended. */
    double prod_strength;

    /* Amount of time (from start) that simple downward load force should be
    applied for (to nodes with isLoadForceEnabled = true, set via input data file).*/
    double load_force_time;

    /* Factor determining strength of load force, which has magnitude
    LoadStrength*ShearModulus*ApproxMinInitElemSize*SheetThickness. */
    double load_strength;

    /* Bool determining whether to perturb node positions with a small amount of
    random noise just before dynamical evolution begins, to 'break away' from
    an initial planar geometry, for example. */
    bool is_perturbation_of_initial_positions_enabled;

    /* Bool determining whether to clamp whole boundary in addition to the clamp
    indicators given in the input data file. */
    bool is_boundary_clamped;

    // Mechanical Poisson ratio of sheet material.
    double poisson_ratio;

    /* Approximate length scale of expected minimum wavevector in the system,
    e.g. the longest length scale in the sample. */
    double sample_char_length;

    /* Approximate scales of stretching energy and energy density.*/
    double char_stretch_energy_density_scale;
    // double charBendEnergyDensityScale;
    double char_stretch_energy_scale;
    // double charBendEnergyScale;

    // Initial sheet thickness, which determines the balance between the
    //stretching and bending energies.
    double sheet_thickness;

    // Shear modulus.
    double shear_modulus;

    // Young's Modulus.
    double youngs_modulus;

    // Initial( reference) state density.
    double init_density;

    /* Square root of the smallest initial element area, used in calculating time step
    and characteristic force scale for testing equilibrium. */
    double approx_min_init_elem_size;

    // Number of boundary nodes.
    int num_boundary_nodes;

    /* Characteristic force scale in the simulation, determined from other
    parameters (in main.cpp). */
    double char_force_scale;

    /* Factor to determine the threshold for the allowed condition number of the
    matrices which are inverted in finding the second fundamental form
    coefficients. Any tempMatForPatchSecDerivs with a condition number greater
    than PatchMatrixDimensionlessConditioningThreshold divided by the patch's
    approximate linear size^2 will not be used, and an alternative will be found
    instead. Worse meshes may require higher values to accept any of the
    possible patches. In general, choose as low a value as possible that
    allows all triangles to find a patch successfully.
    */
    double patch_matrix_dimensionless_conditioning_threshold;

    double * getParameterAddressDouble(const std::string &parameterName);

    int * getParameterAddressInt(const std::string &parameterName);

    bool * getParameterAddressBool(const std::string &parameterName);

    void SetupDialIn(CustomOutStreamClass &logStream);

//    long int &getParameterAddressLongInt(const std::string &parameterName);

    void SetupStepTime(CustomOutStreamClass &logStream);

    void SetupPrintFrequency(CustomOutStreamClass &logStream);

//    void SetupSmallestElements(CustomOutStreamClass &logStream, std::vector<Triangle> &triangles,
//                               std::vector<std::vector<double>> &programmed_taus);
};

#endif
