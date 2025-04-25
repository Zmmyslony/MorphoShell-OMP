//
// Created by Michał Zmyślony on 12/02/2024.
//

#ifndef SHELLOMORPH_SIMULATION_H
#define SHELLOMORPH_SIMULATION_H

#include <iostream>
#include <fstream>

#include "Node.hpp"
#include "Triangle.hpp"
#include "Edge.hpp"
#include "SimulationStatus.hpp"
#include "settings_new.h"

template<class CharT, class Traits = std::char_traits<CharT> >

struct teestream : std::basic_streambuf<CharT, Traits> {

private:
    std::basic_streambuf<CharT, Traits> *m_rdbuf1;
    std::basic_streambuf<CharT, Traits> *m_rdbuf2;

public:
    teestream(std::basic_streambuf<CharT, Traits> *rdbuf1, std::basic_streambuf<CharT, Traits> *rdbuf2)
            : m_rdbuf1(rdbuf1), m_rdbuf2(rdbuf2) {}

    ~teestream() {
        m_rdbuf1->pubsync();
        m_rdbuf2->pubsync();
    }

protected:
    typename std::basic_streambuf<CharT, Traits>::int_type overflow(typename std::basic_streambuf<CharT, Traits>::int_type ch = Traits::eof()) override {
        typename std::basic_streambuf<CharT, Traits>::int_type result = m_rdbuf1->sputc(ch);
        if (result != Traits::eof()) {
            result = m_rdbuf2->sputc(ch);
        }
        return result;
    }

    virtual int sync() override {
        int result = m_rdbuf1->pubsync();
        if (result == 0) {
            result = m_rdbuf2->pubsync();
        }
        return result;
    }

};

class Simulation {
    std::ofstream out;

    std::string init_string;
    std::string output_dir_name;
//    Settings settings;
    SettingsNew settings_new;
    std::vector<Slide> slides;
    std::vector<Cone> cones;
    std::streambuf *cout_buf = std::cout.rdbuf();

    int num_nodes = 0;
    int num_triangles = 0;
    int num_edges = 0;
    int step_count = 0;
    double damping_power_loss = 0;
    std::vector<Eigen::Vector3d> forcesForEachTriangle;

    double time_global = 0;
    // Reset to zero every time equilibrium is checked, or every time a new DialInFactor value is reached.
    double time_equilibriation = 0;
    // Reset to zero every time equilibrium is reached and a new "dialling in" phase starts.
    double time_phase = 0;
    // Keeps track of which DialInFactor value in dial_in_phases was last dialled *from* (not held at, which is the next value along in the list)
    std::size_t phase_counter = 0;
    // Keeps track of current value
    double dial_in_factor = 0;
    // Smallest altitude of all triangles.
    double characteristic_short_length = 0;
    // Perimeter of the object
    double characteristic_long_length = 0;
    // Smallest altitude of all triangles divided over sqrt(tau)
    double characteristic_length_over_tau = 0;

    // For LCE mode.
    double lambda = 1;

    std::string settings_filename;
    std::string initialisation_filename;
    std::string ansatz_filename = "no_ansatz_file";
    std::string log_filename;

    std::vector<Node> nodes;
    std::vector<Triangle> triangles;
    std::vector<Edge> edges;
    int stage_count = 0;

    SimulationStatus simulation_status = WaitingForEquilibrium;

    std::size_t initial_stage = 0;
    // This is the equivalent variable for dial_in_factor.
    double dialInFactorToStartFrom = 0.0;
    // Create container to store node ansatz positions, if given.
    std::vector<Eigen::Vector3d> nodeAnsatzPositions;
    bool is_equilibrium_seeked = false;

    // Vector holding which triangles correspond with which nodes.
    std::vector<std::vector<std::pair<int, int>>> correspondingTrianglesForNodes;

    // Dial-in factors which divide the simulation into phases between which the equilibrium is reached.
    std::vector<double> dial_in_phases;

    std::vector<double> gaussCurvatures;
    std::vector<double> meanCurvatures;
    std::vector<double> stretchEnergies;
    std::vector<double> bendEnergies;
    std::vector<double> stretchEnergyDensities;
    std::vector<double> bendEnergyDensities;
    std::vector<double> kineticEnergies;
    std::vector<double> strainMeasures;
    std::vector<Eigen::Vector2d> cauchyStressEigenvals;
    std::vector<Eigen::Matrix<double, 3, 2>> cauchyStressEigenvecs;

    std::vector<double> angleDeficits;
    std::vector<double> interiorNodeAngleDeficits;
    std::vector<double> boundaryNodeAngleDeficits;


    void setupLogging();

    void setup_filenames(int argc, char *argv[]);

//    void read_settings();

    void read_vtk_data(const CoreConfig &config);

    void configure_nodes() const;

    /** Remaining geometry/topology/mesh-related things. */
    void configure_topological_properties();

//    void count_boundary_nodes();

    void configure_triangles();

    /* Now calculate a 3x6 matrix for each triangle that can be stored and used
    repeatedly to obtain the elements of the Second Fundamental Form for each
    triangle during the simulation. The 2nd F.F. is an important extrinsic geometric
    object, which completely specifies a surface (up to rigid motions) when combined
    with a metric tensor (aka 1st F.F.). This is called the Bonnet Theorem. An
    'intended' 2nd F.F. can be imposed with a nematic director that 'twists' through
    the thickness of an LCE sheet, and thus the intended and actual 2nd F.F.s enter
    into the energy that is being minimised in this simulation. Intuitively, the 2nd
    F.F. describes the way the deformed surface is curved in 3D space in terms of
    the directional derivative of the surface normal with respect to the 2D
    parametrisation coordinates on the surface. Its matrix representation is
    symmetric.
    A related object is found by 'raising an index' by left-multiplying by the
    inverse metric, which gives the 'Shape Operator'. This has determinant and trace
    equal to the Gaussian curvature and twice the mean curvature respectively.
    The approach taken to approximate the 2nd F.F. for a triangulated mesh in this
    code is based on 'patch fitting', which is simple, and has been shown to
    converge. Consider a triangle t. It has three vertices, and three other nearby
    vertices are selected, subject to the following quadratic fitting being
    well-conditioned.
    The unique quadratic surface (a function of the 2D parametrisation coords) that
    goes through all 6 of these points in the current state is then found and
    assigned to triangle t. The 2nd F.F. for this 'patch' could be calculated in
    full and averaged over the patch, although this complicates matters. Instead, we
    just take the normal to the triangle face and combine that with the 2nd
    derivatives of the patch surface following a particular expression for the secFF.
    The matrix pre-calculated here is later multiplied onto the current node
    positions of the 6 patch nodes to give patch surface derivatives.
    */
    void set_node_patches();

    /*For the current and initial in-plane coordinate bases to match up, the node
    labels for each triangle must be ordered such that the resulting calculated
    normals for the initial flat state point along the +z axis, not -z. To do this
    I calculate the normals for an arbitrary label order, then loop over triangles,
    shuffling the labels for those with negative-z-component normals, before
    recalculating the initial normals, areas etc. Note, a more complex procedure
    will be required for non-planar initial geometries. The magic numbers in the
    argument list are just debugging values that should not matter here.*/
    void orient_node_labels();

    void set_initial_conditions();

    /* First find the approximate smallest linear size of mesh element, based on
    smallest altitude of each triangle. Find also the smallest value of
    this linear size / sqrt(progTau), which is what really matters for the stretching
    time step. If progTau varied a lot over the programmed_taus, you might
    want to be less wasteful and recalculate the time step each time you move on to
    the next set of progTaus in the sequence, because the time step then might be
    much bigger for other sets in the sequence. If you wanted to get even fancier,
    you could change the time step as prog_Tau is dialled in between two sets in the
    sequence.*/
    void find_smallest_element();

//    void print_total_load_force();

    void setup_characteristic_scales();
    void setup_equilibrium_dial_in_factors();
    void setInitialTriangleElongations();

    /* Create std::vectors (with one element for each triangle) that will be passed
    by reference to functions calculating the Gauss and mean curvatures, and
    stretching and bending energies and energy densities. */
    void initialise_simulation_vectors();

    void run_ansatz(int counter);

    void run_tensor_increment(int stage_counter);

//    void impose_seide_deformation(double s, const std::vector<Eigen::Vector3d> &nodeUnstressedConePosits);
//
//    void update_slide_properties();
//
//    void setup_glass_cones(int highest_node, int lowest_node);

    void setup_imposed_seide_deformations(double &s1, int highest_node, int lowest_node,
                                          std::vector<Eigen::Vector3d> &nodeUnstressedConePosits);

    void first_step_configuration();

    long long int progress_single_step(int counter);

    void update_dial_in_factor();

    void save_and_print_details(int counter, long long int duration_us);

    void begin_equilibrium_search(int counter);

    void check_for_equilibrium();

    void error_large_force(int counter);

    void setup_reached_equilibrium();

    std::string log_prefix() const;

    void init(int argc, char *argv[]);

    void read_settings_new(int argc, char **argv);

    void updateTriangleProperties(int counter);

    void advance_time();

    // Adds damping, additional physics e.g. slides, cones, and applies boundary condition.
    void add_non_elastic_forces();

    // Calculates and adds elastic forces to each
    void add_elastic_forces();

    void advance_physics();

    void setup_tensor_increment(int stage_counter);

    void check_if_equilibrium_search_begun(int stage_counter);

    bool isDataPrinted();

    long long int export_vtk(int counter);

public :

    Simulation(int argc, char *argv[]);

    int run_simulation();

    void equilibriumTest(int stage_counter, long long int duration_us);
};


#endif //SHELLOMORPH_SIMULATION_H
