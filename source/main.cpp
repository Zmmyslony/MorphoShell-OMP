/* Author: Daniel Duffy, University of Cambridge, dld34@cam.ac.uk

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
along with Shellmorph.  If not, see <https://www.gnu.org/licenses/>.`
/////////////////////////////////////////////////////

Libconfig++ is distributed under the Lesser GPL license (2.1 or later), and copyright is held by Mark A Lindner (at least in large part).
Version 1.7.2 was used for this code, which can be accessed via: http://hyperrealm.github.io/libconfig/

Eigen is distributed under the Mozilla Public License v. 2.0 , and copyright is held by Gael Guennebaud and Benoit Jacob (at least in large part).
Version 3.3.7 was used for this code, which can be accessed via: http://eigen.tuxfamily.org/

//////////////////////////////////////////////////////////////////////////

Main file for Shellmorph, the main code developed for my PhD.

This code simulates a thin shape-morphing 2D sheet with a programmed metric and
second fundamental form.
The main object in the simulations is an unstructured triangulated mesh, with
data stored in three std:vectors: 'nodes', 'triangles' and 'edges', each element
of which holds a class, in turn holding the data for that node/triangle/edge.
So for example nodes[0] is a Node instance holding the data for the first
node, which is labelled zero, matching the indexing of 'nodes'. The initial node
and triangle data is read in from a VTK Legacy file, which is also the output
format. Paraview is recommended for visualising these files.

The code executable takes either two or three command line arguments. The first
is a settings file, and the second is an input data file containing the planar
2D initial (reference) state mesh and the programmed quantities. The optional
third argument is an `ansatz' file (also .vtk), which can be an output file from
a previous simulation. If this is provided, the nodes will be moved to the
positions specified in the ansatz file before dynamics begins from that state.
Thus you can terminate a simulation and restart it from where it left off with
a new settings file (though the velocities start from zero upon restarting).

An elastic energy that penalises stretching and bending deviations from the
programmed metric and secondfundamental form is minimised (using analytical
derivatives) via either linearly-damped Newtonian dynamics or gradient descent.
The former is on the whole recommended, while the latter can be useful for
smoothing out meshes that have significant angles between triangles i.e. small
wavelength oscillations. Equilibrium is deemed achieved when two dimensionless
numbers fall below thresholds specified in the settings file. These numbers are
maxima (over the mesh) of non-dimensionalised node speed and elastic force.
*/

// Turn Eigen bounds checking off for speed (after running with checks naturally).
#ifndef EIGEN_NO_DEBUG
#define EIGEN_NO_DEBUG
#endif

#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#define _USE_MATH_DEFINES

#include <iomanip>

#include "simulation.h"

std::pair<double, double> find_timing_fit(const std::vector<long long>& thread_timings)
{
    /** Fits a function: t_constant + t_parallel / num_threads, to a list of
     * timings.
     **/

    double sum = 0;
    double sum_inverse = 0;
    double nInverse = 0;
    double nInverseSquared = 0;
    auto n = static_cast<double>(thread_timings.size());
    for (int i = 0; i < thread_timings.size() - 1; i++)
    {
        sum += static_cast<double>(thread_timings[i]);
        sum_inverse += static_cast<double>(thread_timings[i] / (i + 1));
        nInverse += 1. / (i + 1);
        nInverseSquared += 1. / pow(i + 1, 2);
    }
    double denominator = n * nInverseSquared - nInverse * nInverse;
    double constant_duration = (nInverseSquared * sum - nInverse * sum_inverse) / denominator;
    double parallel_duration = (-nInverse * sum + n * sum_inverse) / denominator;

    return {constant_duration, parallel_duration};
}

void benchmark_multithreading(int argc, char* argv[], unsigned short max_thread_count)
{
    std::vector<long long> full_simulation_durations;
    std::vector<long long> mechanics_simulation_durations;
    std::vector<long long> exporting_simulation_durations;
    std::vector<long long> unaccounted_simulation_durations;

    for (int i = 1; i <= max_thread_count; i++)
    {
        Simulation simulation = Simulation(argc, argv, i);
        simulation.run_simulation();
        full_simulation_durations.push_back(simulation.benchmarking_full_duration / 1000);
        mechanics_simulation_durations.push_back(simulation.benchmarking_mechanics_duration / 1000);
        exporting_simulation_durations.push_back(simulation.benchmarking_export_duration / 1000);
        unaccounted_simulation_durations.push_back(
            (simulation.benchmarking_full_duration - simulation.benchmarking_mechanics_duration - simulation.
                benchmarking_export_duration) / 1000);
    }


    std::pair<double, double> durations_full = find_timing_fit(full_simulation_durations);
    std::pair<double, double> duration_mechanics = find_timing_fit(mechanics_simulation_durations);
    std::pair<double, double> duration_exporting = find_timing_fit(exporting_simulation_durations);
    std::pair<double, double> duration_unaccounted = find_timing_fit(unaccounted_simulation_durations);


    std::cout << std::endl;
    std::cout << std::defaultfloat << std::setprecision(2);
    std::cout << "========================================" << std::endl;
    std::cout << "======= Benchmark results below ========" << std::endl;
    std::cout << "========================================" << std::endl;

    double mechanics_ratio = (duration_mechanics.first + duration_mechanics.second) / (durations_full.first +
        durations_full.second) * 100;
    double export_ratio = (duration_exporting.first + duration_exporting.second) / (durations_full.first +
        durations_full.second)* 100;
    double unaccounted_ratio = (duration_unaccounted.first + duration_unaccounted.second) / (durations_full.first +
        durations_full.second) * 100;

    std::cout << "Category \tamount [%] \toverhead [ms] \tt_parallel [ms]\tparallelisation [%]" << std::endl;
    std::printf("Total \t\t%.2f \t\t%.0f \t\t\t%.0f \t\t\t%.2f\n", 100., durations_full.first, durations_full.second, (durations_full.second) / (durations_full.first + durations_full.second) * 100);
    std::printf("Mechanics \t%.2f \t\t%.0f \t\t\t%.0f \t\t\t%.2f\n", mechanics_ratio, duration_mechanics.first, duration_mechanics.second, (duration_mechanics.second) / (duration_mechanics.first + duration_mechanics.second) * 100);
    std::printf("Exporting \t%.2f \t\t%.0f \t\t\t%.0f \t\t\t%.2f\n", export_ratio, duration_exporting.first, duration_exporting.second, (duration_exporting.second) / (duration_exporting.first + duration_exporting.second) * 100);
    std::printf("Unaccounted\t%.2f \t\t%.0f \t\t\t%.0f \t\t\t%.2f\n", unaccounted_ratio, duration_unaccounted.first, duration_unaccounted.second, (duration_unaccounted.second) / (duration_unaccounted.first + duration_unaccounted.second) * 100);
}

int main(int argc, char* argv[])
{
    // Simulation simulation(argc, argv, -1);
    // return simulation.run_simulation();
    benchmark_multithreading(argc, argv, 8);
}
