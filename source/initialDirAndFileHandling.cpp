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

Function to set up data directories and copy input files into them to retain
a record of which settings and data were used for the run. If a directory
already exists with the same name as the output directory being set up, the old
one is renamed to avoid accidental deletion. Any previous existing version of
this 'previous output' directory is deleted. */

#include <iostream>

#include <string>
#include <boost/filesystem/operations.hpp>
#include "initialDirAndFileHandling.hpp"

fs::path
directory_setup(std::string &initWriteToLogStr, const std::vector<fs::path> &config_paths,
                const std::vector<fs::path> &vtk_paths) {
    
    fs::path cwd = fs::current_path();
    initWriteToLogStr += "Current working directory: " + cwd.string() + "\n";
    
    fs::path simulation_vtk = vtk_paths.front();
    fs::path output_directory = cwd / (simulation_vtk.stem().string() + "_output");
    
    if (fs::exists(output_directory)) {
        fs::path prev_output_directory = cwd / (simulation_vtk.stem().string() + "_previous_output");
        fs::remove_all(prev_output_directory);
        fs::rename(output_directory, prev_output_directory);
    }

    fs::path input_files_dir(output_directory / "input_files");

    if (!fs::create_directory(output_directory)) {throw std::runtime_error("Error creating the output directory.");}
    if (!fs::create_directory(input_files_dir)) {throw std::runtime_error("Error creating the input files directory.");}

    for (auto config : config_paths) {
        if (!fs::exists(config)) {throw (std::runtime_error(config.string() + " does not exist."));}
        if (!fs::copy_file(config, input_files_dir / config.filename())) {
            throw std::runtime_error("Error while copying " + config.string());
        }
    }

    for (auto vtk : vtk_paths) {
        if (!fs::exists(vtk)) {throw (std::runtime_error(vtk.string() + " does not exist."));}
        if (!fs::copy_file(vtk, input_files_dir / vtk.filename())) {
            throw std::runtime_error("Error while copying " + vtk.string());
        }
    }
    return output_directory;
}
