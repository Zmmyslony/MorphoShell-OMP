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


//    std::string outputDirName("outputfiles_using_" + settings_file_name_str_final_piece);
//    std::string prevOutputDirName("prev_outputfiles_using_" + settings_file_name_str_final_piece);
//
//    /* If previous version of outputfiles_using... directory exists for the
//    current settings file name, save (rename) it as prev_outputfiles...
//    (**OVERWRITTEN EVERY TIME THIS IS DONE!**).*/
//    fs::path toRename(cwd / outputDirName);
//
//    if (fs::exists(toRename)) {
//
//        /* Delete any existing prev_outputfiles_using_... directory corresponding to
//        the current settings file (i.e. if the code saved you from accidentally
//        overwriting outputfiles_using_... the first time, but you forgot to rename
//        the prev_outputfiles_using_... directory that saved you: unlucky!).*/
//        fs::path toDelete(cwd / prevOutputDirName); // '/' means append.
//        if (fs::exists(toDelete)) {
//            std::uintmax_t numObjDeleted = fs::remove_all(toDelete);
//            if (numObjDeleted) {}; // Stops compiler complaining about unused variable.
//            initWriteToLogStr +=
//                    "Deleted existing output directory from run-before-last: " + prevOutputDirName + "\n";
//        }
//
//        // Do renaming of outputfiles_using...
//        fs::rename(toRename, cwd / prevOutputDirName);
//        initWriteToLogStr += "Previous outputfiles directory still present.\n";
//        initWriteToLogStr +=
//                "The previous version was saved this time, but you may want to move it "
//                "out of the working directory, because it will be overwritten next time "
//                "the code is run with this settings file.\n";
//    } else {
//        initWriteToLogStr += "No previous outputfiles directory was found that needed saving.\n";
//    }


    /* Create new directory to hold output .vtk files, and copy settings and data
    files into a subdirectory of that directory.*/
//    fs::path outputDir(cwd / outputDirName);
//    fs::path inputFilesDir(outputDir / "input_files_used");
//    fs::path settingsFile(cwd / settings_file_name_str);
//    fs::path initDataFile(cwd / initial_data_file_name_str);
//    fs::path ansatzDataFile(cwd / ansatz_data_file_name_str);
//
//
//    bool isOutputDirMade = fs::create_directory(outputDir);
//    bool isInputDirMade = fs::create_directory(input_files_dir);
//
//    bool isSettingsCopied = fs::copy_file(settingsFile,
//                                          input_files_dir / settings_file_name_str_final_piece);
//    bool isInitCopied = fs::copy_file(initDataFile,
//                                      input_files_dir / initial_data_file_name_str_final_piece);
//    bool isAnsatzCopied = true; // Default true so no error thrown if no ansatz given.
//    if (argc == 4) {
//        isAnsatzCopied = fs::copy_file(ansatzDataFile,
//                                       input_files_dir / ansatz_data_file_name_str_final_piece);
//    }
//
//    if (!(isOutputDirMade && isInputDirMade && isSettingsCopied && isInitCopied && isAnsatzCopied)) {
//        // This exception on top of the fs ones may be overkill.
//        throw std::runtime_error("Error: setting up output files directory failed");
//    }

    return output_directory;
}
