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
#include <filesystem> // Note this requires g++ version >=8
#include "initialDirAndFileHandling.hpp"

std::string initialDirAndFileHandling(
        const std::string &settings_file_name_str,
        const std::string &settings_file_name_str_final_piece,
        const std::string &initial_data_file_name_str,
        const std::string &initial_data_file_name_str_final_piece,
        const std::string &ansatz_data_file_name_str,
        const std::string &ansatz_data_file_name_str_final_piece,
        const int &argc,
        std::string &initWriteToLogStr
) {

    // Get full path of current working directory.
    std::filesystem::path cwd = std::filesystem::current_path();
    initWriteToLogStr += "Current working directory: " + cwd.string() + "\n";

    /* Create strings containing names for the output directory, and a potential
    backup for the previous one. */
    std::string outputDirName("outputfiles_using_" + settings_file_name_str_final_piece);
    std::string prevOutputDirName("prev_outputfiles_using_" + settings_file_name_str_final_piece);

    /* If previous version of outputfiles_using... directory exists for the
    current settings file name, save (rename) it as prev_outputfiles...
    (**OVERWRITTEN EVERY TIME THIS IS DONE!**).*/
    std::filesystem::path toRename(cwd / outputDirName);

    if (std::filesystem::exists(toRename)) {

        /* Delete any existing prev_outputfiles_using_... directory corresponding to
        the current settings file (i.e. if the code saved you from accidentally
        overwriting outputfiles_using_... the first time, but you forgot to rename
        the prev_outputfiles_using_... directory that saved you: unlucky!).*/
        std::filesystem::path toDelete(cwd / prevOutputDirName); // '/' means append.
        if (std::filesystem::exists(toDelete)) {
            std::uintmax_t numObjDeleted = std::filesystem::remove_all(toDelete);
            if (numObjDeleted) {}; // Stops compiler complaining about unused variable.
            initWriteToLogStr +=
                    "Deleted existing output directory from run-before-last: " + prevOutputDirName + "\n";
        }

        // Do renaming of outputfiles_using...
        std::filesystem::rename(toRename, cwd / prevOutputDirName);
        initWriteToLogStr += "Previous outputfiles directory still present.\n";
        initWriteToLogStr +=
                "The previous version was saved this time, but you may want to move it "
                "out of the working directory, because it will be overwritten next time "
                "the code is run with this settings file.\n";
    } else {
        initWriteToLogStr += "No previous outputfiles directory was found that needed saving.\n";
    }


    /* Create new directory to hold output .vtk files, and copy settings and data
    files into a subdirectory of that directory.*/
    std::filesystem::path outputDir(cwd / outputDirName);
    std::filesystem::path inputFilesDir(outputDir / "input_files_used");
    std::filesystem::path settingsFile(cwd / settings_file_name_str);
    std::filesystem::path initDataFile(cwd / initial_data_file_name_str);
    std::filesystem::path ansatzDataFile(cwd / ansatz_data_file_name_str);


    bool isOutputDirMade = std::filesystem::create_directory(outputDir);
    bool isInputDirMade = std::filesystem::create_directory(inputFilesDir);

    bool isSettingsCopied = std::filesystem::copy_file(settingsFile,
                                                       inputFilesDir / settings_file_name_str_final_piece);
    bool isInitCopied = std::filesystem::copy_file(initDataFile,
                                                   inputFilesDir / initial_data_file_name_str_final_piece);
    bool isAnsatzCopied = true; // Default true so no error thrown if no ansatz given.
    if (argc == 4) {
        isAnsatzCopied = std::filesystem::copy_file(ansatzDataFile,
                                                    inputFilesDir / ansatz_data_file_name_str_final_piece);
    }

    if (!(isOutputDirMade && isInputDirMade && isSettingsCopied && isInitCopied && isAnsatzCopied)) {
        // This exception on top of the std::filesystem ones may be overkill.
        throw std::runtime_error("Error: setting up output files directory failed");
    }

    return outputDirName;
}
