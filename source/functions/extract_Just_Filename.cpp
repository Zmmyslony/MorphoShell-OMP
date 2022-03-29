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

Function to extract a filename from a full path, or just return the file name
if just a file name is given. The input is a std::string, and the function works
by taking just the piece of the string after the last / character, if any / is
present. */

#include <cstddef>
#include <string>

std::string extract_Just_Filename(const std::string &inputString) {

    const size_t lastSlashLocation = inputString.find_last_of("/");

    //If no slash found, or if it is the last character just return the input
    //string
    if (lastSlashLocation == std::string::npos || lastSlashLocation == inputString.size() - 1) {
        return inputString;
    }
        //Otherwise return the part of the string after the last slash
    else {
        return inputString.substr(lastSlashLocation + 1);
    }
}
