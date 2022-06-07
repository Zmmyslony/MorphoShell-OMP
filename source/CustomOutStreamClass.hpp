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

This is based on the code given in an answer to Stackoverflow 14155364.
See CustomOutStreamClass.cpp for function details.

Class to be used like std::cout << "text here" << std::endl etc, replacing
std::cout with an instance of CustomOutStreamClass that you've constructed:
someInstance.

The CustomOutStreamClass constructor used to construct someInstance should be
passed the path and name of a text file. Then using :
someInstance <<
instead of:
std::cout <<
will write whatever follows BOTH to std::cout AND to the chosen file. This is
useful if you want to write a log file but also want everything in the log
file to get written to std::cout e.g. for easy viewing in a terminal.
*/

#ifndef _CUSTOM_OUT_STREAM_CLASS_TAG_
#define CUSTOM_OUT_STREAM_CLASS_TAG_

#include <iostream>
#include <fstream>
#include <string>

class CustomOutStreamClass {

public:

    // Constructor (just default).
    CustomOutStreamClass() = default;

    /* Set path+name of file which will be printed to (e.g. a log file, in the
    principal use case). */
    void set_outputFileName(const std::string &);

    // Open or create output file corresponding to fileName member data.
    void open();

    // Close output file.
    void close();

    // Getter function to return a reference to the fileName.
    [[nodiscard]] const std::string &get_outputFileName() const;

    // Getter function to return a reference to the outputFileStream.
    [[nodiscard]] const std::ofstream &get_outputFileStream() const;

    /* For regular output of variables etc.
    This is defined here since the function is templated so defining it just in
    the .cpp file leads to problems with compilation/linking!*/
    template<typename T>
    CustomOutStreamClass &operator<<(const T &something) {

        std::cout << something;
        outputFileStream << something;

        return *this;
    }

    ///////////////////////////////////////////////////////////////////////////
    // For outputting manipulators such as std::endl.
    typedef std::ostream &(*streamFunction)(std::ostream &);

    CustomOutStreamClass &operator<<(streamFunction);
    //////////////////////////////////////////////////////////////////////////

private:
    // File path+name the file stream will correspond to.
    std::string fileName;

    // Actual output file stream
    std::ofstream outputFileStream;

};

#endif
