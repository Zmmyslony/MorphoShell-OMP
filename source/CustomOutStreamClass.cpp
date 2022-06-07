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

Non-templated member function definitions for CustomOutStreamClass.
See CustomOutStreamClass.hpp for details of the class.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "CustomOutStreamClass.hpp"


/* Set which file will be printed to (the log file, in the principal use
case). */
void CustomOutStreamClass::set_outputFileName(const std::string &outputFilePathAndName) {

    fileName = outputFilePathAndName;
}


// Open or create log file corresponding to fileName member data.
void CustomOutStreamClass::open() {

    // Open file in append mode.
    outputFileStream.open(fileName, std::ofstream::app);

    if (!outputFileStream) {
        throw std::runtime_error("Error: Problem creating or opening log file.");
    }
}

// Close log file
void CustomOutStreamClass::close() {
    outputFileStream.close();
}

// Getter function to return a reference to the fileName.
const std::string &CustomOutStreamClass::get_outputFileName() const {
    return fileName;
}

// Getter function to return a reference to the outputFileStream.
const std::ofstream &CustomOutStreamClass::get_outputFileStream() const {
    return outputFileStream;
}


///////////////////////////////////////////////////////////////////////////

// For outputting manipulators such as std::endl.
/*
How this works: we first for convenience define a function pointer typedef so
that streamFunction is an alias for: "a pointer to a functions that takes
a std::ostream& as input and returns a std::ostream&". std::endl is an
example of such a function. Then we overload the binary operator << to take
any such function, and use it as a SECOND argument. No first argument is
required because we're implementing this as a member function, so the first
argument is taken to be the calling object automatically!
Hence, some_my_ostream << std::endl takes some instance of my_ostream and feeds
that to its overloaded << operator, with the function (pointer I guess?)
std::endl as the second. What the << operator then does is apply the function
it took as an argument (std::endl here) and apply it in the way it usually
would be, first to std::cout, then to outputFileStream. std::endl takes a stream
(outputFileStream, say) and does its newline and flush thing, and then returns
a reference to the stream so that you can then chain, i.e. put another << or
similar straight after. We want this behaviour too, which is why our
overloaded << operator has a return type and returns a reference to the
calling object (some_my_ostream); so that chaining can occur.

See:
https://www.austincc.edu/comer/ds12s/binaryop.htm
Stackoverflow 4295432
Stackoverflow 3674200
Stackoverflow 14155364
*/

typedef std::ostream &(*streamFunction)(std::ostream &);

CustomOutStreamClass &CustomOutStreamClass::operator<<(const streamFunction func) {
    std::cout << std::scientific << std::setprecision(2);
    func(std::cout);
    func(outputFileStream);
    return *this;
}

////////////////////////////////////////////////////////////////////////////////
