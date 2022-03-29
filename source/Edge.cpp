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

This file defines the member functions for the Edge class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there*/

#include <iostream>

#include "Edge.hpp"


//This is a debugging tool to display the edge's data
void Edge::display() {
    edgeLogStream.open();
    edgeLogStream << "-----------------------------" << std::boolalpha << std::endl;
    edgeLogStream << "Edge " << label << ":" << std::endl;
    edgeLogStream << "Node labels " << nodeLabels.transpose() << std::endl;
    edgeLogStream << "Adjacent triangle labels: " << adjTriLabels.transpose() << std::endl;
    //edgeLogStream << "Left and Right 'other' (non-edge) node labels: " << otherNodeLabel_L << ", " << otherNodeLabel_R << std::endl;
    edgeLogStream << "Boundary indicator: " << isOnBoundary << std::endl;
    edgeLogStream << "-----------------------------" << std::endl;
    edgeLogStream.close();
}
