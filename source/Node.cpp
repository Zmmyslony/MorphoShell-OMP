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

This file defines the member functions for the Node class that holds
the data for each node. The class constructors are the exception; they are
left in the header file for clarity there.*/

#include <iostream>
#include <iomanip>

#include "Node.hpp"

//This is a debugging tool to display the node's data
void Node::display() {
    nodeLogStream.open();
    nodeLogStream << "-----------------------------" << std::setprecision(15) << std::boolalpha << std::endl;
    nodeLogStream << "Node " << label << ":" << std::endl;
    nodeLogStream << "Labels of incident triangles: " << "\n" << incidentTriLabels << std::endl;
    nodeLogStream << "neighbourNodeLabels = " << "\n" << neighbourNodeLabels << std::endl;
    nodeLogStream << "Position = " << "\n" << pos << std::endl;
    nodeLogStream << "Velocity = " << "\n" << vel << std::endl;
    nodeLogStream << "Force = " << "\n" << force << std::endl;
    nodeLogStream << "Mass = " << mass << std::endl;
    nodeLogStream << "Boundary indicator: " << isOnBoundary << std::endl;
    nodeLogStream << "Clamp indicator: " << isClamped << std::endl;
    nodeLogStream << "Load indicator: " << isLoadForceEnabled << std::endl;
    nodeLogStream << "-----------------------------" << std::endl;
    nodeLogStream.close();
}
