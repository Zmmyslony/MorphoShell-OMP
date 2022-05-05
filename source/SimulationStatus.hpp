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

Header file for enumeration of 'status', which refers to whether the
simulation is currently in a 'dialling in' phase, or is not dialling and is
instead waiting for equilibrium, or whether equilibrium has been reached but
the next dialling phase has not yet begun.*/

#ifndef _STATUS_ENUM_TAG_
#define _STATUS_ENUM_TAG_ 1

enum SimulationStatus {
    dialling, waitingForEquilibrium, equilibriumReached
};

#endif
