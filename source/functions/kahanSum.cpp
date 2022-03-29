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

Function to take a std::vector of doubles, and sum its elements using Kahan
summation. This has a much lower error than a naive summation, which can be
extremely bad in some situations. Note this is much more expensive however, so
think carefully if you are tempted to use this to sum long vectors many times
over in a code.*/

#include <vector>
#include "kahanSum.hpp"

double kahanSum(const std::vector<double> &vec) {

    /* On the Kahan summation wiki page, compensation is called 'c',
    compensated value is called 'y', and temp is called 't'.*/
    double sum = 0;
    double compensation = 0;

    for (auto &element: vec){
        double compensatedValue = element - compensation;
        double temp = sum + compensatedValue;
        /* Note, Nick Maclaren recommends static_cast<double>(temp-sum) here,
        because some compilers may if you are very unlucky, actual do the wrong
        here, completely trashing the Kahan summation and essentially doing a
        naive sum instead. Something to bear in mind in case something very
        weird goes wrong, but I haven't put the cast in because it's probably
        an extremely unlikely situation, and nobody else seems to worry about
        that kind of thing.*/
        compensation = (temp - sum) - compensatedValue;
        sum = temp;
    }
    return sum;
}
