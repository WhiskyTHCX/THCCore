//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2016, David Radice <dradice@caltech.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#define SQ(X) ((X)*(X))

extern "C" void THC_SG_CalcVelDens(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_VelDens");
    }

    int const siz = UTILS_GFSIZE(cctkGH);
    for(int i = 0; i < siz; ++i) {
        veldens[i + 0*siz] = volform[i] * vel[i + 0*siz];
        veldens[i + 1*siz] = volform[i] * vel[i + 1*siz];
        veldens[i + 2*siz] = volform[i] * vel[i + 2*siz];
    }
}
