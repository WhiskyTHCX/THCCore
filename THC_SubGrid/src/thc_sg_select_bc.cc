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


#include <cassert>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void THC_SG_SelectBC(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_SelectBC");
    }

    char const * mybc;
    if(CCTK_Equals(bc_type, "flat")) {
        mybc = "flat";
    }
    else {
        mybc = "none";
    }

    assert(cctk_nghostzones[0] == cctk_nghostzones[1]);
    assert(cctk_nghostzones[1] == cctk_nghostzones[2]);

    int ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES,
            cctk_nghostzones[0], -1, "THC_SubGrid::turb_stress", mybc);
    if(ierr) {
        CCTK_ERROR("Failed to select the BCs");
    }
}
