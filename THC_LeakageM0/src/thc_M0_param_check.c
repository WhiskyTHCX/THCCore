//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2015, David Radice <dradice@caltech.edu>
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


#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void THC_M0_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int ntheta = round(sqrt(nray/2));
    if(2*ntheta*ntheta != nray) {
        CCTK_PARAMWARN("nray must be equal to 2*ntheta*ntheta for some "
                "integer ntheta");
    }

    if(bns_sep_threshold > 0) {
        if(!CCTK_IsThornActive("BNSTracker") &&
                !CCTK_IsThornActive("BNSTrackerGen")) {
            CCTK_PARAMWARN("You need to have either BNSTracker or BNSTrackerGen "
                    "activated to use the bns_sep_threshold switch");
        }
    }
}
