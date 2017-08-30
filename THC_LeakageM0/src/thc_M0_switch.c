//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2014, David Radice <dradice@caltech.edu>
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


#include <assert.h>
#include <stdbool.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void THC_M0_Switch(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if((cctk_iteration-1) % compute_every != 0) {
        return;
    }
    if(verbose) {
        CCTK_INFO("THC_M0_Switch");
    }

    if(wait_until_time > 0 && cctk_time < wait_until_time) {
        return;
    }

    bool thc_leakage_M0_was_on = *thc_leakage_M0_is_on;

    if(bns_sep_threshold > 0) {
        assert(CCTK_IsThornActive("BNSTrackerGen"));
        CCTK_REAL const * bns_sep_tot = CCTK_VarDataPtr(cctkGH, 0,
                "BNSTrackerGen::bns_sep_tot");
        assert(bns_sep_tot);

        if(*bns_sep_tot < bns_sep_threshold) {
            *thc_leakage_M0_is_on = true;
        }
        else {
            *thc_leakage_M0_is_on = false;
        }
    }
    else {
        *thc_leakage_M0_is_on = true;
    }

    if((!thc_leakage_M0_was_on && *thc_leakage_M0_is_on)) {
        *thc_leakage_M0_time = cctk_time;
        CCTK_INFO("THC_M0_Switch: leakage is on");
    }
    else if(thc_leakage_M0_was_on && !*thc_leakage_M0_is_on) {
        /* Cleanup everything if we are switching the M0 off */
        THC_M0_InitData(CCTK_PASS_CTOC);
        CCTK_INFO("THC_M0_Switch: leakage is off");
    }
    else if(verbose) {
        if(*thc_leakage_M0_is_on) {
            CCTK_INFO("THC_M0_Switch: leakage is on");
        }
        else {
            CCTK_INFO("THC_M0_Switch: leakage is off");
        }
    }
}
