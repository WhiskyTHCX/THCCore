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


#include <cstdio>
#include <cstring>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void THC_RF_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS

    if(!refluxing) {
        CCTK_PARAMWARN("HRSCCore::refluxing is set to \"no\", but "
            "\"THC_Refluxing\" has been activated.");
    }

    int expected_nvars = -1;
    if(CCTK_Equals(eos_type, "barotropic") ||
       CCTK_Equals(eos_type, "ultrarelativistic")) {
        expected_nvars = 4;
    }
    else if(CCTK_Equals(eos_type, "ideal")) {
        expected_nvars = 5;
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        expected_nvars = 6;
    }
    else {
        CCTK_ERROR("This is a bug in THC_Refluxing");
    }

    if(expected_nvars != nvars) {
        char msg[512];
        snprintf(msg, 512, "THC_Refluxing::nvars should be %d, but is %d",
            expected_nvars, nvars);
        CCTK_PARAMWARN(msg);
    }
}
