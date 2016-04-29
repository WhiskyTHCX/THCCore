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

#include "thc_macro.hh"
#include "utils.hh"

extern "C" void THC_SG_AddToTmunu(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_AddToTmunu");
    }

    CCTK_INT const * bitmask = static_cast<CCTK_INT *>(
            utils::cctk::var_data_ptr(cctkGH, 0, "THC_Core::bitmask"));

    int const gfsiz = UTILS_GFSIZE(cctkGH);
    for(int i = 0; i < gfsiz; ++i) {
        CCTK_REAL const atmof = UTILS_BITMASK_CHECK_FLAG(bitmask[i],
                THC_FLAG_ATMOSPHERE) ? 0.0 : 1.0;
        eTxx[i] += atmof * tau_xx[i];
        eTxy[i] += atmof * tau_xy[i];
        eTxz[i] += atmof * tau_xz[i];
        eTyy[i] += atmof * tau_yy[i];
        eTyz[i] += atmof * tau_yz[i];
        eTzz[i] += atmof * tau_zz[i];
    }
}
