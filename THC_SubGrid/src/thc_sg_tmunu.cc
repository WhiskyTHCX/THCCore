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
#include <cmath>

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

    for(int k = cctk_nghostzones[2]; k < cctk_lsh[2]-cctk_nghostzones[2]; ++k)
    for(int j = cctk_nghostzones[1]; j < cctk_lsh[1]-cctk_nghostzones[1]; ++j)
    for(int i = cctk_nghostzones[0]; i < cctk_lsh[0]-cctk_nghostzones[0]; ++i) {
        int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL const atmof = UTILS_BITMASK_CHECK_FLAG(bitmask[ijk],
                THC_FLAG_ATMOSPHERE) ? 0.0 : 1.0;
        assert(std::isfinite(tau_xx[ijk]));
        assert(std::isfinite(tau_xy[ijk]));
        assert(std::isfinite(tau_xz[ijk]));
        assert(std::isfinite(tau_yy[ijk]));
        assert(std::isfinite(tau_yz[ijk]));
        assert(std::isfinite(tau_zz[ijk]));
        eTxx[ijk] += atmof * tau_xx[ijk];
        eTxy[ijk] += atmof * tau_xy[ijk];
        eTxz[ijk] += atmof * tau_xz[ijk];
        eTyy[ijk] += atmof * tau_yy[ijk];
        eTyz[ijk] += atmof * tau_yz[ijk];
        eTzz[ijk] += atmof * tau_zz[ijk];
    }
}
