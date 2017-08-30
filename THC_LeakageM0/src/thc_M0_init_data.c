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


#include <stdbool.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void THC_M0_InitData(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_M0_InitData");
    }

    int group_id = CCTK_GroupIndex("THC_LeakageM0::thc_leakage_vars");
    cGroupDynamicData group_data;
    int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
    assert(!ierr);
    assert(group_data.lsh[0] == nrad);
    assert(group_data.lbnd[0] == 0);

    size_t const lsiz = group_data.lsh[0]*group_data.lsh[1];

    *thc_leakage_M0_is_on = false;

    *thc_M0_nue_num_flux = 0;
    *thc_M0_nua_num_flux = 0;
    *thc_M0_nux_num_flux = 0;
    *thc_M0_nue_ene_flux = 0;
    *thc_M0_nua_ene_flux = 0;
    *thc_M0_nux_ene_flux = 0;

    memset(thc_M0_N_nue, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_N_nua, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_N_nux, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nue, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nua, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nux, 0, lsiz*sizeof(CCTK_REAL));

    memset(thc_M0_N_nue_old, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_N_nua_old, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_N_nux_old, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nue_old, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nua_old, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_E_nux_old, 0, lsiz*sizeof(CCTK_REAL));

    memset(thc_M0_ndens_nue, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_ndens_nua, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_ndens_nux, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_eave_nue, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_eave_nua, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_eave_nux, 0, lsiz*sizeof(CCTK_REAL));

    memset(thc_M0_abs_number, 0, lsiz*sizeof(CCTK_REAL));
    memset(thc_M0_abs_energy, 0, lsiz*sizeof(CCTK_REAL));

    memset(thc_M0_mask, 0, lsiz*sizeof(CCTK_INT));

    *thc_leakage_M0_time = cctk_time;
}
