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


#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_M0_kernel.h"

void THC_M0_SetupGrid(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_M0_SetupGrid");
    }

    int ntheta = round(sqrt(nray/2));
    int nphi = 2*ntheta;
    assert(ntheta*nphi == nray);

    CCTK_REAL const zero[3] = {0,0,0};
    thc_sph_grid_init(&M0Grid, zero, rmax, nrad, ntheta, nphi, true);

    int group_id = CCTK_GroupIndex("THC_LeakageM0::thc_leakage_vars");
    cGroupDynamicData group_data;
    int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
    assert(!ierr);
    assert(group_data.lsh[0] == nrad);
    assert(group_data.lbnd[0] == 0);

    for(int iray = group_data.lbnd[1]; iray <= group_data.ubnd[1]; ++iray)
    for(int irad = group_data.lbnd[0]; irad <= group_data.ubnd[0]; ++irad) {
        int const ij = THC_M0_INDEX(group_data, irad, iray);
        thc_sph_grid_get_x_y_z(M0Grid, irad, iray, &thc_M0_x[ij], &thc_M0_y[ij],
                &thc_M0_z[ij]);
    }
}

void THC_M0_RecoverGrid(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_M0_RecoverGrid");
    }

    int ntheta = round(sqrt(nray/2));
    int nphi = 2*ntheta;
    assert(ntheta*nphi == nray);

    CCTK_REAL const zero[3] = {0,0,0};
    thc_sph_grid_init(&M0Grid, zero, rmax, nrad, ntheta, nphi, true);
}

void THC_M0_FreeGrid(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS
    if(verbose) {
        CCTK_INFO("THC_M0_FreeGrid");
    }
    thc_sph_grid_free(M0Grid);
    M0Grid = NULL;
}
