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


#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

#define POW2(X) ((X)*(X))

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

void THC_M0_TestInitHydro(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int const siz = UTILS_GFSIZE(cctkGH);
    CCTK_REAL * velx = &vel[0*siz];
    CCTK_REAL * vely = &vel[1*siz];
    CCTK_REAL * velz = &vel[2*siz];

    UTILS_LOOP3(thc_M0_test_init_hydro,
            k, 0, cctk_lsh[2],
            j, 0, cctk_lsh[1],
            i, 0, cctk_lsh[0]) {
        int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

        CCTK_REAL const rp = sqrt(POW2(x[ijk]) + POW2(y[ijk]) + POW2(z[ijk]));
        rho[ijk] = rp <= 1.0 ? 1.0 : 0.0;
        velx[ijk] = 0.0;
        vely[ijk] = 0.0;
        velz[ijk] = 0.0;
        eps[ijk]  = 1.0;
    } UTILS_ENDLOOP3(thc_M0_test_init_hydro);
}
