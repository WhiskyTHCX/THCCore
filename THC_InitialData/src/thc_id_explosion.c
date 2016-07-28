//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
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

#include "utils_macro.h"

void THC_ID_Explosion(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_Explosion");
    }

    CCTK_REAL const Gamma = 1.4;

    CCTK_REAL const rho_1 = 1.0;
    CCTK_REAL const press_1 = 1.0;
    CCTK_REAL const eps_1 = press_1 / ((Gamma-1.0)*rho_1);

    CCTK_REAL const rho_2 = 0.125;
    CCTK_REAL const press_2 = 0.1;
    CCTK_REAL const eps_2 = press_2 / ((Gamma-1.0)*rho_2);


#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_explosion,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if(r[ijk] < 0.4) {
                rho[ijk] = rho_1;
                eps[ijk] = eps_1;
            }
            else {
                rho[ijk] = rho_2;
                eps[ijk] = eps_2;
            }
            velx[ijk] = 0;
            vely[ijk] = 0;
            velz[ijk] = 0;
        } UTILS_ENDLOOP3(thc_id_explosion);
    }
}
