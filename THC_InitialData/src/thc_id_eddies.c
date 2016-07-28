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


#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

void THC_ID_Eddies(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_Eddies");
    }

    int const siz = UTILS_GFSIZE(cctkGH);

    CCTK_REAL const vnorm = eddies_max_velocity/sqrt(3.0);
    CCTK_REAL const kappa = eddies_wave_number;

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_taylor_green,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            rho[ijk]  = 1.0;
            velx[ijk] = - vnorm*cos(2*M_PI*kappa*x[ijk])*
                sin(2*M_PI*kappa*y[ijk])*cos(2*M_PI*kappa*z[ijk]);
            vely[ijk] = vnorm*sin(2*M_PI*kappa*x[ijk])*
                cos(2*M_PI*kappa*y[ijk])*cos(2*M_PI*kappa*z[ijk]);
            velz[ijk] = - vnorm*sin(2*M_PI*kappa*x[ijk])*
                sin(2*M_PI*kappa*y[ijk])*sin(2*M_PI*kappa*z[ijk]);
            eps[ijk]  = 1.0;

            if(ntracers > 0) {
                for(int e = 0; e < ntracers; ++e) {
                    tracer[ijk + e*siz] = 0;
                }
            }
        } UTILS_ENDLOOP3(thc_id_taylor_green);
    }
}
