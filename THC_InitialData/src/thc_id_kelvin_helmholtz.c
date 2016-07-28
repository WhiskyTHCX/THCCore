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

#define SQ(X) ((X)*(X))

void THC_ID_KelvinHelmholtz(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_KelvinHelmholtz");
    }

    CCTK_REAL const Gamma = 4.0 / 3.0;

    CCTK_REAL const rho0 = 0.5 + 0.005;
    CCTK_REAL const rho1 = 0.5 - 0.005;
    CCTK_REAL const Vshear = 0.5;

    CCTK_REAL const a = 0.01;
    CCTK_REAL const A0 = 0.1;
    CCTK_REAL const sigma = 0.1;

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_kelvin_helmholtz,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if(y[ijk] > 0) {
                rho[ijk]  = rho1 * tanh((y[ijk] - 0.5) / a) + rho0;
                velx[ijk] = Vshear * tanh((y[ijk] - 0.5) / a);
                vely[ijk] = A0 * Vshear * sin(2*M_PI*x[ijk]) *
                    exp(-SQ(y[ijk] - 0.5)/SQ(sigma));
                if(CCTK_Equals(kelvin_helmholtz_case, "2D")) {
                    velz[ijk] = 0.0;
                }
                else if(CCTK_Equals(kelvin_helmholtz_case, "3D")) {
                    velz[ijk] = 0.01*UTILS_RANDOM(0, 1);
                }
                eps[ijk]  = 1.0 / ((Gamma-1.0)*rho[ijk]);
            }
            else {
                rho[ijk]  = - rho1 * tanh((y[ijk] + 0.5) / a) + rho0;
                velx[ijk] = - Vshear * tanh((y[ijk] + 0.5) / a);
                vely[ijk] = - A0 * Vshear * sin(2*M_PI*x[ijk]) *
                    exp(-SQ(y[ijk] + 0.5)/SQ(sigma));
                if(CCTK_Equals(kelvin_helmholtz_case, "2D")) {
                    velz[ijk] = 0.0;
                }
                else if(CCTK_Equals(kelvin_helmholtz_case, "3D")) {
                    velz[ijk] = 0.01*UTILS_RANDOM(0, 1);
                }
                eps[ijk]  = 1.0 / ((Gamma-1.0)*rho[ijk]);
            }

            if(ntracers > 0) {
                tracer[ijk] = rho[ijk];
            }
        } UTILS_ENDLOOP3(thc_id_kelvin_helmholtz);
    }
}
