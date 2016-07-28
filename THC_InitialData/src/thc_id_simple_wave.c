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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

static double soundspeed(double kappa, double gamma, double rho) {
    double const eps = (kappa*pow(rho, gamma - 1)) / (gamma - 1.0);
    double const enthalpy = 1.0 + eps + kappa*pow(rho, gamma - 1);
    return sqrt(kappa*gamma*pow(rho, gamma-1)/enthalpy);
}

static double riemann_m(double kappa, double gamma,
        double rho, double vel) {
    double const cs = soundspeed(kappa, gamma, rho);
    double const tmp = sqrt(gamma - 1.0);
    return 0.5*log((1.0 + vel)/(1.0 - vel)) -
            1.0/tmp*log((tmp + cs)/(tmp - cs));
}

void THC_ID_SimpleWave(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_SimpleWave");
    }

    CCTK_REAL const Jm = riemann_m(simple_wave_kappa,
            simple_wave_gamma, 1.0, 0.0);
    CCTK_REAL const tmp = sqrt(simple_wave_gamma - 1.0);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_simple_wave,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if(fabs(x[ijk]) < 0.3) {
//                rho[ijk] = 1.0 + pow(pow(x[ijk]/0.3, 2) - 1.0, 4);
                rho[ijk] = 1.0 + exp(-1.0/(1.0 - pow(x[ijk]/0.3, 2)));
            }
            else {
                rho[ijk] = 1.0;
            }

            CCTK_REAL const cs = soundspeed(simple_wave_kappa,
                    simple_wave_gamma, rho[ijk]);
            CCTK_REAL const frho = Jm + 1.0/tmp*log((tmp + cs)/(tmp - cs));
            CCTK_REAL const efrho = pow(exp(frho), 2);
            velx[ijk] = (efrho - 1.0)/(efrho + 1.0);
            vely[ijk] = 0;
            velz[ijk] = 0;

            eps[ijk] = simple_wave_kappa*pow(rho[ijk],
                    simple_wave_gamma - 1.0) / (simple_wave_gamma - 1.0);
        } UTILS_ENDLOOP3(thc_id_simple_wave);
    }
}
