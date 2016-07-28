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

void THC_ID_Vortex(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_ID_Vortex");
    }

    CCTK_REAL const gamma_eps_inf = vortex_gamma * vortex_eps_inf;
    CCTK_REAL const press_inf = (vortex_gamma - 1.0)*vortex_eps_inf;
    CCTK_REAL const rcav = (1.0 + gamma_eps_inf)/
        sqrt(gamma_eps_inf*(gamma_eps_inf + 2));

    int const siz = UTILS_GFSIZE(cctkGH);

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_vortex,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            CCTK_REAL const rcyl = sqrt(POW2(x[ijk]) + POW2(y[ijk]));
            CCTK_REAL phi  = atan2(y[ijk], x[ijk]);
            if(phi < 0) {
                phi += 2*M_PI;
            }

            if(ntracers > 0) {
                tracer[ijk] = phi;
                if(ntracers > 1) {
                    tracer[ijk + siz] = rcyl;
                    if(ntracers > 2) {
                        tracer[ijk + 2*siz] = z[ijk];
                    }
                }
            }

            if(rcyl <= rcav) {
                rho[ijk]  = 0;
                velx[ijk] = 0;
                vely[ijk] = 0;
                velz[ijk] = 0;
                eps[ijk]  = 0;
            }
            else {
                CCTK_REAL const xi = 1.0/gamma_eps_inf*(
                        (1 + gamma_eps_inf)*sqrt(1 - POW2(1/rcyl)) - 1);
                assert(xi > 0);
                CCTK_REAL const u_phi2 = 1/(POW2(rcyl) - 1);
                CCTK_REAL v_phi  = sqrt(u_phi2/(1 + u_phi2));
                if(vortex_cutoff_radius > 0 && rcyl > vortex_cutoff_radius) {
                    v_phi = v_phi*exp(-(rcyl/vortex_cutoff_radius - 1));
                }

                rho[ijk]  = pow(xi, 1.0/(vortex_gamma - 1));
                velx[ijk] = - v_phi * sin(phi);
                vely[ijk] = v_phi * cos(phi);
                velz[ijk] = 0;
                eps[ijk]  = 1/(vortex_gamma - 1) * press_inf *
                    pow(rho[ijk], vortex_gamma - 1);
            }
        } UTILS_ENDLOOP3(thc_id_vortex);
    }
}
