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
#include <gsl/gsl_rng.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"
#define SQ(X) ((X)*(X))

void THC_ID_RandomEddies(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_RandomEddies");
    }

    int const siz = UTILS_GFSIZE(cctkGH);
    int const ksiz = pow(random_eddies_max_wave_number + 1, 3);

    CCTK_REAL * cux = malloc(2*ksiz*sizeof(*cux));
    CCTK_REAL * cuy = malloc(2*ksiz*sizeof(*cuy));
    CCTK_REAL * cuz = malloc(2*ksiz*sizeof(*cuz));

    CCTK_REAL * Pcux = malloc(2*ksiz*sizeof(*Pcux));
    CCTK_REAL * Pcuy = malloc(2*ksiz*sizeof(*Pcuy));
    CCTK_REAL * Pcuz = malloc(2*ksiz*sizeof(*Pcuz));

    gsl_rng * rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, 1);

    CCTK_REAL const strength = random_eddies_strength;

    cux[0] = cuy[0] = cuz[0] = 0;
    cux[1] = cuy[1] = cuz[1] = 0;
    for(int kk = 2; kk < 2*ksiz; ++kk) {
        cux[kk] = strength * (2*gsl_rng_uniform_pos(rng)-1.0) / ksiz;
        cuy[kk] = strength * (2*gsl_rng_uniform_pos(rng)-1.0) / ksiz;
        cuz[kk] = strength * (2*gsl_rng_uniform_pos(rng)-1.0) / ksiz;
    }

    gsl_rng_free(rng);

    CCTK_REAL const lambda = random_eddies_spectral_projector_lambda;
    for(int ki = 0; ki <= random_eddies_max_wave_number; ++ki)
    for(int kj = 0; kj <= random_eddies_max_wave_number; ++kj)
    for(int kk = 0; kk <= random_eddies_max_wave_number; ++kk) {
        int const kijk = 2*kk + 2*kj*(random_eddies_max_wave_number + 1) +
            2*ki*SQ(random_eddies_max_wave_number + 1);
        int const k2  = SQ(ki) + SQ(kj) + SQ(kk);
        if(k2 > 0) {
            CCTK_REAL const ik2 = 1.0/k2;

            Pcux[kijk] = lambda*cux[kijk] + ik2*(1 - 2*lambda)*(
                    ki*ki*cux[kijk] + ki*kj*cuy[kijk] + ki*kk*cuz[kijk]);
            Pcux[kijk+1] = lambda*cux[kijk+1] + ik2*(1 - 2*lambda)*(
                    ki*ki*cux[kijk+1] + ki*kj*cuy[kijk+1] + ki*kk*cuz[kijk+1]);
            Pcuy[kijk] = lambda*cuy[kijk] + ik2*(1 - 2*lambda)*(
                    kj*ki*cux[kijk] + kj*kj*cuy[kijk] + kj*kk*cuz[kijk]);
            Pcuy[kijk+1] = lambda*cuy[kijk+1] + ik2*(1 - 2*lambda)*(
                    kj*ki*cux[kijk+1] + kj*kj*cuy[kijk+1] + kj*kk*cuz[kijk+1]);
            Pcuz[kijk] = lambda*cuz[kijk] + ik2*(1 - 2*lambda)*(
                    kk*ki*cux[kijk] + kk*kj*cuy[kijk] + kk*kk*cuz[kijk]);
            Pcuz[kijk+1] = lambda*cuz[kijk+1] + ik2*(1 - 2*lambda)*(
                    kk*ki*cux[kijk+1] + kk*kj*cuy[kijk+1] + kk*kk*cuz[kijk+1]);
        }
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_random_eddies,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            rho[ijk] = 1.0;

            CCTK_REAL ux = 0;
            CCTK_REAL uy = 0;
            CCTK_REAL uz = 0;
            for(int ki = 0; ki <= random_eddies_max_wave_number; ++ki)
            for(int kj = 0; kj <= random_eddies_max_wave_number; ++kj)
            for(int kk = 0; kk <= random_eddies_max_wave_number; ++kk) {
                int const kijk = 2*kk + 2*kj*(random_eddies_max_wave_number + 1) +
                    2*ki*SQ(random_eddies_max_wave_number + 1);
                ux += Pcux[kijk]*cos(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk])) +
                    Pcux[kijk+1]*sin(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk]));
                uy += Pcuy[kijk]*cos(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk])) +
                    Pcuy[kijk+1]*sin(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk]));
                uz += Pcuz[kijk]*cos(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk])) +
                    Pcuz[kijk+1]*sin(2*M_PI*(ki*x[ijk]+kj*y[ijk]+kk*z[ijk]));
            }

            CCTK_REAL const W2 = 0.5*(1 + sqrt(1 + 4*(SQ(ux)+SQ(uy)+SQ(uz))));

            velx[ijk] = ux / W2;
            vely[ijk] = uy / W2;
            velz[ijk] = uz / W2;

            eps[ijk] = 1.0;

            for(int e = 0; e < ntracers; ++e) {
                int const kappa = random_eddies_tracer_wavenumber[e];
                tracer[ijk + e*siz] = sin(2*M_PI*kappa*x[ijk])*
                    sin(2*M_PI*kappa*y[ijk])*sin(2*M_PI*kappa*z[ijk]);
            }
        } UTILS_ENDLOOP3(thc_id_random_eddies);
    }

    free(cux);
    free(cuy);
    free(cuz);

    free(Pcux);
    free(Pcuy);
    free(Pcuz);
}
