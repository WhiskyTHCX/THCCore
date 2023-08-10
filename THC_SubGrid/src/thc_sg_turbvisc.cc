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


#include <algorithm>
#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

static CCTK_REAL kiuchi_lmix_fit(
        CCTK_REAL const rho,
        CCTK_REAL const stretch = 1.0,
        CCTK_REAL const ampl = 1.0) {
    // Fitting coefficients
    CCTK_REAL const lrho_0 = -8.49697276793;
    CCTK_REAL const a      = 0.0151145023166;
    CCTK_REAL const b      = -0.425383267966;

    CCTK_REAL const xi = (std::log10(rho) - lrho_0)/stretch;
    if(xi < 0) {
        return 0.0;
    }
    else {
        return a*ampl*xi*std::exp(-std::pow(std::abs(xi*b), 2.5));
    }
}

static CCTK_REAL kiuchi2_lmix_fit(
        CCTK_REAL const rho,
        CCTK_REAL const lrho0,
        CCTK_REAL const lrho1,
        CCTK_REAL const a,
        CCTK_REAL const b,
        CCTK_REAL const c) {
    CCTK_REAL const lrho = std::max(std::log10(rho), lrho0);
    assert(std::isfinite(lrho));
    if (lrho <= lrho1) {
        return std::pow(10.0, a + b*lrho);
    }
    else {
        return std::pow(10.0, a + b*lrho1 + c*(lrho - lrho1));
    }
}

extern "C" void THC_SG_CalcTurbVisc(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_CalcTurbVisc");
    }

    int const siz = UTILS_GFSIZE(cctkGH);

    if(CCTK_Equals(viscosity, "alpha")) {
        for(int i = 0; i < siz; ++i) {
            nu_turb[i] = lmix * csound[i];
            assert(std::isfinite(nu_turb[i]));
        }
    }
    else if(CCTK_Equals(viscosity, "const")) {
        for(int i = 0; i < siz; ++i) {
            nu_turb[i] = lmix;
            assert(std::isfinite(nu_turb[i]));
        }
    }
    else if(CCTK_Equals(viscosity, "kiuchi")) {
#pragma omp parallel
        {
            UTILS_LOOP3(thc_sg_kiuchi,
                    k, 0, cctk_lsh[2],
                    j, 0, cctk_lsh[1],
                    i, 0, cctk_lsh[0]) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                CCTK_REAL const lmix = kiuchi_lmix_fit(rho[ijk],
                        kiuchi_stretch, kiuchi_ampl);
                nu_turb[ijk] = lmix * csound[ijk];
                assert(std::isfinite(lmix));
                assert(std::isfinite(nu_turb[ijk]));
            } UTILS_ENDLOOP3(thc_sg_kiuchi);
        }
    }
    else if(CCTK_Equals(viscosity, "kiuchi2")) {
#pragma omp parallel
        {
            UTILS_LOOP3(thc_sg_kiuchi2,
                    k, 0, cctk_lsh[2],
                    j, 0, cctk_lsh[1],
                    i, 0, cctk_lsh[0]) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                CCTK_REAL const lmix = kiuchi2_lmix_fit(rho[ijk],
                        kiuchi2_lrho0, kiuchi2_lrho1,
                        kiuchi2_a, kiuchi2_b, kiuchi2_c);
                nu_turb[ijk] = lmix * csound[ijk];
                assert(std::isfinite(lmix));
                assert(std::isfinite(nu_turb[ijk]));
            } UTILS_ENDLOOP3(thc_sg_kiuchi2);
        }
    }
    else {
        CCTK_ERROR("Unknown viscosity prescription");
    }

#pragma omp parallel
    {
        UTILS_LOOP3(thc_sg_lapse_cut,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            if (alp[ijk] < lapse_cut) {
                nu_turb[ijk] = 0.0;
            }
        } UTILS_ENDLOOP3(thc_sg_lapse_cut);
    }

    if(nu_turb_limit >= 0.0) {
        CCTK_REAL const dx = CCTK_DELTA_SPACE(0);
        CCTK_REAL const dt = CCTK_DELTA_TIME;
#pragma omp parallel
        {
            UTILS_LOOP3(thc_sg_limit,
                    k, 0, cctk_lsh[2],
                    j, 0, cctk_lsh[1],
                    i, 0, cctk_lsh[0]) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                nu_turb[ijk] = std::min(nu_turb_limit*(dx/dt)*dx, nu_turb[ijk]);
            } UTILS_ENDLOOP3(thc_sg_limit);
        }
    }
}
