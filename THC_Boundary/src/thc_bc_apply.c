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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

#define BC_NONE         0
#define BC_FLAT         1
#define BC_NOSLIP       2
#define BC_WALL         3


void THC_BC_Apply(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_BC_Apply");
    }

    bool const eos_barotropic = CCTK_Equals(eos_type, "barotropic");
    bool const eos_ideal      = CCTK_Equals(eos_type, "ideal");
    bool const eos_nuclear    = CCTK_Equals(eos_type, "nuclear");

#define SELECT_BC(BND) \
    int bc_ ## BND ## _t; \
    CCTK_REAL sign_ ## BND ## _n, sign_ ## BND ## _t; \
    if(CCTK_Equals(bc_ ## BND, "none")) { \
        bc_ ## BND ## _t = BC_NONE; \
        sign_ ## BND ## _n = 1.0; \
        sign_ ## BND ## _t = 1.0; \
    } \
    else if(CCTK_Equals(bc_ ## BND, "flat")) { \
        bc_ ## BND ## _t   = BC_FLAT; \
        sign_ ## BND ## _n = 1.0; \
        sign_ ## BND ## _t = 1.0; \
    } \
    else if(CCTK_Equals(bc_ ## BND, "noslip")) { \
        bc_ ## BND ## _t   = BC_NOSLIP; \
        sign_ ## BND ## _n = -1.0; \
        sign_ ## BND ## _t = -1.0; \
    } \
    else if(CCTK_Equals(bc_ ## BND, "wall")) { \
        bc_ ## BND ## _t   = BC_WALL; \
        sign_ ## BND ## _n = -1.0; \
        sign_ ## BND ## _t = 1.0; \
    } \
    else { \
        CCTK_ERROR("Unknown BC type"); \
    }

    SELECT_BC(xmin);
    SELECT_BC(xmax);
    SELECT_BC(ymin);
    SELECT_BC(ymax);
    SELECT_BC(zmin);
    SELECT_BC(zmax);

#undef SELECT_BC

    CCTK_REAL * sconx = &scon[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * scony = &scon[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * sconz = &scon[2*UTILS_GFSIZE(cctkGH)];

    if(cctk_bbox[0] && bc_xmin_t != BC_NONE) {
        for(int k = 0; k < cctk_lsh[2]; ++k)
        for(int j = 0; j < cctk_lsh[1]; ++j)
        for(int i = 0; i < cctk_nghostzones[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_xmin_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, cctk_nghostzones[0], j, k);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                        2*cctk_nghostzones[0] - i - 1, j, k);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_xmin_n * sconx[ijk_r];
            scony[ijk_p] = sign_xmin_t * scony[ijk_r];
            sconz[ijk_p] = sign_xmin_t * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p] = tau[ijk_r];
            }
        }
    }

    if(cctk_bbox[1] && bc_xmax_t != BC_NONE) {
        for(int k = 0; k < cctk_lsh[2]; ++k)
        for(int j = 0; j < cctk_lsh[1]; ++j)
        for(int i = cctk_lsh[0] - cctk_nghostzones[0]; i < cctk_lsh[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_xmax_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, cctk_lsh[0] -
                        cctk_nghostzones[0] - 1, j, k);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                    2*(cctk_lsh[0] - cctk_nghostzones[0]) - i - 1, j, k);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_xmax_n * sconx[ijk_r];
            scony[ijk_p] = sign_xmax_t * scony[ijk_r];
            sconz[ijk_p] = sign_xmax_t * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p]   = tau[ijk_r];
            }
        }
    }

    if(cctk_bbox[2] && bc_ymin_t != BC_NONE) {
        for(int k = 0; k < cctk_lsh[2]; ++k)
        for(int j = 0; j < cctk_nghostzones[1]; ++j)
        for(int i = 0; i < cctk_lsh[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_ymin_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, i, cctk_nghostzones[1], k);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                    i, 2*cctk_nghostzones[1] - j - 1, k);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_ymin_t * sconx[ijk_r];
            scony[ijk_p] = sign_ymin_n * scony[ijk_r];
            sconz[ijk_p] = sign_ymin_t * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p] = tau[ijk_r];
            }
        }
    }

    if(cctk_bbox[3] && bc_ymax_t != BC_NONE) {
        for(int k = 0; k < cctk_lsh[2]; ++k)
        for(int j = cctk_lsh[1] - cctk_nghostzones[1]; j < cctk_lsh[1]; ++j)
        for(int i = 0; i < cctk_lsh[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_ymax_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, i, cctk_lsh[1] -
                        cctk_nghostzones[1] - 1, k);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                    i, 2*(cctk_lsh[1] - cctk_nghostzones[1]) - j - 1, k);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_ymax_t * sconx[ijk_r];
            scony[ijk_p] = sign_ymax_n * scony[ijk_r];
            sconz[ijk_p] = sign_ymax_t * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p] = tau[ijk_r];
            }
        }
    }

    if(cctk_bbox[4] && bc_zmin_t != BC_NONE) {
        for(int k = 0; k < cctk_nghostzones[2]; ++k)
        for(int j = 0; j < cctk_lsh[1]; ++j)
        for(int i = 0; i < cctk_lsh[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_zmin_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, i, j, cctk_nghostzones[2]);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                    i, j, 2*cctk_nghostzones[2] - k - 1);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_zmin_t * sconx[ijk_r];
            scony[ijk_p] = sign_zmin_t * scony[ijk_r];
            sconz[ijk_p] = sign_zmin_n * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p] = tau[ijk_r];
            }
        }
    }

    if(cctk_bbox[5] && bc_zmax_t != BC_NONE) {
        for(int k = cctk_lsh[2] - cctk_nghostzones[2]; k < cctk_lsh[2]; ++k)
        for(int j = 0; j < cctk_lsh[1]; ++j)
        for(int i = 0; i < cctk_lsh[0]; ++i) {
            int const ijk_p = CCTK_GFINDEX3D(cctkGH, i, j, k);
            int ijk_r;
            if(bc_zmax_t == BC_FLAT) {
                ijk_r = CCTK_GFINDEX3D(cctkGH, i, j, cctk_lsh[2] -
                        cctk_nghostzones[2] - 1);
            }
            else {
                ijk_r = CCTK_GFINDEX3D(cctkGH,
                    i, j, 2*(cctk_lsh[2] - cctk_nghostzones[2]) - k - 1);
            }

            if(eos_barotropic || eos_ideal) {
                dens[ijk_p]  = dens[ijk_r];
            }
            if(eos_nuclear) {
                densxn[ijk_p] = densxn[ijk_r];
                densxp[ijk_p] = densxp[ijk_r];
            }

            sconx[ijk_p] = sign_zmax_t * sconx[ijk_r];
            scony[ijk_p] = sign_zmax_t * scony[ijk_r];
            sconz[ijk_p] = sign_zmax_n * sconz[ijk_r];

            if(!eos_barotropic) {
                tau[ijk_p] = tau[ijk_r];
            }
        }
    }
}
