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


#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#define SQ(X) ((X)*(X))

extern "C" void THC_SG_RHS(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_RHS");
    }

    // Grid data
    int const gfsiz = UTILS_GFSIZE(cctkGH);
    int const stride[3] = {
        CCTK_GFINDEX3D(cctkGH, 1, 0, 0) - CCTK_GFINDEX3D(cctkGH, 0, 0, 0),
        CCTK_GFINDEX3D(cctkGH, 0, 1, 0) - CCTK_GFINDEX3D(cctkGH, 0, 0, 0),
        CCTK_GFINDEX3D(cctkGH, 0, 0, 1) - CCTK_GFINDEX3D(cctkGH, 0, 0, 0),
    };
    CCTK_REAL const idelta[3] = {
        1.0/CCTK_DELTA_SPACE(0),
        1.0/CCTK_DELTA_SPACE(1),
        1.0/CCTK_DELTA_SPACE(2),
    };

    // Slicing geometry
    utils::tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx,
            gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);
    // Conserved momentum, contravariant
    utils::tensor::generic<CCTK_REAL *, 3, 1> S_u;
    S_u[0] = &sconup[0*gfsiz];
    S_u[1] = &sconup[1*gfsiz];
    S_u[2] = &sconup[2*gfsiz];
    // Subgrid stress tensor, mixed, densitized and multiplied by the lapse
    utils::tensor::generic<CCTK_REAL *, 3, 2> tau_du;
    tau_du[0] = &scratch[0*gfsiz];
    tau_du[1] = &scratch[1*gfsiz];
    tau_du[2] = &scratch[2*gfsiz];
    tau_du[3] = &scratch[3*gfsiz];
    tau_du[4] = &scratch[4*gfsiz];
    tau_du[5] = &scratch[5*gfsiz];
    tau_du[6] = &scratch[6*gfsiz];
    tau_du[7] = &scratch[7*gfsiz];
    tau_du[8] = &scratch[8*gfsiz];
    // Momentum sources
    CCTK_REAL * rhs_sconx = static_cast<CCTK_REAL *>(
            utils::cctk::var_data_ptr(cctkGH, 0, "THC_Core::rhs_scon[0]"));
    CCTK_REAL * rhs_scony = static_cast<CCTK_REAL *>(
            utils::cctk::var_data_ptr(cctkGH, 0, "THC_Core::rhs_scon[1]"));
    CCTK_REAL * rhs_sconz = static_cast<CCTK_REAL *>(
            utils::cctk::var_data_ptr(cctkGH, 0, "THC_Core::rhs_scon[2]"));
    utils::tensor::generic<CCTK_REAL *, 3, 1> dot_S_d;
    dot_S_d[0] = rhs_sconx;
    dot_S_d[1] = rhs_scony;
    dot_S_d[2] = rhs_sconz;

#pragma omp parallel
    {
        for(int sign = -1; sign <= 1; sign += 2) {
            UTILS_LOOP3(thc_sg_tau_du,
                    k, cctk_nghostzones[2] - 1,
                       cctk_lsh[2] - cctk_nghostzones[2] + 1,
                    j, cctk_nghostzones[1] - 1,
                       cctk_lsh[1] - cctk_nghostzones[1] + 1,
                    i, cctk_nghostzones[0] - 1,
                       cctk_lsh[0] - cctk_nghostzones[0] + 1) {
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

                // Geometry on the slice
                utils::tensor::metric<3> g_dd;
                geom.get_metric(ijk, &g_dd);
                utils::tensor::inv_metric<3> g_uu;
                geom.get_inv_metric(ijk, &g_uu);

                // Gradient of the momentum, mixed components
                utils::tensor::generic<CCTK_REAL, 3, 2> DS_du;
                for(int a = 0; a < 3; ++a)
                for(int b = 0; b < 3; ++b) {
                    DS_du(a,b) = sign * idelta[0]*
                        (S_u(b)[ijk] - S_u(b)[ijk - sign*stride[a]]);
                }

                // Compute the subgrid stress tensor
                for(int a = 0; a < 3; ++a)
                for(int b = 0; b < 3; ++b) {
                    tau_du(a,b)[ijk] = DS_du(a,b);
                    for(int c = 0; c < 3; ++c)
                    for(int d = 0; d < 3; ++d) {
                        tau_du(a,b)[ijk] += g_dd(a,c) * g_uu(d,b) * DS_du(d,c);
                    }
                    tau_du(a,b)[ijk] *= (-nu_turb[ijk] * alp[ijk]);
                }
            } UTILS_ENDLOOP3(thc_sg_tau_du);
#pragma omp barrier
            UTILS_LOOP3(thc_sg_rhs,
                    k, cctk_nghostzones[2], cctk_lsh[2] - cctk_nghostzones[2],
                    j, cctk_nghostzones[1], cctk_lsh[1] - cctk_nghostzones[1],
                    i, cctk_nghostzones[0], cctk_lsh[0] - cctk_nghostzones[0]) {
                // Discrete divergence of the subgrid stress tensor
                int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
                for(int a = 0; a < 3; ++a)
                for(int b = 0; b < 3; ++b) {
                    dot_S_d(a)[ijk] -= 0.5 * sign * idelta[b] *
                        (tau_du(a,b)[ijk + sign*stride[a]] - tau_du(a,b)[ijk]);
                }
            } UTILS_ENDLOOP3(thc_sg_rhs);
        }
    }
}
