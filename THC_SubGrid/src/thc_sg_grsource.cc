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

#include "finite_difference.h"
#include "thc_macro.hh"
#include "utils.hh"

extern "C" void THC_SG_GRSource(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_GRSource");
    }

    // Grid data
    CCTK_REAL const idelta[3] = {
        1.0/CCTK_DELTA_SPACE(0),
        1.0/CCTK_DELTA_SPACE(1),
        1.0/CCTK_DELTA_SPACE(2),
    };

    // Slicing geometry
    utils::tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);
    // Subgrid stress tensor, covariant, undensitized
    utils::tensor::symmetric2<CCTK_REAL *, 3, 2> tau_dd;
    tau_dd[0] = tau_xx;
    tau_dd[1] = tau_xy;
    tau_dd[2] = tau_xz;
    tau_dd[3] = tau_yy;
    tau_dd[4] = tau_yz;
    tau_dd[5] = tau_zz;
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
    // Energy source
    CCTK_REAL * rhs_tau = static_cast<CCTK_REAL *>(
            utils::cctk::var_data_ptr(cctkGH, 0, "THC_Core::rhs_tau"));

#pragma omp parallel
    {
        UTILS_LOOP3(thc_sg_grsource,
                k, cctk_nghostzones[2], cctk_lsh[2]-cctk_nghostzones[2],
                j, cctk_nghostzones[1], cctk_lsh[1]-cctk_nghostzones[1],
                i, cctk_nghostzones[0], cctk_lsh[0]-cctk_nghostzones[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Geometry on the slice
            utils::tensor::metric<3> g_dd;
            geom.get_metric(ijk, &g_dd);
            utils::tensor::inv_metric<3> g_uu;
            geom.get_inv_metric(ijk, &g_uu);
            utils::tensor::symmetric2<CCTK_REAL, 3, 2> k_dd;
            geom.get_extr_curv(ijk, &k_dd);
            CCTK_REAL const vol4 = alp[ijk]*volform[ijk];

            // Contravariant components of the subgrid stress tensor
            utils::tensor::symmetric2<CCTK_REAL, 3, 2> tau_uu;
            for(int a = 0; a < 3; ++a)
            for(int b = a; b < 3; ++b) {
                tau_uu(a,b) = 0.0;
                for(int c = 0; c < 3; ++c)
                for(int d = 0; d < 3; ++d) {
                    tau_uu(a,b) += g_uu(a,c) * g_uu(d,b) * tau_dd(d,c)[ijk];
                }
            }

            // Metric derivatives
            utils::tensor::symmetric2<CCTK_REAL, 3, 3> dg_ddd;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
            for(int c = b; c < 3; ++c) {
                CCTK_REAL const * gbc = geom.get_space_metric_comp(b, c);
                dg_ddd(a,b,c) = idelta[a] *
                    cdiff_1(cctkGH, gbc, i, j, k, a, fd_order);
            }

            // Energy source term
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b) {
                rhs_tau[ijk] += vol4 * k_dd(a,b) * tau_uu(a,b);
            }

            // Momentum source terms
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
            for(int c = 0; c < 3; ++c) {
                dot_S_d(a)[ijk] += 0.5 * vol4 * tau_uu(b,c) * dg_ddd(a,b,c);
            }
        } UTILS_ENDLOOP3(thc_sg_grsource);
    }
}
