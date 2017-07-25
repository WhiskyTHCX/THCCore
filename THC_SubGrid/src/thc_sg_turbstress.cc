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


#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "finite_difference.h"
#include "utils.hh"

#define SQ(X) ((X)*(X))

using namespace utils;

extern "C" void THC_SG_CalcSubgridTensor(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_SG_CalcSubgridTensor");
    }

    // Grid data
    int const gfsiz = UTILS_GFSIZE(cctkGH);
    CCTK_REAL const idelta[3] = {
        1.0/CCTK_DELTA_SPACE(0),
        1.0/CCTK_DELTA_SPACE(1),
        1.0/CCTK_DELTA_SPACE(2),
    };

    // Slicing geometry
    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);
    // Three velocity, contravariant
    tensor::generic<CCTK_REAL *, 3, 1> v_u;
    v_u[0] = &vel[0*gfsiz];
    v_u[1] = &vel[1*gfsiz];
    v_u[2] = &vel[2*gfsiz];
    // Subgrid stress tensor, covariant, undensitized
    tensor::symmetric2<CCTK_REAL *, 3, 2> tau_dd;
    tau_dd[0] = tau_xx;
    tau_dd[1] = tau_xy;
    tau_dd[2] = tau_xz;
    tau_dd[3] = tau_yy;
    tau_dd[4] = tau_yz;
    tau_dd[5] = tau_zz;
#pragma omp parallel
    {
        UTILS_LOOP3(thc_sg_tau_dd,
                k, cctk_nghostzones[2], cctk_lsh[2] - cctk_nghostzones[2],
                j, cctk_nghostzones[1], cctk_lsh[1] - cctk_nghostzones[1],
                i, cctk_nghostzones[0], cctk_lsh[0] - cctk_nghostzones[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Geometry on the slice
            tensor::metric<3> g_dd;
            geom.get_metric(ijk, &g_dd);
            tensor::inv_metric<3> g_uu;
            geom.get_inv_metric(ijk, &g_uu);

            // Metric derivatives
            utils::tensor::symmetric2<CCTK_REAL, 3, 3> dg_ddd;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
            for(int c = b; c < 3; ++c) {
                CCTK_REAL const * gbc = geom.get_space_metric_comp(b, c);
                dg_ddd(a,b,c) = idelta[a]*cdiff_1(cctkGH, gbc, i, j, k, a, 2);
            }

            // Christoffel symbols
            utils::tensor::symmetric2<CCTK_REAL, 3, 3> Gamma_udd;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
            for(int c = b; c < 3; ++c) {
                Gamma_udd(a,b,c) = 0.0;
                for(int d = 0; d < 3; ++d) {
                    Gamma_udd(a,b,c) += 0.5 * g_uu(a,d) *
                        (dg_ddd(c,d,b) + dg_ddd(b,c,d) - dg_ddd(d,b,c));
                }
            }

            // Gradient of the velocity, mixed components
            tensor::generic<CCTK_REAL, 3, 2> Dv_du;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b) {
                Dv_du(a,b) = idelta[a]*cdiff_1(cctkGH, v_u(b), i, j, k, a, 2);
                for(int c = 0; c < 3; ++c) {
                    Dv_du(a,b) += Gamma_udd(b,a,c) * v_u(c)[ijk];
                }
            }

            // Trace of the gradient of the velocity
            CCTK_REAL Tr_Dv = 0;
            for(int a = 0; a < 3; ++a) {
                Tr_Dv += Dv_du(a,a);
            }

            // Subgrid stress tensor, covariant components
            for(int a = 0; a < 3; ++a)
            for(int b = a; b < 3; ++b) {
                tau_dd(a,b)[ijk] = - 1.0/3.0 * Tr_Dv * g_dd(a,b);
                for(int c = 0; c < 3; ++c) {
                    tau_dd(a,b)[ijk] += 0.5*(
                        g_dd(c,b)*Dv_du(a,c) + g_dd(a,c)*Dv_du(b,c));
                }
                tau_dd(a,b)[ijk] *= (-2.0 * nu_turb[ijk]);
                tau_dd(a,b)[ijk] *= (rho[ijk] * (1.0 + eps[ijk]) +
                            press[ijk]) * SQ(w_lorentz[ijk]);
                assert(std::isfinite(tau_dd(a,b)[ijk]));
            }
        } UTILS_ENDLOOP3(thc_sg_tau_dd);
    }
}
