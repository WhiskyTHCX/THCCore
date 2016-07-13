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
    tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);
    // Densitized velocity, contravariant
    tensor::generic<CCTK_REAL *, 3, 1> v_u;
    v_u[0] = &veldens[0*gfsiz];
    v_u[1] = &veldens[1*gfsiz];
    v_u[2] = &veldens[2*gfsiz];
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
            tensor::metric<3> g_dd;
            geom.get_metric(ijk, &g_dd);

            // Gradient of the velocity, mixed components
            tensor::generic<CCTK_REAL, 3, 2> Dv_du;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b) {
                Dv_du(a,b) = 0.5 * idelta[a] *
                    (v_u(b)[ijk + stride[a]] - v_u(b)[ijk - stride[a]]);
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
                tau_dd(a,b)[ijk] *= (-2.0 * nu_turb[ijk] / volform[ijk]);
                tau_dd(a,b)[ijk] *= (rho[ijk] * (1.0 + eps[ijk]) +
                            press[ijk]) * SQ(w_lorentz[ijk]);
                assert(std::isfinite(tau_dd(a,b)[ijk]));
            }
        } UTILS_ENDLOOP3(thc_sg_tau_dd);
    }
}
