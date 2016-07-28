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


#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include "cctk.h"
#include "cctk_Parameters.h"

#include "thc_M0_kernel.h"

#include "utils.hh"

#define POW2(X) ((X)*(X))

#ifndef THC_LK_GRID_EPS
#define THC_LK_GRID_EPS 0.05
#endif

#define thc_M0_assert(e) \
    if(!e) { \
        std::fprintf(stderr, "%s:%u: failed assertion `%s'\n", __FILE__, \
                __LINE__, #e); \
        goto fail; \
    }

#define thc_M0_debug_print(var) \
    std::fprintf(stderr, "%s = %.19e\n", #var, var)

namespace {

void __rad_null_flat(CCTK_REAL const r, CCTK_REAL const theta,
        CCTK_REAL * kt, CCTK_REAL * kr, CCTK_REAL * chi,
        CCTK_REAL * sqrt_det_g) {
    *kt = 1.0;
    *kr = 1.0;
    *chi = 1.0;
    *sqrt_det_g = r*r*sin(theta);
}

}

extern "C" {

SphericalGrid * M0Grid = NULL;

void thc_M0_rad_null(
        int const irad,
        int const iray,
        CCTK_REAL const alp,
        CCTK_REAL const betax,
        CCTK_REAL const betay,
        CCTK_REAL const betaz,
        CCTK_REAL const gxx,
        CCTK_REAL const gxy,
        CCTK_REAL const gxz,
        CCTK_REAL const gyy,
        CCTK_REAL const gyz,
        CCTK_REAL const gzz,
        CCTK_REAL const zvecx,
        CCTK_REAL const zvecy,
        CCTK_REAL const zvecz,
        CCTK_INT  * mask,
        CCTK_REAL * kt,
        CCTK_REAL * kr,
        CCTK_REAL * chi,
        CCTK_REAL * sqrt_det_g) {
    /* Initialize some stuff */
    CCTK_REAL r, theta, phi;
    CCTK_REAL Jac[9];
    CCTK_REAL InvJac[9];
    CCTK_REAL gamma_c[9];
    CCTK_REAL gamma_s[9];
    CCTK_REAL shift_c[3] = {betax, betay, betaz};
    CCTK_REAL shift_s[3];
    CCTK_REAL g_s[16];
    CCTK_REAL zvec_c[3] = {zvecx, zvecy, zvecz};
    CCTK_REAL W2 = 1.0;
    CCTK_REAL w_lorentz;
    CCTK_REAL vel_c[3] = {0, 0, 0};
    CCTK_REAL vel_s[3] = {0, 0, 0};
    CCTK_REAL u_s[4];
    CCTK_REAL partial_r[4] = {0, 1, 0, 0};
    CCTK_REAL e_r[4];
    CCTK_REAL xi = 0;
    CCTK_REAL kappa[4];
    CCTK_REAL partial_t[4] = {1, 0, 0, 0};
    CCTK_REAL dir;

    /* Get coordinates and excise coordinate singularities */
    thc_sph_grid_get_r_theta_phi(M0Grid, irad, iray, &r, &theta, &phi);
    r = std::max(r, THC_LK_GRID_EPS*thc_sph_grid_get_dr(M0Grid));
    theta = std::max(theta, THC_LK_GRID_EPS*thc_sph_grid_get_dtheta(M0Grid));
    theta = std::min(theta, M_PI - THC_LK_GRID_EPS*
            thc_sph_grid_get_dtheta(M0Grid));

    /* Use flat spacetime quantities in the excision region */
    if(*mask) {
        __rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
        return;
    }

    /* Compute Jacobians */
    thc_coord_sph_to_cart_J(r, theta, phi, Jac);
    thc_coord_cart_to_sph_J(r, theta, phi, InvJac);

    /* Spatial metric in spherical coordinates */
    utils::metric::space(gxx, gxy, gxz, gyy, gyz, gzz, gamma_c);
    for(int a1 = 0; a1 < 3; ++a1)
    for(int b1 = 0; b1 < 3; ++b1) {
        gamma_s[a1*3 + b1] = 0;
        for(int a2 = 0; a2 < 3; ++a2)
        for(int b2 = 0; b2 < 3; ++b2) {
            gamma_s[a1*3 + b1] += Jac[a2*3 + a1] * Jac[b2*3 + b1] *
                                  gamma_c[a2*3 + b2];
        }
    }

    /* Transform shift vector */
    shift_c[0] = betax;
    shift_c[1] = betay;
    shift_c[2] = betaz;
    utils::gemv<CCTK_REAL, false, 3, 3>::simple(InvJac, shift_c, shift_s);

    /* Spacetime metric in spherical coordinates */
    utils::metric::spacetime(alp, shift_s[0], shift_s[1], shift_s[2],
            gamma_s[0], gamma_s[1], gamma_s[2], gamma_s[4], gamma_s[5],
            gamma_s[8], g_s);

    /* Transform 3-velocity and compute 4-velocity in sph. coordinates */
    zvec_c[0] = zvecx;
    zvec_c[1] = zvecy;
    zvec_c[2] = zvecz;
    W2 = utils::gdot<CCTK_REAL, 3>::simple(gamma_c, zvec_c, zvec_c) + 1;
    w_lorentz = std::sqrt(W2);
    /* First assert: all of the initialization must happen before this point */
    thc_M0_assert(std::isfinite(w_lorentz));

    vel_c[0] = zvec_c[0]/w_lorentz;
    vel_c[1] = zvec_c[1]/w_lorentz;
    vel_c[2] = zvec_c[2]/w_lorentz;
    utils::gemv<CCTK_REAL, false, 3, 3>::simple(InvJac, vel_c, vel_s);
    utils::valencia::uvel(alp, shift_s[0], shift_s[1], shift_s[2],
            w_lorentz, vel_s[0], vel_s[1], vel_s[2], u_s);

    /* Construct a radial vector ortogonal to the velocity */
    xi = utils::gdot<CCTK_REAL, 4>::simple(g_s, partial_r, u_s);
    for(int a = 0; a < 4; ++a) {
        e_r[a] = partial_r[a] + xi * u_s[a];
    }
    dir = utils::gdot<CCTK_REAL, 4>::simple(g_s, partial_r, e_r);
    xi = sqrt(utils::gdot<CCTK_REAL, 4>::simple(g_s, e_r, e_r));
    thc_M0_assert(std::isfinite(xi));
    for(int a = 0; a < 4; ++a) {
        e_r[a] = e_r[a]/xi * utils::sign(dir);
    }

    /* Construct the radial null vector from the diad u, e_r */
    kappa[0] = u_s[0] + e_r[0];
    kappa[1] = u_s[1] + e_r[1];
    kappa[2] = u_s[2] + e_r[2];
    kappa[3] = u_s[3] + e_r[3];
    *kt = kappa[0];
    *kr = kappa[1];
    thc_M0_assert((*kt > 0));

    /* Compute chi */
    *chi = - utils::gdot<CCTK_REAL, 4>::simple(g_s, kappa, partial_t);

    /* Handle forming horizons */
    if(*kr <= 0 || *chi <= 0) {
        char msg[512];
        snprintf(msg, 512, "Point (r,theta) = (%g, %g) inside the horizon "
                "but not excised!", r, theta);
#pragma omp critical
        CCTK_WARN(CCTK_WARN_COMPLAIN, msg);
        *mask = true;
        __rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
        return;
    }

    /* Compute sqrt_det_g */
    *sqrt_det_g = alp*std::sqrt(utils::metric::spatial_det(gamma_s));
    thc_M0_assert(std::isfinite(*sqrt_det_g));

    return;

fail:
    thc_M0_debug_print(r);
    thc_M0_debug_print(theta);
    thc_M0_debug_print(phi);
    thc_M0_debug_print(Jac[0]);
    thc_M0_debug_print(Jac[1]);
    thc_M0_debug_print(Jac[2]);
    thc_M0_debug_print(Jac[3]);
    thc_M0_debug_print(Jac[4]);
    thc_M0_debug_print(Jac[5]);
    thc_M0_debug_print(Jac[6]);
    thc_M0_debug_print(Jac[7]);
    thc_M0_debug_print(Jac[8]);
    thc_M0_debug_print(InvJac[0]);
    thc_M0_debug_print(InvJac[1]);
    thc_M0_debug_print(InvJac[2]);
    thc_M0_debug_print(InvJac[3]);
    thc_M0_debug_print(InvJac[4]);
    thc_M0_debug_print(InvJac[5]);
    thc_M0_debug_print(InvJac[6]);
    thc_M0_debug_print(InvJac[7]);
    thc_M0_debug_print(InvJac[8]);
    thc_M0_debug_print(gamma_c[0]);
    thc_M0_debug_print(gamma_c[1]);
    thc_M0_debug_print(gamma_c[2]);
    thc_M0_debug_print(gamma_c[3]);
    thc_M0_debug_print(gamma_c[4]);
    thc_M0_debug_print(gamma_c[5]);
    thc_M0_debug_print(gamma_c[6]);
    thc_M0_debug_print(gamma_c[7]);
    thc_M0_debug_print(gamma_c[8]);
    thc_M0_debug_print(gamma_s[0]);
    thc_M0_debug_print(gamma_s[1]);
    thc_M0_debug_print(gamma_s[2]);
    thc_M0_debug_print(gamma_s[3]);
    thc_M0_debug_print(gamma_s[4]);
    thc_M0_debug_print(gamma_s[5]);
    thc_M0_debug_print(gamma_s[6]);
    thc_M0_debug_print(gamma_s[7]);
    thc_M0_debug_print(gamma_s[8]);
    thc_M0_debug_print(alp);
    thc_M0_debug_print(shift_c[0]);
    thc_M0_debug_print(shift_c[1]);
    thc_M0_debug_print(shift_c[2]);
    thc_M0_debug_print(shift_s[0]);
    thc_M0_debug_print(shift_s[1]);
    thc_M0_debug_print(shift_s[2]);
    thc_M0_debug_print(g_s[0]);
    thc_M0_debug_print(g_s[1]);
    thc_M0_debug_print(g_s[2]);
    thc_M0_debug_print(g_s[3]);
    thc_M0_debug_print(g_s[4]);
    thc_M0_debug_print(g_s[5]);
    thc_M0_debug_print(g_s[6]);
    thc_M0_debug_print(g_s[7]);
    thc_M0_debug_print(g_s[8]);
    thc_M0_debug_print(g_s[9]);
    thc_M0_debug_print(g_s[10]);
    thc_M0_debug_print(g_s[11]);
    thc_M0_debug_print(g_s[12]);
    thc_M0_debug_print(g_s[13]);
    thc_M0_debug_print(g_s[14]);
    thc_M0_debug_print(g_s[15]);
    thc_M0_debug_print(w_lorentz);
    thc_M0_debug_print(vel_c[0]);
    thc_M0_debug_print(vel_c[1]);
    thc_M0_debug_print(vel_c[2]);
    thc_M0_debug_print(vel_s[0]);
    thc_M0_debug_print(vel_s[1]);
    thc_M0_debug_print(vel_s[2]);
    thc_M0_debug_print(u_s[0]);
    thc_M0_debug_print(u_s[1]);
    thc_M0_debug_print(u_s[2]);
    thc_M0_debug_print(u_s[3]);
    thc_M0_debug_print(e_r[0]);
    thc_M0_debug_print(e_r[1]);
    thc_M0_debug_print(e_r[2]);
    thc_M0_debug_print(e_r[3]);
    thc_M0_debug_print(xi);
    thc_M0_debug_print(kappa[0]);
    thc_M0_debug_print(kappa[1]);
    thc_M0_debug_print(kappa[2]);
    thc_M0_debug_print(kappa[3]);
    thc_M0_debug_print(*chi);
    thc_M0_debug_print(*sqrt_det_g);
    std::fflush(stderr);
    abort();
}

void thc_M0_evol_density(
        CCTK_REAL const dt,
        CCTK_INT  const * mask,
        CCTK_REAL const * theta,
        CCTK_REAL const * eta,
        CCTK_REAL const * mu,
        CCTK_REAL const * N_p,
        CCTK_REAL * N) {
    int const siz = thc_sph_grid_get_nrad(M0Grid);
    CCTK_REAL const dr = thc_sph_grid_get_dr(M0Grid);
    CCTK_REAL const lambda = dt/dr;

    /* Boundary condition: assume very optically thick */
    N[0] = 0;

    /* Rest of the domain */
    for(int i = 1; i < siz; ++i) {
        if(mask[i]) {
            N[i] = 0;
        }
        else {
            CCTK_REAL const a = 1 + lambda*theta[i] + dt*mu[i];
            CCTK_REAL const b = N_p[i] + lambda*theta[i-1]*N[i-1] + dt*eta[i];
            N[i] = b/a;
            if(!std::isfinite(N[i])) {
                std::fprintf(stderr, "%s:%u: N[i] is not finite!\n",
                        __FILE__, __LINE__);
                std::fprintf(stderr, "i = %d\n", i);
                thc_M0_debug_print(N[i]);
                thc_M0_debug_print(N[i-1]);
                thc_M0_debug_print(N_p[i]);
                thc_M0_debug_print(dr);
                thc_M0_debug_print(lambda);
                thc_M0_debug_print(theta[i]);
                thc_M0_debug_print(theta[i-1]);
                thc_M0_debug_print(mu[i]);
                thc_M0_debug_print(eta[i]);
                thc_M0_debug_print(a);
                thc_M0_debug_print(b);
                std::fflush(stderr);
                abort();
            }
        }
    }
}

void thc_M0_evol_energy_ave(
        CCTK_REAL const dt,
        CCTK_INT  const * mask,
        CCTK_REAL const * n,
        CCTK_REAL const * theta,
        CCTK_REAL const * eta,
        CCTK_REAL const * mu,
        CCTK_REAL const * E_p,
        CCTK_REAL * E) {
    int const siz = thc_sph_grid_get_nrad(M0Grid);
    CCTK_REAL const dr = thc_sph_grid_get_dr(M0Grid);
    CCTK_REAL const lambda = dt/dr;

    /* Boundary condition: assume equilibrium */
    if(!mask[0] && mu[0] > std::numeric_limits<double>::epsilon()) {
        E[0] = eta[0]/mu[0];
    }
    else {
        E[0] = 0;
    }

    /* Rest of the domain */
    for(int i = 1; i < siz; ++i) {
        if(mask[i]) {
            E[i] = 0;
        }
        else {
            CCTK_REAL const a = n[i]*(1.0 + lambda*theta[i]) + dt*mu[i];
            CCTK_REAL const b = n[i]*(E_p[i] + lambda*theta[i]*E[i-1]) +
                dt*eta[i];
            /* E is not well defined in vacuum and in absence of sources */
            if(std::abs(a) > std::numeric_limits<double>::epsilon()) {
                E[i] = b/a;
            }
            else {
                E[i] = E[i-1];
            }
            assert(std::isfinite(E[i]));
        }
    }
}

} // extern "C"
