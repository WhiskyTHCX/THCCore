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


#ifndef THC_LEAKAGEM0_KERNEL_H
#define THC_LEAKAGEM0_KERNEL_H

#include "cctk.h"
#include "thc_sph_grid.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Spherical grid used by the M0 scheme */
extern SphericalGrid * M0Grid;

/* Indexing used to access the M0 arrays */
#define THC_M0_INDEX(group_data, irad, iray)                                  \
    ((irad - group_data.lbnd[0]) +                                            \
     group_data.ash[0]*(iray - group_data.lbnd[1]))

/*
 * Computes a null radial outgoing vector, normalized so that
 *   k_a u^a = -1
 * Where u^a is the four-velocity of the fluid.
 *
 * k^a is computed as
 *   k^a = u^a + r^a
 * where r^a is a unit-radial vector orthogonal to u^a
 *
 * This also computes the quantity
 *   \chi = - k_a t^a
 * Which is useful to compute the redshift of photons traveling along the null
 * ray, since their energy E in the fluid rest-frame satisfies
 *   k^a \partial_a (\chi E) = \chi Q - \sigma_abs E
 * If t^a is a Killing vector.
 *
 * Finally we also compute \sqrt{det(g)} (in spherical coordinates!), as well
 * as the proper volume of the two sphere Vol2
 */
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
        CCTK_REAL * sqrt_det_g);

/*
 * Solves the radial continuity equation
 *   \partial_t N + \partial_r [ N \theta ] = \eta - \mu N
 * with \theta, \eta and \mu >= 0, using a 1st order implicit upwind
 * method
 */
void thc_M0_evol_density(
        CCTK_REAL const dt,
        CCTK_INT  const * mask,
        CCTK_REAL const * theta,
        CCTK_REAL const * eta,
        CCTK_REAL const * mu,
        CCTK_REAL const * N_p,
        CCTK_REAL * N);

/*
 * Solves the radial advection problem
 *    n [ \partial_t E + \theta \partial_r E ] = \eta - \mu E
 * with \theta, \eta and \mu >= 0, using a 1st order implicit upwind
 * method.
 */
void thc_M0_evol_energy_ave(
        CCTK_REAL const dt,
        CCTK_INT  const * mask,
        CCTK_REAL const * n,
        CCTK_REAL const * theta,
        CCTK_REAL const * eta,
        CCTK_REAL const * mu,
        CCTK_REAL const * E_p,
        CCTK_REAL * E);

#ifdef __cplusplus
}
#endif

#endif
