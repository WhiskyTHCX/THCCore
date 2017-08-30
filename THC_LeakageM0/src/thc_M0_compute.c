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


#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "thc_lk_rates.h"
#include "thc_M0_kernel.h"

#define SQ(X) ((X)*(X))

void THC_M0_Compute(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if((cctk_iteration-1) % compute_every != 0) {
        return;
    }
    if(!*thc_leakage_M0_is_on) {
        return;
    }

    if(verbose) {
        CCTK_INFO("THC_M0_Compute");
    }

    CCTK_REAL const dt = cctk_time - (*thc_leakage_M0_time);
    assert(dt >= 0);
    CCTK_REAL const dtheta = thc_sph_grid_get_dtheta(M0Grid);
    CCTK_REAL const dphi = thc_sph_grid_get_dphi(M0Grid);

    int group_id = CCTK_GroupIndex("THC_LeakageM0::thc_leakage_vars");
    cGroupDynamicData group_data;
    int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
    assert(!ierr);
    assert(group_data.lsh[0] == nrad);
    assert(group_data.lbnd[0] == 0);

    CCTK_REAL const mb = AverageBaryonMass();
    assert(isfinite(mb));

    /* Initialize fluxes variables */
    CCTK_REAL my_M0_nue_num_flux = 0;
    CCTK_REAL my_M0_nua_num_flux = 0;
    CCTK_REAL my_M0_nux_num_flux = 0;
    CCTK_REAL my_M0_nue_ene_flux = 0;
    CCTK_REAL my_M0_nua_ene_flux = 0;
    CCTK_REAL my_M0_nux_ene_flux = 0;

    /* Rotate timelevels if using M0 */
    int const lsiz = group_data.lsh[0]*group_data.lsh[1];
    memcpy(thc_M0_N_nue_old, thc_M0_N_nue, lsiz*sizeof(CCTK_REAL));
    memcpy(thc_M0_N_nua_old, thc_M0_N_nua, lsiz*sizeof(CCTK_REAL));
    memcpy(thc_M0_N_nux_old, thc_M0_N_nux, lsiz*sizeof(CCTK_REAL));
    memcpy(thc_M0_E_nue_old, thc_M0_E_nue, lsiz*sizeof(CCTK_REAL));
    memcpy(thc_M0_E_nua_old, thc_M0_E_nua, lsiz*sizeof(CCTK_REAL));
    memcpy(thc_M0_E_nux_old, thc_M0_E_nux, lsiz*sizeof(CCTK_REAL));

    /* Compute everything ray-by-ray */
#pragma omp parallel
    {
        CCTK_REAL * kt           = malloc(nrad*sizeof(CCTK_REAL));
        CCTK_REAL * kr           = malloc(nrad*sizeof(CCTK_REAL));
        CCTK_REAL * chi          = malloc(nrad*sizeof(CCTK_REAL));
        CCTK_REAL * sqrt_det_g   = malloc(nrad*sizeof(CCTK_REAL));

        CCTK_REAL * theta = malloc(nrad*sizeof(CCTK_REAL));
        CCTK_REAL * eta   = malloc(nrad*sizeof(CCTK_REAL));
        CCTK_REAL * mu    = malloc(nrad*sizeof(CCTK_REAL));

#pragma omp for reduction(+ : my_M0_nue_num_flux, my_M0_nua_num_flux, \
        my_M0_nux_num_flux, my_M0_nue_ene_flux, my_M0_nua_ene_flux, \
        my_M0_nux_ene_flux)
        for(int iray = group_data.lbnd[1]; iray <= group_data.ubnd[1]; ++iray) {
            int const offset = THC_M0_INDEX(group_data, 0, iray);
            int const itheta = thc_sph_grid_get_itheta(M0Grid, iray);
            int const iphi   = thc_sph_grid_get_iphi(M0Grid, iray);

            /* 1st step: initialize excision mask */
            if(excision) {
                CCTK_REAL excision_rad = -1.0;
                if(sf_active[excision_surface]) {
                    int const sn = excision_surface;

                    CCTK_REAL const theta = thc_sph_grid_get_theta(M0Grid,iray);
                    CCTK_REAL const phi   = thc_sph_grid_get_phi(M0Grid, iray);

                    /* Copied from CarpetMask */
                    CCTK_REAL const theta0 = sf_origin_theta[sn];
                    CCTK_REAL const phi0   = sf_origin_phi  [sn];
                    CCTK_REAL const dtheta = sf_delta_theta[sn];
                    CCTK_REAL const dphi   = sf_delta_phi  [sn];

                    int const a = floor((theta - theta0) / dtheta + 0.5);
                    int const b = floor((phi   - phi0  ) / dphi   + 0.5);
                    CCTK_REAL const dr =
                        sf_radius[a + maxntheta * (b + maxnphi * sn)];
                    /* End of part copied from CarpetMask */

                    excision_rad = dr;
                }
                for(int irad = 0; irad < nrad; ++irad) {
                    int const ijk = THC_M0_INDEX(group_data, irad, iray);
                    CCTK_REAL const rad = thc_sph_grid_get_r(M0Grid, irad);
                    thc_M0_mask[ijk] = (rad < excision_rad) ||
                                       (thc_M0_alp[ijk] < excision_lapse);
                }
            }
            else {
                memset(&thc_M0_mask[offset], 0, nrad*sizeof(CCTK_INT));
            }

            /* 2nd step: compute radial null vector */
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);

                thc_M0_rad_null(irad, iray,
                        thc_M0_alp[ijk],
                        thc_M0_betax[ijk], thc_M0_betay[ijk], thc_M0_betaz[ijk],
                        thc_M0_gxx[ijk], thc_M0_gxy[ijk], thc_M0_gxz[ijk],
                        thc_M0_gyy[ijk], thc_M0_gyz[ijk], thc_M0_gzz[ijk],
                        thc_M0_zvecx[ijk], thc_M0_zvecy[ijk], thc_M0_zvecz[ijk],
                        &thc_M0_mask[ijk],
                        &kt[irad], &kr[irad], &chi[irad], &sqrt_det_g[irad]);

                thc_M0_flux_fac[ijk] = kr[irad]/kt[irad];
            }

            /* 3rd step: compute absorption opacities */
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);

                if(thc_M0_mask[ijk]) {
                    /* Also set the sources to zero in the excision region */
                    thc_M0_R_nue[ijk] = 0;
                    thc_M0_R_nua[ijk] = 0;
                    thc_M0_R_nux[ijk] = 0;
                    thc_M0_Q_nue[ijk] = 0;
                    thc_M0_Q_nua[ijk] = 0;
                    thc_M0_Q_nux[ijk] = 0;

                    thc_M0_abs_nue[ijk] = 0;
                    thc_M0_abs_nua[ijk] = 0;

                    thc_M0_N_nue_old[ijk] = 0;
                    thc_M0_N_nua_old[ijk] = 0;
                    thc_M0_N_nux_old[ijk] = 0;
                    thc_M0_E_nue_old[ijk] = 0;
                    thc_M0_E_nua_old[ijk] = 0;
                    thc_M0_E_nux_old[ijk] = 0;
                }
                else {
                    thc_M0_mask[ijk] = false;

                    /* Compute the free emissivities */
                    int ierr = NeutrinoEmission(
                            thc_M0_rho[ijk], thc_M0_temp[ijk], thc_M0_Y_e[ijk],
                            &thc_M0_R_nue[ijk], &thc_M0_R_nua[ijk],
                            &thc_M0_R_nux[ijk],
                            &thc_M0_Q_nue[ijk], &thc_M0_Q_nua[ijk],
                            &thc_M0_Q_nux[ijk]);
                    assert(!ierr);
                    assert(isfinite(thc_M0_R_nue[ijk]));
                    assert(isfinite(thc_M0_R_nua[ijk]));
                    assert(isfinite(thc_M0_R_nux[ijk]));
                    assert(isfinite(thc_M0_Q_nue[ijk]));
                    assert(isfinite(thc_M0_Q_nua[ijk]));
                    assert(isfinite(thc_M0_Q_nux[ijk]));

                    /* Get equilibrium neutrino density and energy density */
                    CCTK_REAL num_nue, num_nua, num_nux;
                    CCTK_REAL ene_nue, ene_nua, ene_nux;
                    ierr = NeutrinoDensity(
                            thc_M0_rho[ijk], thc_M0_temp[ijk], thc_M0_Y_e[ijk],
                            &num_nue, &num_nua, &num_nux,
                            &ene_nue, &ene_nua, &ene_nux);
                    assert(!ierr);
                    assert(isfinite(num_nue));
                    assert(isfinite(num_nua));
                    assert(isfinite(num_nux));
                    assert(isfinite(ene_nue));
                    assert(isfinite(ene_nua));
                    assert(isfinite(ene_nux));

                    /* Compute local opacities */
                    CCTK_REAL kappa_0_nue, kappa_0_nua, kappa_0_nux;
                    CCTK_REAL kappa_1_nue, kappa_1_nua, kappa_1_nux;
                    ierr = NeutrinoOpacity(
                            thc_M0_rho[ijk], thc_M0_temp[ijk], thc_M0_Y_e[ijk],
                            &kappa_0_nue, &kappa_0_nua, &kappa_0_nux,
                            &kappa_1_nue, &kappa_1_nua, &kappa_1_nux);
                    assert(!ierr);
                    assert(isfinite(kappa_0_nue));
                    assert(isfinite(kappa_0_nua));
                    assert(isfinite(kappa_0_nux));
                    assert(isfinite(kappa_1_nue));
                    assert(isfinite(kappa_1_nua));
                    assert(isfinite(kappa_1_nux));

                    /* Compute effective rates */
                    thc_M0_R_nue[ijk] = thc_lk_calc_eff_rate(thc_M0_R_nue[ijk],
                            num_nue, kappa_0_nue, thc_M0_optd_0_nue[ijk]);
                    thc_M0_R_nua[ijk] = thc_lk_calc_eff_rate(thc_M0_R_nua[ijk],
                            num_nua, kappa_0_nua, thc_M0_optd_0_nua[ijk]);
                    thc_M0_R_nux[ijk] = thc_lk_calc_eff_rate(thc_M0_R_nux[ijk],
                            num_nux, kappa_0_nux, thc_M0_optd_0_nux[ijk]);
                    thc_M0_Q_nue[ijk] = thc_lk_calc_eff_rate(thc_M0_Q_nue[ijk],
                            ene_nue, kappa_1_nue, thc_M0_optd_1_nue[ijk]);
                    thc_M0_Q_nua[ijk] = thc_lk_calc_eff_rate(thc_M0_Q_nua[ijk],
                            ene_nua, kappa_1_nua, thc_M0_optd_1_nua[ijk]);
                    thc_M0_Q_nux[ijk] = thc_lk_calc_eff_rate(thc_M0_Q_nux[ijk],
                            ene_nux, kappa_1_nux, thc_M0_optd_1_nux[ijk]);
                    assert(isfinite(thc_M0_R_nue[ijk]));
                    assert(isfinite(thc_M0_R_nua[ijk]));
                    assert(isfinite(thc_M0_R_nux[ijk]));
                    assert(isfinite(thc_M0_Q_nue[ijk]));
                    assert(isfinite(thc_M0_Q_nua[ijk]));
                    assert(isfinite(thc_M0_Q_nux[ijk]));

                    /* Compute effective absorption rates */
                    CCTK_REAL abs_0_nux, abs_1_nue, abs_1_nua, abs_1_nux;
                    ierr = NeutrinoAbsorptionRate(
                            thc_M0_rho[ijk], thc_M0_temp[ijk], thc_M0_Y_e[ijk],
                            &thc_M0_abs_nue[ijk], &thc_M0_abs_nua[ijk], &abs_0_nux,
                            &abs_1_nue, &abs_1_nua, &abs_1_nux);
                    assert(!ierr);
                    assert(isfinite(thc_M0_abs_nue[ijk]));
                    assert(isfinite(thc_M0_abs_nua[ijk]));
                    assert(isfinite(abs_0_nux));
                    assert(isfinite(abs_1_nue));
                    assert(isfinite(abs_1_nua));
                    assert(isfinite(abs_1_nux));

                    /* Scale neutrino opacity with the energy */
                    if(use_enedep_opacity) {
                        CCTK_REAL const abs_fac_nue = (ene_nue > 0 ?
                            (num_nue/ene_nue)*(thc_M0_E_nue[ijk]/chi[irad]) : 1);
                        assert(isfinite(abs_fac_nue));
                        thc_M0_abs_nue[ijk] *= SQ(abs_fac_nue);

                        CCTK_REAL const abs_fac_nua = (ene_nua > 0 ?
                            (num_nua/ene_nua)*(thc_M0_E_nua[ijk]/chi[irad]) : 1);
                        assert(isfinite(abs_fac_nua));
                        thc_M0_abs_nua[ijk] *= SQ(abs_fac_nua);
                    }
                    if(use_reduced_opacity) {
                        if(thc_M0_optd_0_nue[ijk] > FLT_EPSILON) {
                            thc_M0_abs_nue[ijk] *= exp(-thc_M0_optd_0_nue[ijk]);
                        }
                        if(thc_M0_optd_0_nua[ijk] > FLT_EPSILON) {
                            thc_M0_abs_nua[ijk] *= exp(-thc_M0_optd_0_nua[ijk]);
                        }
                    }
                }
            }

            /* 4th step: transport neutrinos  */
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                theta[irad] = kr[irad]/kt[irad];
                eta[irad] = thc_M0_R_nue[ijk]*sqrt_det_g[irad];
                mu[irad] = thc_M0_abs_nue[ijk]/kt[irad];
            }
            thc_M0_evol_density(dt, &thc_M0_mask[offset], theta, eta, mu,
                    &thc_M0_N_nue_old[offset], &thc_M0_N_nue[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_ndens_nue[ijk] = thc_M0_N_nue[ijk]/
                    (sqrt_det_g[irad]*kt[irad]);
            }

            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                eta[irad] = thc_M0_R_nua[ijk]*sqrt_det_g[irad];
                mu[irad] = thc_M0_abs_nua[ijk]/kt[irad];
            }
            thc_M0_evol_density(dt, &thc_M0_mask[offset], theta, eta, mu,
                    &thc_M0_N_nua_old[offset], &thc_M0_N_nua[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_ndens_nua[ijk] = thc_M0_N_nua[ijk]/
                    (sqrt_det_g[irad]*kt[irad]);
            }

            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                eta[irad] = thc_M0_R_nux[ijk]*sqrt_det_g[irad];
                mu[irad] = 0;
            }
            thc_M0_evol_density(dt, &thc_M0_mask[offset], theta, eta, mu,
                    &thc_M0_N_nux_old[offset], &thc_M0_N_nux[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_ndens_nux[ijk] = thc_M0_N_nux[ijk]/
                    (sqrt_det_g[irad]*kt[irad]);
            }

            /* 5th step: compute average energies */
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                eta[irad] = thc_M0_Q_nue[ijk]*(chi[irad]/kt[irad]);
                mu[irad] = thc_M0_R_nue[ijk]/kt[irad];
                assert(isfinite(eta[irad]));
                assert(isfinite(mu[irad]));
            }
            thc_M0_evol_energy_ave(dt, &thc_M0_mask[offset],
                    &thc_M0_ndens_nue[offset], theta, eta, mu,
                    &thc_M0_E_nue_old[offset], &thc_M0_E_nue[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_eave_nue[ijk] = thc_M0_E_nue[ijk]/chi[irad];
            }

            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                eta[irad] = thc_M0_Q_nua[ijk]*(chi[irad]/kt[irad]);
                mu[irad] = thc_M0_R_nua[ijk]/kt[irad];
                assert(isfinite(eta[irad]));
                assert(isfinite(mu[irad]));
            }
            thc_M0_evol_energy_ave(dt, &thc_M0_mask[offset],
                    &thc_M0_ndens_nua[offset], theta, eta, mu,
                    &thc_M0_E_nua_old[offset], &thc_M0_E_nua[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_eave_nua[ijk] = thc_M0_E_nua[ijk]/chi[irad];
            }

            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                eta[irad] = thc_M0_Q_nux[ijk]*(chi[irad]/kt[irad]);
                mu[irad] = thc_M0_R_nux[ijk]/kt[irad];
                assert(isfinite(eta[irad]));
                assert(isfinite(mu[irad]));
            }
            thc_M0_evol_energy_ave(dt, &thc_M0_mask[offset],
                    &thc_M0_ndens_nux[offset], theta, eta, mu,
                    &thc_M0_E_nux_old[offset], &thc_M0_E_nux[offset]);
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);
                thc_M0_eave_nux[ijk] = thc_M0_E_nux[ijk]/chi[irad];
            }

            /* 6th step: compute absorption / heating */
            for(int irad = 0; irad < nrad; ++irad) {
                int const ijk = THC_M0_INDEX(group_data, irad, iray);

                CCTK_REAL eff_abs_nue = thc_M0_abs_nue[ijk];
                CCTK_REAL eff_abs_nua = thc_M0_abs_nua[ijk];

                thc_M0_abs_number[ijk] = mb*(
                    eff_abs_nue * thc_M0_ndens_nue[ijk] -
                    eff_abs_nua * thc_M0_ndens_nua[ijk]);
                thc_M0_abs_energy[ijk] =
                    eff_abs_nue * thc_M0_ndens_nue[ijk] * thc_M0_eave_nue[ijk] +
                    eff_abs_nua * thc_M0_ndens_nua[ijk] * thc_M0_eave_nua[ijk];
            }

            /* 7th step: update fluxes at the outer boundary */
            int const ijk = THC_M0_INDEX(group_data, nrad-1, iray);
            /* Note: we have duplicated points on the axis and for phi = 2*pi */
            CCTK_REAL dS = dtheta*dphi;
            if(iphi == 0 || iphi == thc_sph_grid_get_nphi(M0Grid)  - 1) {
                dS = dS/2.0;
            }
            if(itheta == 0 || itheta == thc_sph_grid_get_ntheta(M0Grid) - 1) {
                dS = dS/thc_sph_grid_get_nphi(M0Grid);
            }

            CCTK_REAL wt = thc_M0_flux_fac[ijk] * dS;

            my_M0_nue_num_flux += thc_M0_N_nue[ijk] * wt;
            my_M0_nua_num_flux += thc_M0_N_nua[ijk] * wt;
            my_M0_nux_num_flux += thc_M0_N_nux[ijk] * wt;
            my_M0_nue_ene_flux += thc_M0_N_nue[ijk] * thc_M0_eave_nue[ijk] * wt;
            my_M0_nua_ene_flux += thc_M0_N_nua[ijk] * thc_M0_eave_nua[ijk] * wt;
            my_M0_nux_ene_flux += thc_M0_N_nux[ijk] * thc_M0_eave_nux[ijk] * wt;
        }
        free(mu);
        free(eta);
        free(theta);

        free(sqrt_det_g);
        free(chi);
        free(kr);
        free(kt);
    }

    /* Reduce fluxes */
    ierr = MPI_Allreduce(&my_M0_nue_num_flux, thc_M0_nue_num_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);
    ierr = MPI_Allreduce(&my_M0_nua_num_flux, thc_M0_nua_num_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);
    ierr = MPI_Allreduce(&my_M0_nux_num_flux, thc_M0_nux_num_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);
    ierr = MPI_Allreduce(&my_M0_nue_ene_flux, thc_M0_nue_ene_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);
    ierr = MPI_Allreduce(&my_M0_nua_ene_flux, thc_M0_nua_ene_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);
    ierr = MPI_Allreduce(&my_M0_nux_ene_flux, thc_M0_nux_ene_flux, 1,
            MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); assert(!ierr);

    *thc_leakage_M0_time = cctk_time;
}
