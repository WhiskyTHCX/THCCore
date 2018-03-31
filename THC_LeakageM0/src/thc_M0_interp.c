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
#include <stdlib.h>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <util_Table.h>

#include "utils_macro.h"

#include "thc_M0_kernel.h"

#define length(X) (sizeof((X))/sizeof(*(X)))

#define MIN(X,Y) ((X)<(Y)?(X):(Y))

void THC_LK_NoAbsorption(CCTK_ARGUMENTS);

/* Interpolate the grid arrays onto the spherical grid */
void THC_M0_InterpToSph(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if((cctk_iteration-1) % compute_every != 0) {
        return;
    }
    if(!*thc_leakage_M0_is_on) {
        return;
    }

    if(verbose) {
        CCTK_INFO("THC_M0_InterpToSph");
    }

    int group_id = CCTK_GroupIndex("THC_LeakageM0::thc_leakage_vars");
    cGroupDynamicData group_data;
    int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
    assert(!ierr);
    assert(group_data.lsh[0] == nrad);
    assert(group_data.lbnd[0] == 0);

    int const interp_handle = CCTK_InterpHandle(interpolator);
    assert(interp_handle >= 0);
    int const options_handle =
        Util_TableCreateFromString(interpolator_options);
    assert(options_handle >= 0);
    int const coords_handle = CCTK_CoordSystemHandle("cart3d");
    assert(coords_handle >= 0);

    void const * interp_coords[] = {thc_M0_x, thc_M0_y, thc_M0_z};
    int const npoints = group_data.ash[0]*group_data.ash[1];
    CCTK_INT const input_array_indices[] = {
        CCTK_VarIndex("ADMBase::alp"),
        CCTK_VarIndex("ADMBase::betax"),
        CCTK_VarIndex("ADMBase::betay"),
        CCTK_VarIndex("ADMBase::betaz"),
        CCTK_VarIndex("ADMBase::gxx"),
        CCTK_VarIndex("ADMBase::gxy"),
        CCTK_VarIndex("ADMBase::gxz"),
        CCTK_VarIndex("ADMBase::gyy"),
        CCTK_VarIndex("ADMBase::gyz"),
        CCTK_VarIndex("ADMBase::gzz"),
        CCTK_VarIndex("HydroBase::rho"),
        CCTK_VarIndex("HydroBase::temperature"),
        CCTK_VarIndex("HydroBase::Y_e"),
        CCTK_VarIndex("THC_Core::zvec[0]"),
        CCTK_VarIndex("THC_Core::zvec[1]"),
        CCTK_VarIndex("THC_Core::zvec[2]"),
        CCTK_VarIndex("THC_LeakageBase::optd_0_nue"),
        CCTK_VarIndex("THC_LeakageBase::optd_0_nua"),
        CCTK_VarIndex("THC_LeakageBase::optd_0_nux"),
        CCTK_VarIndex("THC_LeakageBase::optd_1_nue"),
        CCTK_VarIndex("THC_LeakageBase::optd_1_nua"),
        CCTK_VarIndex("THC_LeakageBase::optd_1_nux")
    };
    int const ninputs = length(input_array_indices);

    CCTK_INT const output_array_types[] = {
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL
    };
    assert(ninputs == length(output_array_types));

    void * output_arrays[] = {
        thc_M0_alp,
        thc_M0_betax,
        thc_M0_betay,
        thc_M0_betaz,
        thc_M0_gxx,
        thc_M0_gxy,
        thc_M0_gxz,
        thc_M0_gyy,
        thc_M0_gyz,
        thc_M0_gzz,
        thc_M0_rho,
        thc_M0_temp,
        thc_M0_Y_e,
        thc_M0_zvecx,
        thc_M0_zvecy,
        thc_M0_zvecz,
        thc_M0_optd_0_nue,
        thc_M0_optd_0_nua,
        thc_M0_optd_0_nux,
        thc_M0_optd_1_nue,
        thc_M0_optd_1_nua,
        thc_M0_optd_1_nux
    };
    assert(ninputs == length(output_arrays));

    ierr = CCTK_InterpGridArrays(cctkGH, 3, interp_handle, options_handle,
            coords_handle, npoints, CCTK_VARIABLE_REAL, interp_coords, ninputs,
            input_array_indices, ninputs, output_array_types, output_arrays);
    assert(!ierr);

    Util_TableDestroy(options_handle);
}

/*
 * This routine calls CCTK_InterpLocalUniform to interpolate from the
 * spherical grid back onto the cartesian grid.
 *
 * 1. Each process creates a list of coordinates where to interpolate
 * 2. MPI is used to distribute the data on the spherical grid
 * 3. Each process interpolates everything back to the local cartesian grid
 *    at the coordinate locations computed at point 1.
 */
void THC_M0_InterpToCart(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(!*thc_leakage_M0_is_on) {
        THC_LK_NoAbsorption(CCTK_PASS_CTOC);
        return;
    }
    if((cctk_iteration-1) % compute_every != 0) {
        return;
    }

    if(verbose) {
        CCTK_INFO("THC_M0_InterpToCart");
    }

    /* We need CCTK_REAL to be doubles for the MPI calls */
    assert(sizeof(CCTK_REAL) == sizeof(double));

    int group_id = CCTK_GroupIndex("THC_LeakageM0::thc_leakage_vars");
    cGroupDynamicData group_data;
    int ierr = CCTK_GroupDynamicData(cctkGH, group_id, &group_data);
    assert(!ierr);
    assert(group_data.lsh[0] == nrad);
    assert(group_data.lbnd[0] == 0);

    /* Interpolation points */
    assert(cctk_lsh[0] == cctk_ash[0]);
    assert(cctk_lsh[1] == cctk_ash[1]);
    assert(cctk_lsh[2] == cctk_ash[2]);
    int const gfsiz     = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    CCTK_REAL * p_r     = malloc(gfsiz*sizeof(CCTK_REAL));
    CCTK_REAL * p_theta = malloc(gfsiz*sizeof(CCTK_REAL));
    CCTK_REAL * p_phi   = malloc(gfsiz*sizeof(CCTK_REAL));
#pragma omp parallel
    {
        UTILS_LOOP3(thc_M0_calc_interp_points,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
            thc_coord_cart_to_sph(x[ijk], y[ijk], z[ijk], &p_r[ijk],
                    &p_theta[ijk], &p_phi[ijk]);
            /* This effectively ensures all of the points out of the spherical
             * grid are filled with a zeroth order extrapolation */
            p_r[ijk] = MIN(p_r[ijk], rmax);
        } UTILS_ENDLOOP3(thc_M0_calc_interp_points);
    }

    /* Temporary arrays used to collect grid functions on the spherical grid */
    int const agsiz = nrad*nray;
    CCTK_REAL * glob_M0_abs_number = malloc(agsiz*sizeof(CCTK_REAL));
    CCTK_REAL * glob_M0_abs_energy = malloc(agsiz*sizeof(CCTK_REAL));

    /* Send all of the data around */
    int rank;
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); assert(!ierr);
    int nprocs;
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs); assert(!ierr);

    int istart = nrad*group_data.lbnd[1];
    int iend   = (nrad - 1) + nrad*(group_data.ubnd[1]);
    int isiz   = iend - istart + 1;

    int * recvcount = malloc(nprocs*sizeof(int));
    ierr = MPI_Allgather(&isiz, 1, MPI_INT, recvcount, 1, MPI_INT,
            MPI_COMM_WORLD);
    assert(!ierr);
    int * displs = malloc(nprocs*sizeof(int));
    ierr = MPI_Allgather(&istart, 1, MPI_INT, displs, 1, MPI_INT,
            MPI_COMM_WORLD);
    assert(!ierr);
    ierr = MPI_Allgatherv(thc_M0_abs_number, isiz, MPI_DOUBLE,
            glob_M0_abs_number, recvcount, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        assert(!ierr);
    ierr = MPI_Allgatherv(thc_M0_abs_energy, isiz, MPI_DOUBLE,
            glob_M0_abs_energy, recvcount, displs, MPI_DOUBLE, MPI_COMM_WORLD);
        assert(!ierr);

    /* Do the actual interpolation */
    int const interp_handle = CCTK_InterpHandle(interpolator);
    assert(interp_handle >= 0);
    int const options_handle =
        Util_TableCreateFromString(interpolator_options);
    assert(options_handle >= 0);

    CCTK_REAL coord_origin[3] = {0, 0, 0};
    CCTK_REAL coord_delta[3];
    thc_sph_grid_get_delta(M0Grid, &coord_delta[0], &coord_delta[1],
            &coord_delta[2]);
    void const * interp_coords[] = {p_r, p_theta, p_phi};

    void const * input_arrays[] = {
        glob_M0_abs_number,
        glob_M0_abs_energy
    };
    int const ninput = length(input_arrays);

    CCTK_INT const input_array_codes[] = {
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL
    };
    assert(length(input_array_codes) == ninput);

    CCTK_INT const input_array_dims[] = {nrad, thc_sph_grid_get_ntheta(M0Grid),
        thc_sph_grid_get_nphi(M0Grid)};

    CCTK_INT const output_array_codes[] = {
        CCTK_VARIABLE_REAL,
        CCTK_VARIABLE_REAL
    };
    assert(length(output_array_codes) == ninput);

    void * output_arrays[] = {
        abs_number,
        abs_energy
    };
    assert(length(output_arrays) == ninput);

    ierr = CCTK_InterpLocalUniform(3, interp_handle, options_handle,
            coord_origin, coord_delta, gfsiz, CCTK_VARIABLE_REAL,
            interp_coords, ninput, input_array_dims, input_array_codes,
            input_arrays, ninput, output_array_codes, output_arrays);
    assert(!ierr);

    /* Cleanup */
    free(recvcount);
    free(displs);

    free(glob_M0_abs_energy);
    free(glob_M0_abs_number);

    free(p_r);
    free(p_theta);
    free(p_phi);

    Util_TableDestroy(options_handle);
}
