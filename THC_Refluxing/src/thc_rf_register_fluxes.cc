//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2015, David Radice <dradice@caltech.edu>
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

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void THC_RF_RegisterFluxes(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    int ierr = 0;

    // Initialize indices
    CCTK_INT * indices_vec[] = {
        idx_dens, idx_densxn, idx_densxp, idx_sconx, idx_scony, idx_sconz,
        idx_tau, idx_tracer_0, idx_tracer_1, idx_tracer_2, idx_tracer_3, NULL
    };
    for(CCTK_INT ** p = &(indices_vec[0]); *p != NULL; ++p) {
        **p = -1;
    }

    // Register variables with the refluxing thorn
    // WARNING: the registration has to happen in the same order with which the
    // conservative variables are registered in HRSCCore
    if(CCTK_Equals(eos_type, "barotropic")) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::dens"),
                false, idx_dens);
    }
    if(CCTK_Equals(eos_type, "ideal")) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::dens"),
                false, idx_dens);
    }
    if(CCTK_Equals(eos_type, "nuclear")) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::densxn"),
                false, idx_densxn);
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::densxp"),
                false, idx_densxp);
    }
    if(CCTK_Equals(eos_type, "ultrarelativistic")) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::tau"),
                false, idx_tau);
    }

    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::scon[0]"),
            false, idx_sconx);
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::scon[1]"),
            false, idx_scony);
    ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::scon[2]"),
            false, idx_sconz);

    if(!CCTK_Equals(eos_type, "barotropic") &&
       !CCTK_Equals(eos_type, "ultrarelativistic")) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex("THC_Core::tau"),
                false, idx_tau);
    }

    if(ntracers > 0) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex(
                    "THC_Tracer::tracer_dens[0]"), false, idx_tracer_0);
    }
    if(ntracers > 1) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex(
                    "THC_Tracer::tracer_dens[1]"), false, idx_tracer_1);
    }
    if(ntracers > 2) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex(
                    "THC_Tracer::tracer_dens[2]"), false, idx_tracer_2);
    }
    if(ntracers > 3) {
        ierr |= RefluxingRegisterVariable(CCTK_VarIndex(
                    "THC_Tracer::tracer_dens[3]"), false, idx_tracer_3);
    }
    if(ntracers > 4) {
        CCTK_ERROR("Having more than 4 tracers is not currently supported");
    }

    assert(!ierr);
}
