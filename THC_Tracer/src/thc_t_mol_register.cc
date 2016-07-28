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


#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void THC_T_MoLRegister(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_T_MoLRegister");
    }

    int ierr = 0;

    char prim_name[BUFSIZ];
    char cons_name[BUFSIZ];
    char rhs_name[BUFSIZ];

    for(int i = 0; i < ntracers; ++i) {
        std::snprintf(prim_name, BUFSIZ, "THC_Tracer::tracer[%d]", i);
        std::snprintf(cons_name, BUFSIZ, "THC_Tracer::tracer_dens[%d]", i);
        std::snprintf(rhs_name, BUFSIZ, "THC_Tracer::tracer_rhs[%d]", i);
        ierr |= MoLRegisterEvolved(CCTK_VarIndex(cons_name),
                CCTK_VarIndex(rhs_name));
        ierr |= MoLRegisterConstrained(CCTK_VarIndex(prim_name));
    }

    if(ierr) {
        CCTK_ERROR("Could not register with MoL");
    }
}
