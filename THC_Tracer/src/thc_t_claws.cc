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
#include <limits>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#include "thc_t_claws.hh"

static int conserved_t_idx[THC_T_MAX_NTRACERS];
static int primitive_t_idx[THC_T_MAX_NTRACERS];
static int rhs_t_idx[THC_T_MAX_NTRACERS];
static int num_flux_t_idx[3*THC_T_MAX_NTRACERS];

namespace hrscc {

template<>
int CLaw<thc::Tracer>::conserved_idx[1] = {0};

template<>
int CLaw<thc::Tracer>::primitive_idx[1] = {0};

template<>
int CLaw<thc::Tracer>::rhs_idx[1] = {0};

template<>
int CLaw<thc::Tracer>::field_idx[8] = {0, 0, 0, 0, 0, 0, 0, 0};

template<>
int CLaw<thc::Tracer>::bitmask_idx[0] = {};

template<>
int CLaw<thc::Tracer>::num_flux_idx[3] = {0, 0, 0};

template<>
CCTK_REAL CLaw<thc::Tracer>::conserved_lbound[1] = {
    std::numeric_limits<CCTK_REAL>::quiet_NaN()
};

} // namespace hrscc

using namespace hrscc;

extern "C"
void THC_T_HRSCCRegister(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(verbose) {
        CCTK_INFO("THC_T_HRSCCRegister");
    }

    char vname[BUFSIZ];

    // Barotropic and ultrarelativistic
    int n_hydro_equations = 4;
    if(CCTK_Equals(eos_type, "ideal")) {
        n_hydro_equations = 5;
    }
    else if(CCTK_Equals(eos_type, "nuclear")) {
        n_hydro_equations = 6;
    }

    // These indices are different for different tracers
    for(int i = 0; i < ntracers; ++i) {
        std::snprintf(vname, BUFSIZ, "THC_Tracer::tracer_dens[%d]", i);
        conserved_t_idx[i] = utils::cctk::var_index(vname);

        std::snprintf(vname, BUFSIZ, "THC_Tracer::tracer[%d]", i);
        primitive_t_idx[i] = utils::cctk::var_index(vname);

        std::snprintf(vname, BUFSIZ, "THC_Tracer::tracer_rhs[%d]", i);
        rhs_t_idx[i] = utils::cctk::var_index(vname);

        if(refluxing) {
            for(int d = 0; d < 3; ++d) {
                std::snprintf(vname, BUFSIZ, "Refluxing::flux[%d]",
                        3*n_hydro_equations + 3*i + d);
                num_flux_t_idx[3*i + d] = utils::cctk::var_index(vname);
            }
        }
    }

    // These indices are common for all tracers
    CLaw<thc::Tracer>::field_idx[0] = utils::cctk::var_index("ADMBase::alp");
    CLaw<thc::Tracer>::field_idx[1] = utils::cctk::var_index("ADMBase::betax");
    CLaw<thc::Tracer>::field_idx[2] = utils::cctk::var_index("ADMBase::betay");
    CLaw<thc::Tracer>::field_idx[3] = utils::cctk::var_index("ADMBase::betaz");
    CLaw<thc::Tracer>::field_idx[4] = utils::cctk::var_index("THC_Core::dens");
    CLaw<thc::Tracer>::field_idx[5] = utils::cctk::var_index("HydroBase::vel[0]");
    CLaw<thc::Tracer>::field_idx[6] = utils::cctk::var_index("HydroBase::vel[1]");
    CLaw<thc::Tracer>::field_idx[7] = utils::cctk::var_index("HydroBase::vel[2]");
}

namespace thc {

void set_tracer(int const tidx) {
    CLaw<thc::Tracer>::conserved_idx[0] = conserved_t_idx[tidx];
    CLaw<thc::Tracer>::primitive_idx[0] = primitive_t_idx[tidx];
    CLaw<thc::Tracer>::rhs_idx[0]       = rhs_t_idx[tidx];
    CLaw<thc::Tracer>::num_flux_idx[0]  = num_flux_t_idx[3*tidx + 0];
    CLaw<thc::Tracer>::num_flux_idx[1]  = num_flux_t_idx[3*tidx + 1];
    CLaw<thc::Tracer>::num_flux_idx[2]  = num_flux_t_idx[3*tidx + 2];
}

} // namespace thc
