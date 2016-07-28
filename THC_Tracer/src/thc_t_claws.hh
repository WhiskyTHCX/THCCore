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


#ifndef THC_T_CLAWS_HH
#define THC_T_CLAWS_HH

#include "hrscc.hh"

#include "thc_t_macro.hh"
#include "thc_t_traits.hh"

namespace thc {

void set_tracer(int const tidx);

class Tracer: public hrscc::CLaw<Tracer> {
    public:
        typedef hrscc::Observer<hrscc::CLaw<Tracer> > Observer;
        enum {nequations = hrscc::CLaw<Tracer>::nequations};
        enum {nexternal = hrscc::CLaw<Tracer>::nexternal};
        enum {nbitmasks = hrscc::CLaw<Tracer>::nbitmasks};

        inline void prim_to_all(
                Observer & observer
                ) const {
            TRR_MAKE_ALIAS(hrscc::policy::x, observer);
            tracer_dens = tracer * dens;
        }

        inline void cons_to_all(
                Observer & observer
                ) const {
            TRR_MAKE_ALIAS(hrscc::policy::x, observer);
            tracer = tracer_dens / dens;
        }

        template<hrscc::policy::direction_t dir>
        inline void fluxes(
                Observer & observer
                ) const {
            TRR_MAKE_ALIAS(dir, observer);
            flux_tracer = alp*tracer_dens*(velx - betax/alp);
        }

        template<hrscc::policy::direction_t dir>
        inline void eigenvalues(
                Observer & observer
                ) const {
            TRR_MAKE_ALIAS(dir, observer);
            eigenvalue = alp*(velx - betax/alp);
        }

        template<hrscc::policy::direction_t dir>
        inline void eig(
                Observer & observer
                ) const {
            TRR_MAKE_ALIAS(dir, observer);
            eigenvalue   = alp*(velx - betax/alp);
            left_eigvec  = 1.0;
            right_eigvec = 1.0;
        }
};

}

#endif
