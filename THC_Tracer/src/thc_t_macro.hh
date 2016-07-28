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


#ifndef THC_T_MACRO_HH
#define THC_T_MACRO_HH

#include "thc_macro.hh"

#define THC_T_MAX_NTRACERS 4

#define TRR_MAKE_ALIAS(dir,observer)                                        \
    cGH const * cctkGH       = observer.cctkGH;                             \
    CCTK_REAL & tracer_dens  = observer.conserved[0];                       \
    CCTK_REAL & tracer       = observer.primitive[0];                       \
    CCTK_REAL & flux_tracer  = observer.flux[dir][0];                       \
    CCTK_REAL & alp          = observer.field[0];                           \
    CCTK_REAL & betax        = observer.field[1+thc::index<dir>::x];        \
    CCTK_REAL & betay        = observer.field[1+thc::index<dir>::y];        \
    CCTK_REAL & betaz        = observer.field[1+thc::index<dir>::z];        \
    CCTK_REAL & dens         = observer.field[4];                           \
    CCTK_REAL & velx         = observer.field[5 + thc::index<dir>::x];      \
    CCTK_REAL & vely         = observer.field[5 + thc::index<dir>::y];      \
    CCTK_REAL & velz         = observer.field[5 + thc::index<dir>::z];      \
    CCTK_REAL & eigenvalue   = observer.eigenvalue[0];                      \
    CCTK_REAL & left_eigvec  = observer.left_eigenvector[0][0];             \
    CCTK_REAL & right_eigvec = observer.right_eigenvector[0][0];            \
                                                                            \
    UNUSED(cctkGH);                                                         \
    UNUSED(tracer_dens);                                                    \
    UNUSED(tracer);                                                         \
    UNUSED(flux_tracer);                                                    \
    UNUSED(alp);                                                            \
    UNUSED(betax);                                                          \
    UNUSED(betay);                                                          \
    UNUSED(betaz);                                                          \
    UNUSED(dens);                                                           \
    UNUSED(velx);                                                           \
    UNUSED(vely);                                                           \
    UNUSED(velz);                                                           \
    UNUSED(eigenvalue);                                                     \
    UNUSED(left_eigvec);                                                    \
    UNUSED(right_eigvec)


#endif
