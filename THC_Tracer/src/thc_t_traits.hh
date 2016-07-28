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


#ifndef THC_T_TRAITS_HH
#define THC_T_TRAITS_HH

namespace thc {

class Tracer;

} // namespace thc

namespace hrscc {

template<>
class traits<thc::Tracer> {
    public:
        enum {nequations = 1};
        enum {nexternal = 8};
        enum {nbitmasks = 0};
        static bool const pure = true;
};

} // namespace hrscc

#endif
