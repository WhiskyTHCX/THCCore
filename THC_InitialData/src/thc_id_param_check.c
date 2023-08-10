//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2011, David Radice <david.radice@aei.mpg.de>
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


#include <stdio.h>
#include <string.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void THC_ID_ParamCheck(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if(CCTK_Equals(initial_Y_e, "THC_Initial")) {
        if(!CCTK_Equals(id_type, "atmosphere") &&
           !CCTK_Equals(id_type, "shocktube") &&
           !CCTK_Equals(id_type, "static")) {
            char buf[BUFSIZ];
            snprintf(buf, BUFSIZ, "id_type = \"%s\" does not support "
                    "setting Y_e yet!", id_type);
            CCTK_PARAMWARN(buf);
        }
    }
}
