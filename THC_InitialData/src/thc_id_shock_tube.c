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


#include <float.h>
#include <math.h>
#include <stdbool.h>
#define SQ(X) ((X)*(X))

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "utils_macro.h"

void THC_ID_ShockTube(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_REAL * velx = &vel[0*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * vely = &vel[1*UTILS_GFSIZE(cctkGH)];
    CCTK_REAL * velz = &vel[2*UTILS_GFSIZE(cctkGH)];

    if(verbose) {
        CCTK_INFO("THC_ID_ShockTube");
    }

    CCTK_REAL rhoL, velL, veltL, pressL, YeL;
    CCTK_REAL rhoR, velR, veltR, pressR, YeR;
    CCTK_REAL Gamma;

    bool const set_Y_e = CCTK_Equals(initial_Y_e, "THC_Initial");

    // This is for testing the advection of Ye
    YeL = 0.55;
    YeR = 0.45;

    CCTK_REAL const tracerL = 2.0;
    CCTK_REAL const tracerR = 1.0;

    if(CCTK_Equals(shocktube_case, "blast_wave")) {
        Gamma  = 5.0/3.0;

        rhoL   = 1.0e-3;
        velL   = 0.0;
        veltL  = 0.0;
        pressL = 1.0;

        rhoR   = 1.0e-3;
        velR   = 0.0;
        veltR  = 0.0;
        pressR = 1.0e-5;
    }
    else if(CCTK_Equals(shocktube_case, "collision")) {
        Gamma  = 4.0/3.0;

        rhoL   = 0.001;
        velL   = sqrt(1.0 - 1.0/SQ(collision_w_lorentz));
        veltL  = 0.0;
        pressL = 1000.0;

        rhoR   = 0.001;
        velR   = -velL;
        veltR  = 0.0;
        pressR = 1000.0;
    }
    else if(CCTK_Equals(shocktube_case, "contact")) {
        Gamma  = 5.0/3.0;

        rhoL   = 1.0;
        velL   = 0.5;
        veltL  = 0.0;
        pressL = 1.0;

        rhoR   = 0.5;
        velR   = 0.5;
        veltR  = 0.0;
        pressR = 1.0;
    }
    else if(CCTK_Equals(shocktube_case, "sod")) {
        Gamma  = 1.4;

        rhoL   = 1.0;
        velL   = 0.0;
        veltL  = 0.0;
        pressL = 1.0;

        rhoR   = 0.125;
        velR   = 0.0;
        veltR  = 0.0;
        pressR = 0.1;
    }
    else if(CCTK_Equals(shocktube_case, "strong_shock")) {
        Gamma  = 5.0/3.0;

        rhoL   = 10.0;
        velL   = 0.0;
        veltL  = 0.0;
        pressL = 10.0;

        rhoR   = 1.0;
        velR   = 0.0;
        veltR  = 0.0;
        pressR = 1.0e-5;
    }
    else if(CCTK_Equals(shocktube_case, "transverse")) {
        Gamma  = 5.0/3.0;

        rhoL   = 1.0e-3;
        velL   = 0.0;
        veltL  = 0.0;
        pressL = 1.0;

        rhoR   = 1.0e-3;
        velR   = 0.0;
        veltR  = 0.99;          // w_lorentz ~ 7
        pressR = 1.0e-5;
    }
    else if(CCTK_Equals(shocktube_case, "vacuum")) {
        Gamma   = 2.0;

        rhoL    = 1.0e-3;
        velL    = 0.0;
        veltL   = 0.0;
        pressL  = 1.0e-4;

        rhoR    = 0.0;
        velR    = 0.0;
        veltR   = 0.0;
        pressR  = 0.0;
    }
    else {
        abort();
    }

    CCTK_REAL const epsL = pressL/((Gamma - 1)*rhoL);
    CCTK_REAL const epsR = pressR/((Gamma - 1)*rhoR);

#pragma omp parallel
    {
        UTILS_LOOP3(thc_id_shock_tube,
                k, 0, cctk_lsh[2],
                j, 0, cctk_lsh[1],
                i, 0, cctk_lsh[0]) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            if(CCTK_Equals(shocktube_dir, "d")) {
                if(shocktube_width < DBL_EPSILON) {
                    if(x[ijk] + y[ijk] + z[ijk] < 0) {
                        rho[ijk]  = rhoL;
                        velx[ijk] = (velL - veltL) / sqrt(3.0);
                        vely[ijk] = (velL - veltL) / sqrt(3.0);
                        velz[ijk] = (velL - veltL) / sqrt(3.0);
                        eps[ijk]  = epsL;
                        if(set_Y_e) {
                            Y_e[ijk]  = YeL;
                        }
                        if(ntracers) {
                            tracer[ijk] = tracerL;
                        }
                    }
                    else {
                        rho[ijk]  = rhoR;
                        velx[ijk] = (velR - veltR) / sqrt(3.0);
                        vely[ijk] = (velR - veltR) / sqrt(3.0);
                        velz[ijk] = (velR - veltR) / sqrt(3.0);
                        eps[ijk]  = epsR;
                        if(set_Y_e) {
                            Y_e[ijk]  = YeR;
                        }
                        if(ntracers) {
                            tracer[ijk] = tracerR;
                        }
                    }
                }
                else {
                    CCTK_REAL const s = sqrt(SQ(x[ijk]) + SQ(y[ijk]) +
                            SQ(z[ijk])) / shocktube_width;
                    CCTK_REAL const v1 = (velL - veltL) / sqrt(3.0);
                    CCTK_REAL const v2 = (velR - veltR) / sqrt(3.0);

                    rho[ijk]  = 0.5*(rhoL + (2*rhoR - rhoL)*tanh(s));
                    velx[ijk] = 0.5*(v1 + (2*v2 - v1)*tanh(s));
                    vely[ijk] = 0.5*(v1 + (2*v2 - v1)*tanh(s));
                    velz[ijk] = 0.5*(v1 + (2*v2 - v1)*tanh(s));
                    eps[ijk]  = 0.5*(epsL + (2*epsR - epsL)*tanh(s));
                    if(set_Y_e) {
                        Y_e[ijk]  = 0.5*(YeL + (2*YeR - YeL)*tanh(s));
                    }

                    if(ntracers) {
                        tracer[ijk] = 0.5*(tracerL + (2*tracerR - tracerL)*
                                tanh(s));
                    }
                }
            }
            else {
                if(CCTK_Equals(shocktube_dir, "x")) {
                    if(shocktube_width < DBL_EPSILON) {
                        if(x[ijk] < 0) {
                            rho[ijk]  = rhoL;
                            velx[ijk] = velL;
                            vely[ijk] = veltL;
                            velz[ijk] = 0;
                            eps[ijk]  = epsL;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeL;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerL;
                            }
                        }
                        else {
                            rho[ijk]  = rhoR;
                            velx[ijk] = velR;
                            vely[ijk] = veltR;
                            velz[ijk] = 0;
                            eps[ijk]  = epsR;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeR;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerR;
                            }
                        }
                    }
                    else {
                        CCTK_REAL const s = tanh(x[ijk] / shocktube_width);

                        rho[ijk]  = 0.5*(rhoL + rhoR + (rhoR - rhoL)*s);
                        velx[ijk] = 0.5*(velL + velR + (velR - velL)*s);
                        vely[ijk] = 0.5*(veltL + veltR + (veltR - veltL)*s);
                        velz[ijk] = 0;
                        eps[ijk]  = 0.5*(epsL + epsR + (epsR - epsL)*s);
                        if(set_Y_e) {
                            Y_e[ijk]  = 0.5*(YeL + YeR + (YeR - YeL)*s);
                        }

                        if(ntracers) {
                            tracer[ijk] = 0.5*(tracerL + (2*tracerR - tracerL)*
                                    tanh(s));
                        }
                    }
                }
                else if(CCTK_Equals(shocktube_dir, "y")) {
                    if(shocktube_width < DBL_EPSILON) {
                        if(y[ijk] < 0) {
                            rho[ijk]  = rhoL;
                            velx[ijk] = 0;
                            vely[ijk] = velL;
                            velz[ijk] = veltL;
                            eps[ijk]  = epsL;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeL;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerL;
                            }
                        }
                        else {
                            rho[ijk]  = rhoR;
                            velx[ijk] = 0;
                            vely[ijk] = velR;
                            velz[ijk] = veltR;
                            eps[ijk]  = epsR;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeR;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerR;
                            }
                        }
                    }
                    else {
                        CCTK_REAL const s = y[ijk] / shocktube_width;

                        rho[ijk]  = 0.5*(rhoL + (2*rhoR - rhoL)*tanh(s));
                        velx[ijk] = 0;
                        vely[ijk] = 0.5*(velL + (2*velR - velL)*tanh(s));
                        velz[ijk] = 0.5*(veltL + (2*veltR - veltL)*tanh(s));
                        eps[ijk]  = 0.5*(epsL + (2*epsR - epsL)*tanh(s));
                        if(set_Y_e) {
                            Y_e[ijk]  = 0.5*(YeL + YeR + (YeR - YeL)*s);
                        }

                        if(ntracers) {
                            tracer[ijk] = 0.5*(tracerL + (2*tracerR - tracerL)*
                                    tanh(s));
                        }
                    }
                }
                else {
                    if(shocktube_width < DBL_EPSILON) {
                        if(z[ijk] < 0) {
                            rho[ijk]  = rhoL;
                            velx[ijk] = 0;
                            vely[ijk] = - veltL;
                            velz[ijk] = velL;
                            eps[ijk]  = epsL;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeL;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerL;
                            }
                        }
                        else {
                            rho[ijk]  = rhoR;
                            velx[ijk] = 0;
                            vely[ijk] = - veltR;
                            velz[ijk] = velR;
                            eps[ijk]  = epsR;
                            if(set_Y_e) {
                                Y_e[ijk]  = YeR;
                            }
                            if(ntracers) {
                                tracer[ijk] = tracerR;
                            }
                        }
                    }
                    else {
                        CCTK_REAL const s = z[ijk] / shocktube_width;

                        rho[ijk]  = 0.5*(rhoL + (2*rhoR - rhoL)*tanh(s));
                        velx[ijk] = 0;
                        vely[ijk] = - 0.5*(veltL + (2*veltR - veltL)*tanh(s));
                        velz[ijk] = 0.5*(velL + (2*velR - velL)*tanh(s));
                        eps[ijk]  = 0.5*(epsL + (2*epsR - epsL)*tanh(s));
                        if(set_Y_e) {
                            Y_e[ijk]  = 0.5*(YeL + YeR + (YeR - YeL)*s);
                        }

                        if(ntracers) {
                            tracer[ijk] = 0.5*(tracerL + (2*tracerR - tracerL)*
                                    tanh(s));
                        }
                    }
                }
            }
        } UTILS_ENDLOOP3(thc_id_shock_tube);
    }
}
