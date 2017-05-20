
extern int riemann(const double *WL, const double *WR, double *W, double xovert, double tol, double vacuum); 
void anrs(const double *WL, const double *WR, double *ps, double *us); 

void exact_flux(const double *UL, const double *UR, double *F,char direction) {
    /* Flux using exact riemann solver.
     * U = (rho, rho*u, rho*v, rho*w, E, rho*si)
     * F = (rho*u, rho*u*u+p, rho*u*v, rho*u*w,u*(E+p), rho*u*si)
     * G = (rho*v, rho*u*v, rho*v*v+p, rho*v*w,u*(E+p), rho*v*si)
     * H = (rho*w, rho*u*w, rho*v*w, rho*w*w+p, rho*w*w,w*(E+p), rho*w*si)
     */
    double WL[3], WR[3], WS[3];
    /* Convert to primative variables */

    double g = GAMMA;
    WL[0] = UL[0];
    WL[2] = (g-1)*(UL[4] - .5*(UL[1]*UL[1] + UL[2]*UL[2] + UL[3]*UL[3])/UL[0]);
    WR[0] = UR[0];
    WR[2] = (g-1)*(UR[4] - .5*(UR[1]*UR[1] + UR[2]*UR[2] + UR[3]*UR[3])/UR[0]);
    if (direction == 'z') {
        WL[1] = UL[3]/UL[0];
        WR[1] = UR[3]/UR[0];
    }
    else if (direction == 'y') {
        WL[1] = UL[2]/UL[0];
        WR[1] = UR[2]/UR[0];
    }
    else {
        WL[1] = UL[1]/UL[0];
        WR[1] = UR[1]/UR[0];
    }

    int use_left = riemann(WL, WR,WS, 0.0, 1e-6, 1e-6); 


    F[0] = WS[0]*WS[1];
    F[4] = .5*F[0]*(WS[1]*WS[1]) + WS[1]*WS[2]; 
    F[5] = F[0]*UL[5]/UL[0] ? use_left : F[0]*UR[5]/UR[0];
    if (direction == 'z') {
        F[1] = F[0]*UL[1]/UL[0] ? use_left : F[0]*UR[1]/UR[0];
        F[2] = F[0]*UL[2]/UL[0] ? use_left : F[0]*UR[2]/UR[0];
        F[3] = F[0]*WS[1] + WS[2];
        F[4] += .5*F[0]*pow(UL[1]/UL[0],2) + .5*F[0]*pow(UL[2]/UL[0],2) ? use_left :  .5*F[0]*pow(UR[1]/UR[0],2) + .5*F[0]*pow(UR[2]/UR[0],2); 

    }
    else if (direction =='y') {
        F[1] = F[0]*UL[1]/UL[0] ? use_left : F[0]*UR[1]/UR[0];
        F[3] = F[0]*UL[3]/UL[0] ? use_left : F[0]*UR[3]/UR[0];
        F[2] = F[0]*WS[1] + WS[2];
        F[4] += .5*F[0]*pow(UL[1]/UL[0],2) + .5*F[0]*pow(UL[3]/UL[0],2) ? use_left :  .5*F[0]*pow(UR[1]/UR[0],2) + .5*F[0]*pow(UR[3]/UR[0],2); 
    }
    else {
        F[2] = F[0]*UL[2]/UL[0] ? use_left : F[0]*UR[2]/UR[0];
        F[3] = F[0]*UL[3]/UL[0] ? use_left : F[0]*UR[3]/UR[0];
        F[1] = F[0]*WS[1] + WS[2];
        F[4] += .5*F[0]*pow(UL[2]/UL[0],2) + .5*F[0]*pow(UL[3]/UL[0],2) ? use_left :  .5*F[0]*pow(UR[2]/UR[0],2) + .5*F[0]*pow(UR[3]/UR[0],2); 
    }

    return;
}

void hll_flux(const double *UL, const double *UR, double *F, char direction) {
    /* Flux using HLL scheme */
    double WL[3], WR[3], WS[3];
    double us, ps;
    /* Convert to primative variables */

    double g = GAMMA;
    WL[0] = UL[0];
    WL[2] = (g-1)*(UL[4] - .5*(UL[1]*UL[1] + UL[2]*UL[2] + UL[3]*UL[3])/UL[0]);
    WR[0] = UR[0];
    WR[2] = (g-1)*(UR[4] - .5*(UR[1]*UR[1] + UR[2]*UR[2] + UR[3]*UR[3])/UR[0]);
    if (direction == 'z') {
        WL[1] = UL[3]/UL[0];
        WR[1] = UR[3]/UR[0];
    }
    else if (direction == 'y') {
        WL[1] = UL[2]/UL[0];
        WR[1] = UR[2]/UR[0];
    }
    else {
        WL[1] = UL[1]/UL[0];
        WR[1] = UR[1]/UR[0];
    }


    anrs(WL,WR,&ps,&us);
    int use_left = (us > 0);
    
    double SL, SR;
    double rhoL = WL[0];
    double uL = WL[1];
    double pL = WL[2];
    double rhoR = WR[0];
    double uR = WR[1];
    double pR = WR[2];

    double aL = sqrt(g*pL/rhoL);
    double aR = sqrt(g*pR/rhoR);
    double gp = (g+1)/(2*g);

    if (ps <= pL) {
        SL = uL-aL;
    }
    else {
        SL = uL-aL*sqrt(1 + gp*(ps/pL-1));
    }
    if (ps <= pR) {
        SR = uR+aR;
    }
    else {
        SR = uR+aR*sqrt(1 + gp*(ps/pR-1));
    }

    if (SL >= 0) {
        /* Flux is due to left state */
        if (direction == 'z') {
            F[0] = UL[3];
            F[1] = UL[3]*UL[1]/UL[0];
            F[2] = UL[3]*UL[2]/UL[0];
            F[3] = UL[3]*WL[1] + WL[2];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
        else if (direction =='y') {
            F[0] = UL[2];
            F[1] = UL[2]*UL[1]/UL[0];
            F[3] = UL[2]*UL[3]/UL[0];
            F[2] = UL[2]*WL[1] + WL[2];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
        else {
            F[0] = UL[1];
            F[2] = UL[1]*UL[2]/UL[0];
            F[3] = UL[1]*UL[3]/UL[0];
            F[1] = UL[1]*WL[1] + WL[2];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
    }
    else {
        if (SR <= 0) {
            /* Flux due to right state */
            if (direction == 'z') {
                F[0] = UR[3];
                F[1] = UR[3]*UR[1]/UR[0];
                F[2] = UR[3]*UR[2]/UR[0];
                F[3] = UR[3]*WR[1] + WR[2];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
            else if (direction =='y') {
                F[0] = UR[2];
                F[1] = UR[2]*UR[1]/UR[0];
                F[3] = UR[2]*UR[3]/UR[0];
                F[2] = UR[2]*WR[1] + WR[2];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
            else {
                F[0] = UR[1];
                F[2] = UR[1]*UR[2]/UR[0];
                F[3] = UR[1]*UR[3]/UR[0];
                F[1] = UR[1]*WR[1] + WR[2];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
        }
        else {
            /* HLL Flux 
             * (SR*FL - SL*FR + SL*SR*(UR-UL) ) /(SR-SL)
             */
            if (direction == 'z') {
                F[0] = (SR*UL[3] - SL*UR[3] + SL*SR*(UR[0]-UL[0]))/(SR-SL);
                F[1] = ( SR*(UL[3]*UL[1]/UL[0]) - SL*(UR[3]*UR[1]/UR[0]) + SL*SR*(UR[1]-UL[1]))/(SR-SL);
                F[2] = ( SR*(UL[3]*UL[2]/UL[0]) - SL*(UR[3]*UR[2]/UR[0]) + SL*SR*(UR[2]-UL[2]))/(SR-SL);
                F[3] = ( SR*(UL[3]*WL[1] + WL[2]) - SL*(UR[3]*WR[1] + WR[2]) + SL*SR*(UR[3]-UL[3]))/(SR-SL);
                F[4] = ( SR*((UL[4] + WL[2])*WL[1]) - SL*((UR[4] + WR[2])*WR[1])+SL*SR*(UR[4]-UL[4]))/(SR-SL);
                F[5] = ( SR*(UL[5]*WL[1]) - SL*(UR[5]*WR[1]) +  SL*SR*(UR[5]-UL[5]))/(SR-SL);
            }
            else if (direction =='y') {
                F[0] = (SR*UL[2] - SL*UR[2] + SL*SR*(UR[0]-UL[0]))/(SR-SL);
                F[1] = ( SR*(UL[2]*UL[1]/UL[0]) - SL*(UR[2]*UR[1]/UR[0]) + SL*SR*(UR[1]-UL[1]))/(SR-SL);
                F[2] = ( SR*(UL[2]*WL[1] + WL[2]) - SL*(UR[2]*WR[1] + WR[2]) + SL*SR*(UR[2]-UL[2]))/(SR-SL);
                F[3] = ( SR*(UL[2]*UL[3]/UL[0]) - SL*(UR[2]*UR[3]/UR[0]) + SL*SR*(UR[3]-UL[3]))/(SR-SL);
                F[4] = ( SR*((UL[4] + WL[2])*WL[1]) - SL*((UR[4] + WR[2])*WR[1])+SL*SR*(UR[4]-UL[4]))/(SR-SL);
                F[5] = ( SR*(UL[5]*WL[1]) - SL*(UR[5]*WR[1]) +  SL*SR*(UR[5]-UL[5]))/(SR-SL);

            }
            else {
                F[0] = (SR*UL[1] - SL*UR[1] + SL*SR*(UR[0]-UL[0]))/(SR-SL);
                F[1] = ( SR*(UL[1]*WL[1] + WL[2]) - SL*(UR[1]*WR[1] + WR[2]) + SL*SR*(UR[1]-UL[1]))/(SR-SL);
                F[2] = ( SR*(UL[1]*UL[2]/UL[0]) - SL*(UR[1]*UR[2]/UR[0]) + SL*SR*(UR[2]-UL[2]))/(SR-SL);
                F[3] = ( SR*(UL[1]*UL[3]/UL[0]) - SL*(UR[1]*UR[3]/UR[0]) + SL*SR*(UR[3]-UL[3]))/(SR-SL);
                F[4] = ( SR*((UL[4] + WL[2])*WL[1]) - SL*((UR[4] + WR[2])*WR[1])+SL*SR*(UR[4]-UL[4]))/(SR-SL);
                F[5] = ( SR*(UL[5]*WL[1]) - SL*(UR[5]*WR[1]) +  SL*SR*(UR[5]-UL[5]))/(SR-SL);
            }
        }
    }
    
    return;
}

void hllc_flux(const double *UL, const double *UR, double *F, char direction) {
    /* Flux using HLLC scheme */
    double WL[3], WR[3], WS[3];
    double us, ps;
    /* Convert to primative variables */

    double g = GAMMA;
    WL[0] = UL[0];
    WL[2] = (g-1)*(UL[4] - .5*(UL[1]*UL[1] + UL[2]*UL[2] + UL[3]*UL[3])/UL[0]);
    WR[0] = UR[0];
    WR[2] = (g-1)*(UR[4] - .5*(UR[1]*UR[1] + UR[2]*UR[2] + UR[3]*UR[3])/UR[0]);
    if (direction == 'z') {
        WL[1] = UL[3]/UL[0];
        WR[1] = UR[3]/UR[0];
    }
    else if (direction == 'y') {
        WL[1] = UL[2]/UL[0];
        WR[1] = UR[2]/UR[0];
    }
    else {
        WL[1] = UL[1]/UL[0];
        WR[1] = UR[1]/UR[0];
    }


    anrs(WL,WR,&ps,&us);
    int use_left = (us > 0);
    
    double SL, SR, Sstar;
    double rhoL = WL[0];
    double uL = WL[1];
    double pL = WL[2];
    double rhoR = WR[0];
    double uR = WR[1];
    double pR = WR[2];

    double aL = sqrt(g*pL/rhoL);
    double aR = sqrt(g*pR/rhoR);
    double gp = (g+1)/(2*g);
    double ufac,efac;


    if (ps <= pL) {
        SL = uL-aL;
    }
    else {
        SL = uL-aL*sqrt(1 + gp*(ps/pL-1));
    }
    if (ps <= pR) {
        SR = uR+aR;
    }
    else {
        SR = uR+aR*sqrt(1 + gp*(ps/pR-1));
    }
    Sstar = ((pR-pL) + rhoL*uL*(SL-uL) - rhoR*uR*(SR-uR))/(rhoL*(SL-uL) - rhoR*(SR-uR));

    if (SL >= 0) {
        /* Flux is due to left state */
        if (direction == 'z') {
            F[0] = UL[3];
            F[1] = UL[3]*UL[1]/UL[0];
            F[2] = UL[3]*UL[2]/UL[0];
            F[3] = UL[3]*WL[1] + WL[2];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
        else if (direction =='y') {
            F[0] = UL[2];
            F[1] = UL[2]*UL[1]/UL[0];
            F[2] = UL[2]*WL[1] + WL[2];
            F[3] = UL[2]*UL[3]/UL[0];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
        else {
            F[0] = UL[1];
            F[1] = UL[1]*WL[1] + WL[2];
            F[2] = UL[1]*UL[2]/UL[0];
            F[3] = UL[1]*UL[3]/UL[0];
            F[4] = (UL[4] + WL[2])*WL[1];
            F[5] = UL[5]*WL[1];
        }
    }
    else {
        if (SR <= 0) {
            /* Flux due to right state */
            if (direction == 'z') {
                F[0] = UR[3];
                F[1] = UR[3]*UR[1]/UR[0];
                F[2] = UR[3]*UR[2]/UR[0];
                F[3] = UR[3]*WR[1] + WR[2];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
            else if (direction =='y') {
                F[0] = UR[2];
                F[2] = UR[2]*WR[1] + WR[2];
                F[1] = UR[2]*UR[1]/UR[0];
                F[3] = UR[2]*UR[3]/UR[0];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
            else {
                F[0] = UR[1];
                F[1] = UR[1]*WR[1] + WR[2];
                F[2] = UR[1]*UR[2]/UR[0];
                F[3] = UR[1]*UR[3]/UR[0];
                F[4] = (UR[4] + WR[2])*WR[1];
                F[5] = UR[5]*WR[1];
            }
        }
        else {
            if (Sstar >= 0) {
            /* HLLC Flux 
             * FL +SL*(UstarL-UL)
             */
                ufac= WL[0]*(SL - WL[1])/(SL-Sstar);
                efac = (Sstar-WL[1])*(Sstar + pL/(rhoL*(SL-WL[1])));
                if (direction == 'z') {
                    F[0] = UL[3];
                    F[0] += SL * ( ufac - UL[0]);
                    F[1] = UL[3]*UL[1]/UL[0];
                    F[1] += SL * ( ufac*UL[1]/UL[0] - UL[1]);
                    F[2] = UL[3]*UL[2]/UL[0];
                    F[2] += SL * ( ufac*UL[2]/UL[0] - UL[2]);
                    F[3] = UL[3]*WL[1] + WL[2];
                    F[3] += SL * ( ufac*Sstar - UL[3]);
                    F[4] = (UL[4] + WL[2])*WL[1];
                    F[4] += SL * ( ufac*(UL[4]/UL[0] + efac));
                    F[5] = UL[5]*WL[1];
                    F[5] += SL * ( ufac*UL[5]/UL[0] - UL[5]);
                }
                else if (direction =='y') {
                    F[0] = UL[2];
                    F[0] += SL * ( ufac - UL[0]);
                    F[1] = UL[2]*UL[1]/UL[0];
                    F[1] += SL * ( ufac*UL[1]/UL[0] - UL[1]);
                    F[2] = UL[2]*WL[1] + WL[2];
                    F[2] += SL * ( ufac*Sstar - UL[2]);
                    F[3] = UL[2]*UL[3]/UL[0];
                    F[3] += SL * ( ufac*UL[3]/UL[0] - UL[3]);
                    F[4] = (UL[4] + WL[2])*WL[1];
                    F[4] += SL * ( ufac*(UL[4]/UL[0] + efac));
                    F[5] = UL[5]*WL[1];
                    F[5] += SL * ( ufac*UL[5]/UL[0] - UL[5]);
                }
                else {
                    F[0] = UL[1];
                    F[0] += SL * ( ufac - UL[0]);
                    F[1] = UL[1]*WL[1] + WL[2];
                    F[1] += SL * ( ufac*Sstar - UL[1]);
                    F[2] = UL[1]*UL[2]/UL[0];
                    F[2] += SL * ( ufac*UL[2]/UL[0] - UL[2]);
                    F[3] = UL[1]*UL[3]/UL[0];
                    F[3] += SL * ( ufac*UL[3]/UL[0] - UL[3]);
                    F[4] = (UL[4] + WL[2])*WL[1];
                    F[4] += SL * ( ufac*(UL[4]/UL[0] + efac));
                    F[5] = UL[5]*WL[1];
                    F[5] += SL * ( ufac*UL[5]/UL[0] - UL[5]);
                }

            }
            else {
            /* HLLC Flux 
             * FR +SR*(UstarR-UR)
             */
                ufac= WR[0]*(SR - WR[1])/(SR-Sstar);
                efac = (Sstar-WR[1])*(Sstar + pR/(rhoR*(SR-WR[1])));
                if (direction == 'z') {
                    F[0] = UR[3];
                    F[0] += SR * ( ufac - UR[0]);
                    F[1] = UR[3]*UR[1]/UR[0];
                    F[1] += SR * ( ufac*UR[1]/UR[0] - UR[1]);
                    F[2] = UR[3]*UR[2]/UR[0];
                    F[2] += SR * ( ufac*UR[2]/UR[0] - UR[2]);
                    F[3] = UR[3]*WR[1] + WR[2];
                    F[3] += SR * ( ufac*Sstar - UR[3]);
                    F[4] = (UR[4] + WR[2])*WR[1];
                    F[4] += SR * ( ufac*(UR[4]/UR[0] + efac));
                    F[5] = UR[5]*WR[1];
                    F[5] += SR * ( ufac*UR[5]/UR[0] - UR[5]);
                }
                else if (direction =='y') {
                    F[0] = UR[2];
                    F[0] += SR * ( ufac - UR[0]);
                    F[1] = UR[2]*UR[1]/UR[0];
                    F[1] += SR * ( ufac*UR[1]/UR[0] - UR[1]);
                    F[2] = UR[2]*WR[1] + WR[2];
                    F[2] += SR * ( ufac*Sstar - UR[2]);
                    F[3] = UR[2]*UR[3]/UR[0];
                    F[3] += SR * ( ufac*UR[3]/UR[0] - UR[3]);
                    F[4] = (UR[4] + WR[2])*WR[1];
                    F[4] += SR * ( ufac*(UR[4]/UR[0] + efac));
                    F[5] = UR[5]*WR[1];
                    F[5] += SR * ( ufac*UR[5]/UR[0] - UR[5]);
                }
                else {
                    F[0] = UR[1];
                    F[0] += SR * ( ufac - UR[0]);
                    F[1] = UR[1]*WR[1] + WR[2];
                    F[1] += SR * ( ufac*Sstar - UR[1]);
                    F[2] = UR[1]*UR[2]/UR[0];
                    F[2] += SR * ( ufac*UR[2]/UR[0] - UR[2]);
                    F[3] = UR[1]*UR[3]/UR[0];
                    F[3] += SR * ( ufac*UR[3]/UR[0] - UR[3]);
                    F[4] = (UR[4] + WR[2])*WR[1];
                    F[4] += SR * ( ufac*(UR[4]/UR[0] + efac));
                    F[5] = UR[5]*WR[1];
                    F[5] += SR * ( ufac*UR[5]/UR[0] - UR[5]);
                }

            }
        }
    }
    return;
}
