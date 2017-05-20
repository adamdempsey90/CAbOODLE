#define QUSER 2.0


/* Approximate solver functions */
void pvrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Primative variable approximate riemman solver.
   * Used in smooth regions to give estimates for pressure and velocity in star region.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  double g = GAMMA;
  double a_L = sqrt(g*p_L/rho_L);
  double a_R = sqrt(g*p_R/rho_R);
  
  double rho_a_bar = .25*(rho_L + rho_R)*(a_L + a_R);
  
  *ps = .5*(p_L + p_R) - .5*(u_R - u_L) * rho_a_bar;
  *us = .5*(u_L + u_R) - .5*(p_R-p_L)/rho_a_bar;
  
  return;

  
}



void trrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Two rarefaction approximate riemman solver.
   * Used near rarefaction waves to give estimates for pressure and velocity in star region.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  double g = GAMMA;
  double a_L = sqrt(g*p_L/rho_L);
  double a_R = sqrt(g*p_R/rho_R);
   
  double gamma_1 = (g - 1)/2.;
  double inv_gamma_1 = 1./gamma_1;
  double gamma_2 = (g - 1)/(2.*g);
  double inv_gamma_2 =1/gamma_2;
  double plr = pow(p_L/p_R,gamma_2);
  
  
  *ps = pow( ( a_L + a_R - gamma_1*(u_R - u_L))/(a_L*pow(p_L,-gamma_2) + a_R*pow(p_R,-gamma_2) ), inv_gamma_2);
  *us = (plr*u_L/a_L + u_R/a_R + inv_gamma_1 * (plr - 1) )/(plr/a_L + 1/a_R);
  
  
  return;

  
}



void tsrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Two shock approximate riemman solver.
   * Used near shocks waves to give estimates for pressure and velocity in star region.
   * ps contains the pressure estimate from the pvrs function and is overwritten.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  double g = GAMMA;
  
  
  double gamma_1 = 2./(g  + 1);
  double gamma_2 = (g -1)/(g + 1);
  
  double p0 = fmax(0.0, *ps);
  
  double Al = gamma_1/rho_L;
  double Ar = gamma_1/rho_R;
  double Bl = gamma_2 * p_L;
  double Br = gamma_2 * p_R;
  double gl = sqrt( Al/(p0 + Bl));
  double gr = sqrt(Ar/(p0 +Br));
  
  *ps = (gl * p_L + gr * p_R - (u_R-u_L))/(gl + gr);
  *us = .5*(u_L + u_R) + .5*( (*ps - p_R) *gr - (*ps - p_L)*gl);
  
  
  return;

  
}


void anrs(const double *WL, const double *WR, double *ps, double *us) {
  
  
  double pmax = fmax(WL[2],WR[2]);
  double pmin = fmin(WL[2],WR[2]);
  double Q = pmax/pmin;
  double pstar, ustar;
 
  
  pvrs(WL,WR,&pstar,&ustar);
  
  
  if ( (Q < QUSER) && (pstar < pmax) && (pstar > pmin)) {
    
    *ps = pstar;
    *us = ustar;
    
  }
  else {
    
    if ( pstar < pmin ) {
     
      trrs(WL, WR, ps, us);
      
    }
    else {
      tsrs(WL, WR, &pstar, us);
      *ps = pstar;
    }
  }
  
  return;
  
}

/* Exact solver functions */
double fk_prime(const double p, const double pk, const double rhok, const double Ak, const double Bk, const double ak, const double z, const double g) {

    if (p > pk) {
        return ( 1 - (p-pk)/(2*(Bk+p)))*sqrt(Ak/(p+Bk));
    }
    else {
        return 1./(ak*rhok) * pow(p/pk,-z-1./g);
    }

}

double fk(const double p, const double pk, const double Ak, const double Bk, const double ak, const double z, const double g) {

    if (p > pk) {
        return (p-pk)*sqrt(Ak/(p+Bk));
    }
    else {
        return 2*ak/(g-1) *( pow(p/pk,z) - 1);
    }

}

double exact_func(const double p, const double du, const double AL, const double AR, const double BL, const double BR, const double pL, const double pR, const double aL, const double aR, const double z,const double g) {
    /* Evaluate f(p) = fL + fR + du */
    
    return fk(p,pL,AL,BL,aL,z,g) + fk(p,pR,AR,BR,aR,z,g) + du;

}
double exact_func_prime(const double p, const double AL, const double AR, const double BL, const double BR, const double pL, const double pR, const double rhoL, const double rhoR, const double aL, const double aR, const double z,const double g) {
    /* Evaluate f'(p) = fL' + fR'  */
    
    return fk_prime(p,pL,rhoL,AL,BL,aL,z,g) + fk_prime(p,pR,rhoR,AR,BR,aR,z,g);

}

void exact_riemann(const double *WL, const double *WR, double *ps, double *us, double tol) {
    /* Exact Riemann solver */


    double pstar, ustar;

    /* Use ANRS for initial guess */
    anrs(WL,WR,&pstar,&ustar);

    double rho_L = WL[0];
    double u_L = WL[1];
    double p_L = WL[2];
    double rho_R = WR[0];
    double u_R = WR[1];
    double p_R = WR[2];
    double g = GAMMA;
    double du = u_R - u_L;
    
    double a_L = sqrt(g*p_L/rho_L);
    double a_R = sqrt(g*p_R/rho_R);


    double z = (g-1)/(2*g);
    double AL = 2/((g+1)*rho_L);
    double AR = 2/((g+1)*rho_R);
    double BL = (g-1)*p_L/(g+1);
    double BR = (g-1)*p_R/(g+1);
    double dt = u_R - u_L;
    
    /* Newton-Raphson iteration to convergence */
    double resid = 1.0;
    double pold = pstar;
    while (resid > tol) {

        pstar = pold - exact_func(pold,du,AL,AR,BL,BR,p_L,p_R,a_L,a_R,z,g)/exact_func_prime(pold,AL,AR,BL,BR,p_L,p_R,rho_L,rho_R,a_L,a_R,z,g); 

        resid = fabs(pstar-pold)*2/(pold+pstar);
    }


    ustar = .5*( u_L + u_R  + fk(pstar,p_R,AR,BR,a_R,z,g) - fk(pstar,p_L,AL,BL,a_L,z,g) );

    *us = ustar;
    *ps = pstar;

    return; 


}  


int riemann(const double *WL, const double *WR, double *W, double xovert, double tol, double vacuum) {
    /* Full solution to the Riemann problem with exact solver at x/t 
     * includes support for vacuum.
     */

    double ps, us;

    /* Sample solution for left and right density in star region */
    double rhoL = WL[0];
    double uL = WL[1];
    double pL = WL[2];
    double rhoR = WR[0];
    double uR = WR[1];
    double pR = WR[2];
    
    double g = GAMMA;
    double gp = (g+1)/(2*g);
    double gm = (g-1)/(2*g);

    double aL ,aR;

    double rho_f, u_f, p_f;

    double S, St;


    /* Test for vacuum states */
    if ( (rhoL <= vacuum) || (pL <= vacuum)) {
        /* Left vacuum state */


        aR = sqrt(g*pR/rhoR);
        S = uR - 2*aR/(g-1);

        if (uR+aR <= xovert) {
            /* Right of rarefaction */
            u_f = uR;
            p_f = pR;
            rho_f = rhoR;
        }
        else {
            if (S >= xovert) {
                /* Left of rarefaction */
                u_f = S;
                p_f = vacuum;
                rho_f = vacuum;
            }
            else {
                /* In rarefaction fan */

                p_f = 2./(g+1) - (gm/gp)*(uR - xovert)/aR;
                rho_f = rhoR*pow(p_f, 2./(g-1));
                p_f = pR*pow(p_f,1./gm);
                u_f = (2./(g+1))*( -aR + .5*(g-1)*uR + xovert);
            }
        }
        W[0] = rho_f;
        W[1] = u_f;
        W[2]= p_f;

        return (S >= xovert);
    }
    if ( (rhoR <= vacuum) || (pR <= vacuum)) {
        /* Right vacuum state */
        aL = sqrt(g*pL/rhoL);
        S = uL + 2*aL/(g-1);

        if ( uL-aL >= xovert) {
            /* Left of rarefaction */
            u_f = uL;
            p_f = pL;
            rho_f = rhoL;
        }
        else {
            if (S <= xovert) {
                /* Right of rarefaction in vacuum */
                u_f = S;
                p_f = vacuum;
                rho_f = vacuum;
            }
            else {
                /* In rarefaction fan */
                rho_f = 2./(g+1) + (gm/gp)*(uL - xovert)/aL;
                p_f = pL*pow(rho_f, 1./gm);
                rho_f = rhoL*pow(rho_f, 2./(g-1));
                u_f = (aL + .5*(g-1)*uL + xovert) * 2./(g+1);
            }
        }
        W[0] = rho_f;
        W[1] = u_f;
        W[2] = p_f;
        return (S >= xovert);
    }

    /* Left and right stats are not vacuum */
    aL = sqrt(g*pL/rhoL);
    aR = sqrt(g*pR/rhoR);


    /* Test whether a vacuum will be created */

    S = uL + 2*aL/(g-1);
    St = uR - 2*aR/(g-1);

    if ( S <= St ) {
        /* Vacuum condition satisfied */

        if (S >= xovert) {
            /* In left rarefaction */
            if ( uL-aL >= xovert) {
                /* Left of rarefaction */
                u_f = uL;
                p_f = pL;
                rho_f = rhoL;
            }
            else {
                if (S <= xovert) {
                    /* Right of rarefaction in vacuum */
                    u_f = S;
                    p_f = vacuum;
                    rho_f = vacuum;
                }
                else {
                    /* In rarefaction fan */
                    rho_f = 2./(g+1) + (gm/gp)*(uL - xovert)/aL;
                    p_f = pL*pow(rho_f, 1./gm);
                    rho_f = rhoL*pow(rho_f, 2./(g-1));
                    u_f = (aL + .5*(g-1)*uL + xovert) * 2./(g+1);
                }
            }
        }
        else {
            if (St <= xovert) {
                /* In right rarefaction */
                if (uR+aR <= xovert) {
                    /* Right of rarefaction */
                    u_f = uR;
                    p_f = pR;
                    rho_f = rhoR;
                }
                else {
                    if (St >= xovert) {
                        /* Left of rarefaction */
                        u_f = St;
                        p_f = vacuum;
                        rho_f = vacuum;
                    }
                    else {
                        /* In rarefaction fan */
                        p_f = 2./(g+1) - (gm/gp)*(uR - xovert)/aR;
                        rho_f = rhoR*pow(p_f, 2./(g-1));
                        p_f = pR*pow(p_f,1./gm);
                        u_f = (2./(g+1))*( -aR + .5*(g-1)*uR + xovert);
                    }
                }
            }
            else {
                /* In vacuum */
                u_f = S;
                p_f = vacuum;
                rho_f = vacuum;
            }
        }

        return (S >= xovert); 

    }
    /* Pressure positivity condition satisfied */

    /* Get the pressure and velocity in the star region */
    exact_riemann(WL,WR,&ps,&us,tol);

    if (us >= xovert) {
        /* Left side of contact */
        if (ps > pL) {
            /* Left shock wave */

            S = uL - aL*sqrt(gp*(ps/pL) + gm);  
            if (S < xovert) {
                /* Right of shock */
                rho_f = rhoL*( ps/pL + gm/gp)/( (gm/gp)*(ps/pL) + 1);
                u_f = us;
                p_f = ps;
            }
            else {
                /* Left of shock */
                rho_f = rhoL;
                u_f = uL;
                p_f = pL;
            }

        }
        else {
            /* Left rarefaction wave */
            S = uL-aL;
            St = us - aL*pow( ps/pL, gm);
            if (S >= xovert) {
                /* Left of rarefaction head */
                rho_f = rhoL;
                u_f = uL;
                p_f = pL;
            }
            else {
                if (St <= xovert) {
                    /* Right of rarefaction tail */
                    rho_f = rhoL*pow(ps/pL,1./g);
                    u_f = us;
                    p_f = ps;
                }
                else {
                    /* Inside rarefaction fan */

                    rho_f = 2./(g+1) + (gm/gp)*(uL - xovert)/aL;
                    p_f = pL*pow(rho_f, 1./gm);
                    rho_f = rhoL*pow(rho_f, 2./(g-1));
                    u_f = (aL + .5*(g-1)*uL + xovert) * 2./(g+1);

                }
            }
        }
    }
    else {
        /* Right side of contact */

        if (ps > pR) {
            /* Right shock wave */
            S = uR + aR*sqrt( gp*(ps/pR) + gm);

            if (S <= xovert) {
                /* Right of shock */
                p_f =  pR;
                rho_f = rhoR;
                u_f = uR;
            }
            else {
                /* Left of shock */
                u_f = us;
                p_f = ps;
                rho_f= rhoR*( ps/pR + gm/gp)/( (gm/gp)*(ps/pR) + 1);
            }
        }
        else {
            /* Right rarefaction wave */
            S =  uR + aR;
            St = us + aR*pow(ps/pR,gm);
            if ( S <= xovert) {
                /* Right of rarefaction head */
                u_f = uR;
                p_f = pR;
                rho_f = rhoR;
            }
            else {
                if ( St >= xovert) {
                    /* Left of rarefaction tail */
                    u_f = us;
                    p_f = ps;
                    rho_f = rhoR*pow(ps/pR,1./g);
                }
                else {
                    /* In rarefaction fan */
                    p_f = 2./(g+1) - (gm/gp)*(uR - xovert)/aR;
                    rho_f = rhoR*pow(p_f, 2./(g-1));
                    p_f = pR*pow(p_f,1./gm);
                    u_f = (2./(g+1))*( -aR + .5*(g-1)*uR + xovert);
                }
            }

        }
    }


    /* Done */
    W[0] = rho_f;
    W[1] = u_f;
    W[2] = p_f;

    /* For passive scalars and normal velocities
     * we return whether or not to take 
     * the left or right state
     * True = left
     * False = right
     */
    return (us > xovert);

}

