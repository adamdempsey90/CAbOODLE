#include "caboodle.h"


extern double scale_dLx(double, double, double, double, double, double);
extern double scale_dLy(double, double, double, double, double, double);
extern double scale_dLz(double, double, double, double, double, double);
extern double scale_Sxy(double, double, double, double, double, double);
extern double scale_Sxz(double, double, double, double, double, double);
extern double scale_Syz(double, double, double, double, double, double);
extern double scale_vol(double, double, double, double, double, double);

extern void read_param_file(char *, int, char **);

void init_params(int argc, char *argv[]) {
/* Read in the parameters and fill params */
    params = (Parameters *)malloc(sizeof(Parameters));
    char parfile[256];
    if (argc < 2) {
        strcpy(parfile,"in/in.par");
    }
    else {
        strcpy(parfile,argv[1]);
    }
    printf("Reading parameters from %s...\n",parfile);

    if (argc > 2) {
        read_param_file(parfile,argc-2,&argv[2]);
    }
    else {
        read_param_file(parfile,0,NULL);
    }
}

void init_grid(void) {
/* Allocate and set up the grid */
    int  nx = params->nx;
    int  ny = params->ny;
    int  nz = params->nz;
    int  size = nx*ny*nz;

    grid = (Grid *)malloc(sizeof(Grid));
    
    grid->xmin = (double *)malloc(sizeof(double)*(nx+1));
    grid->ymin = (double *)malloc(sizeof(double)*(ny+1));
    grid->zmin = (double *)malloc(sizeof(double)*(nz+1));

    grid->xmed = (double *)malloc(sizeof(double)*nx);
    grid->ymed = (double *)malloc(sizeof(double)*ny);
    grid->zmed = (double *)malloc(sizeof(double)*nz);
    grid->dx = (double *)malloc(sizeof(double)*nx);
    grid->dy = (double *)malloc(sizeof(double)*ny);
    grid->dz = (double *)malloc(sizeof(double)*nz);

    double *xmin = grid->xmin;
    double *ymin = grid->ymin;
    double *zmin = grid->zmin;
    double *xmed = grid->xmed;
    double *ymed = grid->ymed;
    double *zmed = grid->zmed;
    double *dx   = grid->dx ;
    double *dy   = grid->dy ;
    double *dz   = grid->dz ;

    int i,j,k;


    for(i=0;i<nx+1;i++) {
        if (params->logx) {
            xmin[i] = exp( log(params->xmin) + log(params->xmax/params->xmin)/nx); 
        }
        else {
            xmin[i] = params->xmin + (params->xmax - params->xmin)/nx; 
        }
    }
    for(i=0;i<nx;i++) {
        xmed[i] = (xmin[i]+ xmin[i+1])*.5;
        dx[i] = xmin[i+1]-xmin[i];
    }

    for(i=0;i<ny+1;i++) {
        if (params->logy) {
            ymin[i] = exp( log(params->ymin) + log(params->ymax/params->ymin)/ny); 
        }
        else {
            ymin[i] = params->ymin + (params->ymax - params->ymin)/ny; 
        }
    }
    for(i=0;i<ny;i++) {
        ymed[i] = (ymin[i]+ ymin[i+1])*.5;
        dy[i] = ymin[i+1]-ymin[i];
    }

    for(i=0;i<nz+1;i++) {
        if (params->logz) {
            zmin[i] = exp( log(params->zmin) + log(params->zmax/params->zmin)/nz); 
        }
        else {
            zmin[i] = params->zmin + (params->zmax - params->zmin)/nz; 
        }
    }
    for(i=0;i<nz;i++) {
        zmed[i] = (zmin[i]+ zmin[i+1])*.5;
        dz[i] = zmin[i+1]-zmin[i];
    }

    
    grid->Vol = (double *)malloc(sizeof(double)*size);
    grid->invVol = (double *)malloc(sizeof(double)*size);
    grid->Sxy = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));
    grid->Sxz = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));
    grid->Syz = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));
    grid->dSx = (double *)malloc(sizeof(double)*size);
    grid->dSy = (double *)malloc(sizeof(double)*size);
    grid->dSz = (double *)malloc(sizeof(double)*size);

    for(k=0;k<nz;k++) {
        for(j=0;j<ny;j++) {
            for(i=0;i<nx;i++) {
                grid->Vol[indx3] = scale_vol(xmed[i],ymed[j],zmed[k], dx[i],dy[j],dz[k]); 
                grid->invVol[indx3] = 1./grid->invVol[indx3];
                grid->dSx[indx3] = scale_dLx(xmed[i],ymed[j],zmed[k],dx[i],dy[j],dz[k]);
                grid->dSy[indx3] = scale_dLy(xmed[i],ymed[j],zmed[k],dx[i],dy[j],dz[k]);
                grid->dSz[indx3] = scale_dLz(xmed[i],ymed[j],zmed[k],dx[i],dy[j],dz[k]);
            }
        }
    }

    for(k=0;k<nz;k++) {
        for(j=0;j<ny;j++) {
            for(i=0;i<nx;i++) {
                grid->Sxy[indxm3] = scale_Sxy(xmin[i],ymin[j],zmin[k], dx[i],dy[j],dz[k]); 
                grid->Sxz[indxm3] = scale_Sxz(xmin[i],ymin[j],zmin[k], dx[i],dy[j],dz[k]); 
                grid->Syz[indxm3] = scale_Syz(xmin[i],ymin[j],zmin[k], dx[i],dy[j],dz[k]); 
            }
        }
    }


    grid->nx = nx;
    grid->ny = ny;
    grid->nz = nz;
    grid->size= nx*ny*nz;
    return;
}

void init_field(void) {
/* Allocate and set up all fluid variables */
    int  nx = params->nx;
    int  ny = params->ny;
    int  nz = params->nz;
    int  size = nx*ny*nz;

    double *dens = (double *)malloc(sizeof(double)*size);
    double *vx = (double *)malloc(sizeof(double)*size);
    double *vy = (double *)malloc(sizeof(double)*size);
    double *vz = (double *)malloc(sizeof(double)*size);
    double *pres = (double *)malloc(sizeof(double)*size);
    double *energy = (double *)malloc(sizeof(double)*size);
    double *mx = (double *)malloc(sizeof(double)*size);
    double *my = (double *)malloc(sizeof(double)*size);
    double *mz = (double *)malloc(sizeof(double)*size);
    double *scalars = (double *)malloc(sizeof(double)*size*(params->nscalars));
    double *fx = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));
    double *fy = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));
    double *fz = (double *)malloc(sizeof(double)*(nx+1)*(ny+1)*(nz+1));


    int i;
    for(i=0;i<size;i++) {
        dens[i] = 1.0;
        pres[i] = 1.0;
        vx[i] = 0;
        vy[i] = 0;
        vz[i] = 0;
        mx[i] = 0;
        my[i] = 0;
        mz[i] = 0;
        energy[i] = 1.0;
    }
    for(i=0;i<size*(params->nscalars);i++) {
        scalars[i] = 0;
    }
    for(i=0;i<(nx+1)*(ny+1)*(nz+1);i++) {
        fx[i] = 0;
        fy[i] = 0;
        fz[i] = 0;
    }
    fld = (Field *)malloc(sizeof(Field));


    fld->Density = dens;
    fld->Pressure = pres;
    fld->Vx = vx;
    fld->Vy = vy;
    fld->Vz = vz;
    fld->Energy = energy;
    fld->Mx = mx;
    fld->My = my;
    fld->Mz = mz;
    fld->Scalars = scalars;
    fld->Fx = fx;
    fld->Fy = fy;
    fld->Fz = fz;

    fld->nx = nx;
    fld->ny = ny;
    fld->nz = nz;
    fld->size =  size;

    return;
}

void free_params(void) {
    /* Free params struct */
    free(params);
}

void free_grid(void) {
    /* Free grid struct */
    free(grid->xmin  );
    free(grid->ymin );
    free(grid->zmin );
    free(grid->xmed );
    free(grid->ymed );
    free(grid->zmed );

    free(grid->Vol );
    free(grid->invVol );
    free(grid->Sxy );
    free(grid->Sxz );
    free(grid->Syz );
    free(grid->dSx );
    free(grid->dSy );
    free(grid->dSz );
}

void free_field(void) {
    /* Free field struct */

    free(fld->Density);
    free(fld->Pressure );
    free(fld->Vx );
    free(fld->Vy );
    free(fld->Vz );
    free(fld->Energy );
    free(fld->Mx );
    free(fld->My );
    free(fld->Mz );
    free(fld->Scalars );
    free(fld);

}
