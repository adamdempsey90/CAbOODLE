#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define TRUE 1
#define FALSE 0

#define indx1 (i)
#define indx2 (j + i*(params->nx))
#define indx3 (k + (params->ny)*indx2)

#define indxm1 (i)
#define indxm2 (j + i*(params->nx+1))
#define indxm3 (k + (params->ny+1)*indxm2)

typedef struct Grid {
    int nx,ny,nz,size;
    double *xmin, *ymin, *zmin, *xmed, *ymed, *zmed;
    double *dx, *dy, *dz;
    double *Vol, *invVol, *Sxy, *Sxz, *Syz, *dSx, *dSy, *dSz;
} Grid;


typedef struct Field {
    int nx,ny,nz,size,nscalers;
    double *Density, *Vx, *Vy, *Vz, *Pressure, *Energy, *Scalars;
    double *Mx, *My, *Mz;
    double *Fx, *Fy, *Fz;
} Field;

typedef struct Parameters {
    int nx,ny,nz;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    double tend;
    int noutputs;
    int logx, logy, logz; 
    int nscalars;
    char outputdir[1024];

} Parameters;



Parameters *params;
Field *fld;
Grid *grid;



