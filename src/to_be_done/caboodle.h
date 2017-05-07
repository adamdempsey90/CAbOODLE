#include <stdlib.h>
#include <stdio.h>
#include <math.h>

enum{DDD,EEE,UUU,VVV,WWW};

typedef struct Field {

    double *cons_cpu;
    double *prim_cpu;
    double *Ql_cpu;
    double *Qr_cpu;
    double *Fl_cpu;
    double *Fr_cpu;
    double *cons_gpu;
    double *prim_gpu;
    double *Ql_gpu;
    double *Qr_gpu;
    double *Fl_gpu;
    double *Fr_gpu;
    int size_x, size_y, size_z, size;
    int nfields;
    int nscalars;
} Field;

typedef struct Parameters {
    int nx,ny,nz;
    double xmin,ymin,zmin;
    double xmax,ymax,zmax;
} Parameters;

typedef struct Grid {
    double *xm_cpu;
    double *ym_cpu;
    double *zm_cpu;
    double *xmed_cpu;
    double *ymed_cpu;
    double *zmed_cpu;
    double *dx_cpu;
    double *dy_cpu;
    double *dz_cpu;
    double *dAxy_cpu;
    double *dAxz_cpu;
    double *dAyz_cpu;
    double *dV_cpu;
    double *xm_gpu;
    double *ym_gpu;
    double *zm_gpu;
    double *xmed_gpu;
    double *ymed_gpu;
    double *zmed_gpu;
    double *dx_gpu;
    double *dy_gpu;
    double *dz_gpu;
    double *dAxy_gpu;
    double *dAxz_gpu;
    double *dAyz_gpu;
    double *dV_gpu;
} Grid;



Field *Domain;

