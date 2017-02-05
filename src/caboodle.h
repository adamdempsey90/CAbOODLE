#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct Field {

    double *cpu;
    double *gpu;
    int size_x, size_y, size_z, size;

} Field;



Field *dens;
Field *pres;
Field *vy;
Field *vx;
Field *vz;

