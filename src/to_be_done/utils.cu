#include "caboodle.h"
#include <cuda.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
    if (code != cudaSuccess) {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}
extern "C" void dev2host(Field *fld) {
    gpuErrchk( cudaMemcpy( fld->cpu, fld->gpu, sizeof(double)*(fld->size),cudaMemcpyDeviceToHost)) 
}
extern "C" void host2dev(Field *fld) {
    gpuErrchk( cudaMemcpy( fld->gpu, fld->cpu, sizeof(double)*(fld->size),cudaMemcpyHostToDevice)) 
}

extern "C" void all2host(void) {
    dev2host(dens);
    dev2host(pres);
    dev2host(vx);
    dev2host(vy);
    dev2host(vz);
}
extern "C" void all2dev(void) {
    host2dev(dens);
    host2dev(pres);
    host2dev(vx);
    host2dev(vy);
    host2dev(vz);
}

extern "C" void reset_field(Field *fld) {
    int i;
    for(i=0;i<fld->size;i++) fld->cpu[i] = 0;
    host2dev(fld);
}

extern "C" Field *init_field(int size_x, int size_y, int size_z) {
    Field *fld;
    double *arr_cpu;
    double *arr_gpu;
    fld = (Field *)malloc(sizeof(Field));

    int size = size_x*size_y*size_z;

    arr_cpu = (double *)malloc(sizeof(double)*size);
    if (arr_cpu == NULL) {
        printf("Not enough space on CPU\n");
        exit(1);
    }
    /* Should use cudaMallocPitch instead */
    gpuErrchk(cudaMalloc((void **) &arr_gpu, size_x*size_y*size_z*sizeof(double)))
    if (arr_gpu == NULL) {
        printf("Not enough space on GPU\n");
        exit(1);
    }

    fld->size_x = size_x;
    fld->size_y = size_y;
    fld->size_z = size_z;
    fld->size = size;
    fld->cpu = arr_cpu;
    fld->gpu = arr_gpu;
    reset_field(fld);
    return fld;
}
extern "C" void free_field(Field *fld) {
    free(fld->cpu);
    gpuErrchk(cudaFree(fld->gpu));
    free(fld);
}
