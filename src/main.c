#include "caboodle.h"


extern void init_params(int,char **);
extern void init_grid(void);
extern void init_field(void);

extern void algogas(double *, double);
extern void output(int,double);
extern void free_field(void);
extern void free_grid(void);
extern void free_params(void);

int main(int argc, char *argv[]) {
    int i;



    init_params(argc, argv);
    printf("Init grid\n");
    init_grid();
    printf("Init field\n");
    init_field();

    double t = 0;
    double tend = params->tend;;
    double dt_outputs = tend/(params->noutputs);


    for(i=0;i<params->noutputs;i++) {
        tend = t + dt_outputs;
        printf("Algogas\n");
        algogas(&t, tend-t);
        printf("output\n");
        output(i,t);
    }
    
    printf("Done\n");

    free_field();
    free_grid();
    free_params();
    return 1;
}
