#include "caboodle.h"
#include <ctype.h>

#define PRINT_DOUBLE(NAME,VAL) printf("\t%s = %lg\n",NAME,VAL)
#define PRINT_INT(NAME,VAL) printf("\t%s = %d\n",NAME,VAL)
#define PRINT_STR(NAME,VAL) printf("\t%s = %s\n",NAME,VAL)
#define FPRINT_DOUBLE(F,NAME,VAL) fprintf(f,"%s = %lg\n",NAME,VAL)
#define FPRINT_INT(F,NAME,VAL) fprintf(f,"%s = %d\n",NAME,VAL)
#define FPRINT_STR(F,NAME,VAL) fprintf(f,"%s = %s\n",NAME,VAL)

void set_var(char *name,int int_val, double double_val, int bool_val, char *str_val) {
    // Ints
    if (strcmp(name,"nx") == 0) {
        params->nx = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"ny") == 0) {
        params->ny = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"nz") == 0) {
        params->nz = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"noutputs") == 0) {
        params->noutputs = int_val;
        PRINT_INT(name,int_val);
    }
    else if (strcmp(name,"nscalars") == 0) {
        params->nscalars = int_val;
        PRINT_INT(name,int_val);
    }
    // Doubles
    else if (strcmp(name,"tend") == 0) {	
        params->tend = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"xmin") == 0) {	
        params->xmin = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"xmax") == 0) {	
        params->xmax = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"ymin") == 0) {	
        params->ymin = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"ymax") == 0) {	
        params->ymax = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"zmin") == 0) {	
        params->zmin = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    else if (strcmp(name,"zmax") == 0) {	
        params->zmax = double_val;
        PRINT_DOUBLE(name,double_val);

    }
    // Bools
    else if (strcmp(name,"logx") == 0) {	
        params->logx = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"logy") == 0) {	
        params->logy = bool_val;
        PRINT_STR(name,str_val);

    }
    else if (strcmp(name,"logz") == 0) {	
        params->logz = bool_val;
        PRINT_STR(name,str_val);

    }
    // Strings
    else if (strcmp(name,"outputdir") == 0) {	
        sprintf(params->outputdir,"%s",str_val);
        PRINT_STR(name,str_val);
    }

    return;
}

void parse_argument(int argc, char *argv[]) {
    int j;
    unsigned int i;
    char name[100],strval[100];
    double dval;
    int ival;
    int bool_val;
    char testbool;


    for(j=0;j<argc;j++) {
        sscanf(argv[j],"%32[^=]=%s",name,strval);
        dval = atof(strval);
        ival = atoi(strval);
        testbool = toupper(strval[0]);
        if (testbool == 'Y') bool_val = TRUE;
        else bool_val = FALSE;
        for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
        set_var(name,ival,dval,bool_val,strval);
    }



    return;
}

void read_param_file(char *fname, int argc, char *argv[]) {
    FILE *f;

    char tok[20] = "\t :=>";

    char line[100],name[100],strval[100];
    char *data;
    double temp;
    int status;
    int int_val;
    int bool_val;
    char testbool;
    unsigned int i;

    f= fopen(fname,"r");

    while (fgets(line,100,f)) {
       // printf("%s\n",line);
        status = sscanf(line,"%s",name);

      //  printf("%s\n",name);
        if (name[0] != '#' && status == 1) {
        
             data = line + (int)strlen(name);
             sscanf(data + strspn(data,tok),"%lf",&temp);
             sscanf(data + strspn(data,tok),"%s",strval);
             //printf("%lf\t%s\n",temp,strval);
            int_val = (int)temp;
            testbool = toupper(strval[0]);
            if (testbool == 'Y') bool_val = TRUE;
            else bool_val = FALSE;
            
            for (i = 0; i<strlen(name); i++) name[i] = (char)tolower(name[i]);
            
            set_var(name,int_val,temp,bool_val,strval);

        }
    }

    
    if (argc > 0) {
        printf("Redefined on the command line:\n");
        parse_argument(argc,argv);
    }



    return;
}


