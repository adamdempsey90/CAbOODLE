EXECUTABLE=caboodle
SOURCES=read_params.c metric.c output.c init.c algogas.c main.c 
HEADER=caboodle.h

OPENMPLIB=-lgomp
MATHLIB=-lm

LDFLAGS=$(OPENMPLIB) $(MATHLIB)

CFLAGS=-c -fopenmp -Wall -O3 -g 


INCLIB=
LDLIB=


BIN=bin/
SRC=src/


UNAME := $(shell echo $(USER))

ifeq ($(UNAME),apollo)
CC=gcc
endif

ifeq ($(UNAME),jupiter)
CC=gcc-4.9
endif
ifeq ($(UNAME),zeus)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif


ifeq ($(UNAME),helios)
CC=gcc-4.9
INCLIB=-I/usr/local/include/
LDLIB=-L/usr/local/lib/
endif

ifeq ($(UNAME),amd616)
CC=gcc
LDLIB=-L/software/lapack/3.4.0/lib -L/software/gsl/1.16-gcc4.8.3/lib/ -L/software/hdf5/1.8.12-serial/lib/
INCLIB=-I/software/hdf5/1.8.12-serial/include/
endif

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(COBJECTS)
	$(CC)  $(COBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@


clean:
	rm $(COBJECTS) $(EXECUTABLE)
