IDIR =./include
CC=g++

#Use these flags if you don't want to use optimization profiles
CFLAGS=-I$(IDIR) -fopenmp -O3 

# Compile with -fprofile-generate and run with small pop size (R~100).\
This creates optimization profiles. Use second line to use profiles in\
real production run. Speedup is approximately a factor of 2.

#CFLAGS=-I$(IDIR) -fopenmp -O3 -fprofile-generate=./profiles
#CFLAGS=-I$(IDIR) -fopenmp -O3 -fprofile-use=./profiles

ODIR=.
LDIR =.

LIBS=-lm -lgsl -lgslcblas -lstdc++

_DEPS = pa_class.h ising_class.h functions.h params.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o pa_class.o ising_class.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
