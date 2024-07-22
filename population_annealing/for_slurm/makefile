# Include directory
IDIR = ./include
ODIR=.
# LDIR = /home/skeiser_umass_edu/populationannealing/lib
# Use module paths for GSL and LLVM
GSL_DIR = ./lib/gsl# Update this path based on the module system
LLVM_DIR = ./lib/llvm-project# Update this path based on the module system

# Use the correct compiler
CC = gcc  # or g++ or clang++

# Compiler and linker flags
CFLAGS = -I$(IDIR) -I$(GSL_DIR)/include -I$(LLVM_DIR)/llvm/include -I$(LLVM_DIR)/clang/include -fopenmp -O3 
LFLAGS = -L$(GSL_DIR)/lib -L$(LLVM_DIR)/llvm/lib -L$(LLVM_DIR)/clang/lib -L$(LLVM_DIR)/openmp/lib -fopenmp -O3 

# Libraries
LIBS = -lm -lgsl -lgslcblas -lstdc++

# Dependencies
_DEPS = population_class.h lattice_class.h spin_class.h functions.h parameters.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Object files
_OBJ = main.o population_class.o lattice_class.o spin_class.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Compile object files
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CFLAGS)

# Link the final executable
main: $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean

# Clean object files and executable
clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~