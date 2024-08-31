IDIR = ./include
ODIR = .
LLVM_DIR = $(shell brew --prefix llvm)
OMP_DIR = $(shell brew --prefix libomp)
CC = /opt/homebrew/opt/llvm/bin/clang++

# Compiler and linker flags
CPPFLAGS = -I$(IDIR) -I/usr/local/include -I$(LLVM_DIR)/include -I$(OMP_DIR)/include -fopenmp -O3
LDFLAGS = -L/usr/local/lib -L$(LLVM_DIR)/lib -L$(OMP_DIR)/lib -fopenmp

# Libraries
LIBS = -lm -lgsl -lgslcblas -lstdc++ -lfftw3

# Dependencies
_DEPS = population_class.h lattice_class.h spin_class.h functions.h parameters.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Object files
_OBJ = main.o population_class.o lattice_class.o spin_class.o functions.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Compile object files
$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -g -c -o $@ $< $(CPPFLAGS) -fsanitize=address -fno-omit-frame-pointer -Wno-format-security

# Link the final executable
main: $(OBJ)
	$(CC) -g -o $@ $^ $(LDFLAGS) $(LIBS) -fsanitize=address

.PHONY: clean

# Clean object files and executable
clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~