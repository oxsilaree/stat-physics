#!/bin/bash

## export LLVM_DIR=$(brew --prefix llvm)
## export GSL_DIR=$(brew --prefix gsl)

# Set environment variables for LLVM and OpenMP
export PATH="/usr/local/opt/llvm/bin:$PATH"
export MACOSX_DEPLOYMENT_TARGET=12.7.5
## export LDFLAGS="-L$LLVM_DIR/lib -L$GSL_DIR/lib"
## export CPPFLAGS="-I$LLVM_DIR/include -I$GSL_DIR/include"
## export CFLAGS="-I$LLVM_DIR/include -I$GSL_DIR/include"

# Run make command
make

# Remove .o files
# make clean