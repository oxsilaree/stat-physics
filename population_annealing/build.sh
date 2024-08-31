#!/bin/bash

# Quick fix because if we don't do this, apparently libgsl.so.28 cannot be found.
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/skeiser_umass_edu/populationannealing/lib/gsl/lib
export LD_LIBRARY_PATH




# Run make command
make



# Remove .o files
make clean



# Other SLURM Commands
# salloc -p cpu -c 4 --mem=8GB

# srun --output=./slurm_history/slurm_shanekeiser_%j.out --error=./slurm_history/slurm_shanekeiser_%j.err --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8GB ./main 0



