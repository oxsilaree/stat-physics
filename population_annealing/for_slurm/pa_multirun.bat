#!/bin/bash

# sbatch instructions
#SBATCH -c 10                           # number of CPU cores per task
#SBATCH --ntasks=5                      # number of tasks
#SBATCH --mem=4GB                       # requested memory 8GB
#SBATCH -t 04:00:00                     # maximum job run-time
#SBATCH -J skeiser_PA_job                  # job name
#SBATCH -o ./slurm_history/slurm_shanekeiser_%j.out     # output file
#SBATCH --mail-type=BEGIN
/bin/true

module purge                            # removes any existing modules

echo Starting tasks...

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/skeiser_umass_edu/populationannealing/lib/gsl/lib
export LD_LIBRARY_PATH
export OMP_NUM_THREADS=500

srun --exclusive --ntasks=1 ./main ${SLURM_ARRAY_TASK_ID}

echo Done!