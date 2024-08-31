#!/bin/bash

# Define variables
LOCAL_DIR="./for_slurm/*"   # Path to the local directory (use * to copy contents)

# Copy the contents of the local directory to the remote directory
scp -r $LOCAL_DIR unity:~/populationannealing

# Print a message indicating that the operation is complete
echo "Directory contents copied to UNITY"