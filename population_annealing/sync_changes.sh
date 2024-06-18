#!/bin/bash

# Define variables
LOCAL_DIR="/Users/shanekeiser/Documents/Summer 2024/Research/PopulationAnnealing"
REPO_DIR="/Users/shanekeiser/Documents/Summer 2024/Research/stat-physics" # Update this path to where you cloned the repository

# Navigate to the local repository
cd "$REPO_DIR"

# Copy the changes from the local directory to the repository subfolder
cp -R "$LOCAL_DIR"/* population_annealing/

# Add changes to git
git add population_annealing/*

# Commit the changes with a message
git commit -m "Update population annealing files"

# Push the changes to the remote repository
git push origin main