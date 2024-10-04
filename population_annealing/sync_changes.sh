#!/bin/bash

# Define variables
LOCAL_DIR="/Users/shanekeiser/Documents/ANNNI/populationannealing"
REPO_DIR="/Users/shanekeiser/Documents/ANNNI/populationannealing/git-clone/stat-physics" # Update this path to where you cloned the repository

# Function to handle errors
handle_error() {
    echo "Error on line $1"
    exit 1
}

# Trap errors and pass the line number to the error handler
trap 'handle_error $LINENO' ERR

echo "Starting synchronization protocol..."

# Check if local directory exists
if [ ! -d "$LOCAL_DIR" ]; then
    echo "Local directory $LOCAL_DIR does not exist."
    exit 1
fi

# Check if repository directory exists
if [ ! -d "$REPO_DIR" ]; then
    echo "Repository directory $REPO_DIR does not exist."
    exit 1
fi

echo "Navigating to the repository directory..."
cd "$REPO_DIR" || exit 1

echo "Copying changes from local directory to the repository subfolder..."
rsync -av --exclude 'data' "$LOCAL_DIR"/ population_annealing/

echo "Adding changes to git..."
git add population_annealing/*

echo "Committing changes..."
git commit -m "Update population annealing files"

echo "Pushing changes to the remote repository..."
git push origin main

echo "Synchronization complete!"