#!/usr/bin/env bash

# Run run_optical_model.py -p P  with P = 0..199

# Exit immediately if any command fails
set -e

echo "Running with conda env $CONDA_DEFAULT_ENV"

# Create directory if it does not exist
if [ ! -d p0 ]; then
    mkdir -p p0
fi

# Create directory if it does not exist
if [ ! -d indicator ]; then
    mkdir -p indicator
fi

for P in $(seq 0 199); do
    CMD="mpirun -n 32 python -u run_optical_model.py -p $P"

    # Print the command before executing
    echo "$CMD"

    # Execute the command
    eval "$CMD"
done