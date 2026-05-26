#!/bin/bash

# Check if the environment already exists
if conda info --envs | grep -q "viraldb"; then
    echo "Conda environment 'viraldb' already exists. Skipping creation..."
else
    # Create the conda environment from environment.yml
    echo "Creating the conda environment from environment.yml..."
    conda env create -f environment.yml
fi

# Activate the conda environment
echo "Activating the conda environment..."
source ~/opt/anaconda3/etc/profile.d/conda.sh  # Ensure conda is properly sourced
conda activate viraldb

# Run the Python script to import CSVs into the SQLite database
echo "Running the Python script..."
python FIRST_Reassemble.py

echo "Setup complete."
