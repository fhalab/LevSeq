#!/bin/bash
# runminion.sh 
source $(conda info --base)/etc/profile.d/conda.sh

source ./helper.sh

# Activate minION environment for basecalling and demultiplexing
conda activate minion

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <experiment_name> <folder>"
    exit 1
fi

# Arguments

experiment_name="$1"
folder="$2"

# Tmp path to miION
minION_path="/home/emre/github_repo/minION/minION"

# Create Experiment Path 
experiment_path=$(join_paths $minION_path $experiment_name)

echo "Experiment path: $experiment_path"

