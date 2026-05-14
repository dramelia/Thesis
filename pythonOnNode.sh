#!/bin/bash
#SBATCH --job-name=GSEoverlaps
#SBATCH --time=05:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

echo "Running on `hostname`"

python /PATH/SCRIPT.py # change path and script

echo "Finished Job"