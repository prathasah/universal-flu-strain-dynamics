#!/bin/bash
#SBATCH --job-name=cb_1
#SBATCH --ntasks=120 
#SBATCH --output cb_1.txt
#SBATCH --time=240:00:00
#SBATCH --mem-per-cpu=10000

module load Python
mpirun python optimize_cluster.py