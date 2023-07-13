#!/bin/bash
#SBATCH -J Save # Job Name
#SBATCH -o Save.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p normal   # Queue (partition) name -- normal, development, etc.
#SBATCH -t 5:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A GKIMP # Project name

module load python3
python3 plotmaker2.py
