#!/bin/bash
#SBATCH -J Save # Job Name
#SBATCH -o Save.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p normal   # Queue (partition) name -- normal, development, etc.
#SBATCH -t 1:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A GKIMP # Project name

module load python3
python3 /work2/04943/akshukla/stampede2/DNA2MHD/diag/dna2mhd_utils_exe.py '/scratch/04943/akshukla/dna2mhd_output_zero_div'

