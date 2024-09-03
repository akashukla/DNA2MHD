#!/bin/bash
#SBATCH -J DNAMHD       # Job Name
#SBATCH -o DNAHD.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1          # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p development   # Queue (partition) name -- normal, development, etc.
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -t 00:20:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A PHY23037 # Project name
#SBATCH --mail-user=echansen@tacc.utexas.edu

ibrun ../bin2/dna
