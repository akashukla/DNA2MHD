#!/bin/bash
#SBATCH -J DNAMHD       # Job Name
#SBATCH -o DNAHD.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1          # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p vm-small   # Queue (partition) name -- normal, development, etc.
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -t 24:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A PHY23037 # Project name
#SBATCH --mail-user=echansen@tacc.utexas.edu

ibrun /work2/08929/echansen/ls6/DNA2MHD/bin2/dna
