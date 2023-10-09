#!/bin/bash
#SBATCH -J DNAMHD       # Job Name
#SBATCH -o DNAHD.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1          # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p icx-normal   # Queue (partition) name -- normal, development, etc.
#SBATCH -t 30:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A GKIMP # Project name

ibrun /work2/08929/echansen/stampede2/DNA2MHD/bin/dna
