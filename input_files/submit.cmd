#!/bin/bash
#SBATCH -J DNAMHD       # Job Name
#SBATCH -o DNAHD.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1          # Total number of mpi tasks requested
#SBATCH -N 1           # Total number of mpi tasks requested
#SBATCH -p vm-small   # Queue (partition) name -- normal, development, etc.
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -t 10:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A PHY23037 # Project name
#SBATCH --mail-user=echansen@tacc.utexas.edu

python3 runner.py
cd $WORK/DNA2MHD/diag/
python3 convg711.py
