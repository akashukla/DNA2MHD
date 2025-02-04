#!/bin/bash
#SBATCH -J DNAMHD
#SBATCH --qos=regular
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --constraint=cpu
#SBATCH -A m2116
#SBATCH --output=%x.out%j
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

srun -n 4 /global/homes/e/echansen/DNA2MHD/bin2/dna
