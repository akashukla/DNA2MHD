#!/bin/bash
#SBATCH -J DNAMHD
#SBATCH --qos=regular
#SBATCH --time=08:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --constraint=cpu
#SBATCH -A m2116
#SBATCH --output=%x.out%j
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

module load intel
module load forge
srun -n 256 /global/homes/e/echansen/DNA2MHD/bin2/dna
