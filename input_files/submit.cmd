#!/bin/bash
#SBATCH -J DNAMHD
#SBATCH --output=%x.out%j 
#SBATCH --qos=regular
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --constraint=cpu
#SBATCH -A m2116                                                                                                                           
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

srun -n 16 /global/homes/e/echansen/DNA2MHD/bin2/dna
