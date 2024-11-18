#!/bin/bash
#SBATCH -J DNAMHD
#SBATCH --output=%x.out%j 
#SBATCH --qos=regular
#SBATCH --time=360
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=cpu
#SBATCH -A m2116                                                                                                                           
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

srun -n 4 valgrind --num-callers=500 --log-file=vg.%p /global/homes/e/echansen/DNA2MHD/bin2/dna
