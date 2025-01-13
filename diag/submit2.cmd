#!/bin/bash
#SBATCH -J Save
#SBATCH --output=%x.out%j 
#SBATCH --qos=regular
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=cpu
#SBATCH -A m2116                                                            
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

srun -n 1 plotmaker3.py
