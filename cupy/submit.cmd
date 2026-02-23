#!/bin/bash
#SBATCH --qos=regular
#SBATCH --time=03:00:00
#SBATCH --constraint=gpu
#SBATCH --gpus=1
#SBATCH -A m2116
#SBATCH --output=%x.out%j
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

module load nvidia
module load forge
srun -n 1 python3 exporder.py "$1"
