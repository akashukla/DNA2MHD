#!/bin/bash
#SBATCH --qos=preempt
#SBATCH --time=06:00:00
#SBATCH --constraint=gpu
#SBATCH --gpus=1
#SBATCH -A m2116
#SBATCH --output=%x.out%j
#SBATCH --mail-type=begin,end,fail
#SBATCH --mail-user=ehansen99@utexas.edu

module load intel
module load forge
module load python
module list
python3 -m pip list

srun -n 1 python3 restartnlinvtest.py "$1"
