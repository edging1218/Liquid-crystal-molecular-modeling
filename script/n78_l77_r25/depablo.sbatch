#!/bin/sh

#SBATCH --job-name=78_77_25
#SBATCH --partition=depablo
#SBATCH --qos=depablo
#SBATCH --output=drop.out
#SBATCH --error=drop.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:0:0

module load intelmpi/5.1+intel-16.0

mpirun -np 16 /scratch/midway2/yezhou/triple/iso

