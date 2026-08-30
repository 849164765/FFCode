#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/modules/module.sh
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make
mpiexec -n 128 ./RTDeepSeek2D.exe
