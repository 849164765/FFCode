#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/modules/module.sh
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make
# 128x128x512 网格 / BlockCellLen=32 -> 4x4x16=256 blocks, 用 128 进程 (每进程 2 块)
mpiexec -n 128 ./rtMag3d.exe
