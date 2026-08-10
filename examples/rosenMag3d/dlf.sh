#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
source /public3/soft/modules/module.sh
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make
# 128^3 网格 / BlockCellLen=32 -> 4x4x4=64 blocks, 用 64 进程
mpiexec -n 64 ./rosenMag3d.exe
