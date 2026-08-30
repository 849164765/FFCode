#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
source /public3/soft/modules/module.sh
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make
# 512x64x1024 / BlockCellLen=32 -> 16x2x32=1024 blocks, 用 128 进程（8 block/rank）
mpiexec -n 128 ./khMag3dMaxwell.exe
