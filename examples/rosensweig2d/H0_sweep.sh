#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 64
# H0 sweep for the peak-height/valley-depth figure (paper Fig. 23 style).
# Each case runs in its own directory; post-process with:
#   python3 plot_peak_vs_H0.py
source /public3/soft/modules/module.sh
module unload gcc 2>/dev/null
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make

# field strengths (kA/m), 0 .. 15
H0_LIST="0 2 4 6 8 10 12 14 15"
# steps per case: peaks saturate around ~0.454mm (11 cells); 20000 steps is a
# good balance. Increase for stronger fields if still growing at the end.
TOTAL_STEP=20000

for H0 in ${H0_LIST}; do
  DIR="case_H0_${H0}"
  mkdir -p ${DIR}
  cd ${DIR}
  mpiexec -n 128 ../rosensweig2d.exe ../rosensweig2d.ini ${H0} ${TOTAL_STEP}
  cd ..
done
