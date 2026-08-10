#!/bin/bash
#SBATCH -p amd_512
#SBATCH -N 1
#SBATCH -n 128
set -euo pipefail

source /public3/soft/modules/module.sh
module load mpich/3.1.4-gcc8.1.0
module load gcc/12.2

make clean
make

H0_LIST="0 2 4 4.9 6 8 8.2 10 12 14 15"
TOTAL_STEP=20000

for H0 in ${H0_LIST}; do
  DIR="case_H0_${H0}"
  mkdir -p "${DIR}"
  cd "${DIR}"
  mpiexec -n 128 ../rosensweig2d.exe ../rosensweig2d.ini "${H0}" "${TOTAL_STEP}"
  cd ..
done

# One compact time-series file, generated on the cluster without copying VTI.
# Block count is auto-detected per case (works for any MPI rank count);
# 64 workers parallelize over the tens of thousands of VTI files.
python3 extract_rosensweig_stats.py --workers 64
