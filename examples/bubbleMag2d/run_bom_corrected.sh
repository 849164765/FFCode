#!/bin/bash
cd /home/dlf/myCode/Claude0714/examples/bubbleMag2d
make clean 2>/dev/null; make 2>&1 | tail -1
for Bom in 0.13 0.27 0.49 0.76 1.94; do
  echo "=== Running Bom=$Bom (corrected: auto H0, mu_h=2) ==="
  cat > bubbleMag2d.ini << EOF
[workdir]
workdir_ = ./
[parallel]
thread_num = 1
[Mesh]
Ni = 128
Nj = 256
Cell_Len = 1.0
BlockCellLen = 64
[Bubble]
Radius = 25
CenterX = 64
CenterY = 64
[Phase_Field]
Interface_Width = 4.0
Mobility = 0.01
[Two_Phase]
rho_l = 0.001
rho_h = 1.0
Eo = 20
Re = 40
U_g = 0.01
[Magnetic_Field]
chi_l = 0.0
chi_h = 1.0
mu_l = 1.0
mu_h = 2.0
Bom = ${Bom}
[Simulation_Settings]
TotalStep = 3500
OutputStep = 3500
EOF
  rm -rf vtkoutput
  mpirun -n 1 ./bubbleMag2d.exe > run_corr_bom${Bom}.log 2>&1
  mv vtkoutput vtkoutput_corr_Bom${Bom}
  echo "Done Bom=$Bom"
done
echo "All corrected runs complete"
