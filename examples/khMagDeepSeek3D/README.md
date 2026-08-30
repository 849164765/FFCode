# khMagDeepSeek3D

独立重建的 3D Kelvin-Helmholtz 铁磁流体案例。

## 关键设计

1. **磁力模型（仅 KH 案例内补充，src 不变）**
   `src/mfield` 的 Kelvin 力对 bubbleMag3d/rosenMag3d 等法向场案例是正确的，保持不动。
   KH 是水平切向场，本案例在源 Kelvin 力之后补充 KH 专用的 Maxwell 磁压力：

   ```text
   F_p = MagPressure_Factor * (-0.5 * Hx^2 * grad(mu))
   ```
   2D 校准取值 `MagPressure_Factor = 1.5`。

   其中 `grad(mu) = (mu_h - mu_l) * grad(phi)`。

2. **磁势求解**
   使用有限差分 Gauss-Seidel 求解 `div(mu grad psi)=0`，面导磁率取两侧 mu 的调和平均。
   不使用 D3Q7 MRT 伪时间迭代：水平场在 x 方向斜坡+周期缝上会漂移到反号 delta_psi。
   每流体步 `PsiSolver_Iter = 100` 次扫描（与 2D 校准一致）。

3. **解析磁势初值**
   `psi = -H0*x + delta_psi`，delta_psi 为小振幅正弦界面的有限高度解析解。
   界面相位已修正为文献式 `h(x) = z0 - A*cos(2πx/L)`。

4. **两相参数与稳定钳位**
   铁磁流体(φ=1)在下 rho_h=1.0，溶剂(φ=0)在上 rho_l=0.99。
   eta_l 默认钳到 omega=1.98，避免 L=64 下 Re=5000 的 omega~1.997 不稳定加速。

## 时间标定
`t*_plot = 4 * t*_eq = 4*step*sqrt(g/L_ref) = step/800`（与 2D `KH_TSTAR_SCALE=4` 一致）。

## 编译运行

```bash
make
mpiexec -n 64 ./khMagDeepSeek3d.exe khMagDeepSeek3d.ini
# 或集群
sbatch dlf.sh
```

Bom 扫描：修改 ini 中 `Bo_m = 0 / 636 / 1434 / 2547 / 3980`。
