# kh2dStd — 2D Kelvin-Helmholtz 标准对照案例

论文 Sec.F：域 `[0,2L]x[0,2L]`, `L=256`, `Re=5000, We=10000, Fr=1`,
`U=0.02, rho_l/rho_h=0.99, W=4, M=0.01`；磁场 `chi=(0,1), mu=(1,2)`，
`Bom={0,636,1434,2547,3980}`。

## 修改内容（只改本案例，src 不动）

1. **KH 专用磁压力项**
   `src` 的 Kelvin 力 `F_k = chi*|H|*grad|H|` 对 bubble/rosen 法向场正确，
   保持不动。本案例追加与 `khMag3d`/`khMagDeepSeek3d` 相同的水平切向场项：
   ```
   F_p = MagPressure_Factor * (-0.5 * Hx^2 * grad(mu))
   grad(mu) = (mu_h - mu_l) * grad(phi)
   ```
   ini 中 `MagPressure=1` 开启；`MagPressure_Factor` 用于对文献定量校准。

2. **壁面磁力过滤**
   上下壁 Dirichlet 磁势会在最靠壁的 ghost 处产生虚假 `F_k ~ 0.5*H0^2` 尖峰
   （与界面无关）。磁力与诊断跳过上下壁各 16 格子（不超过 10% 域高）；
   左右缝 ghost 由第 4 条的手动同步修复，不再丢失缝上磁力。

3. **相位修正（重要）**
   文献 Fig.30/31 中域内是两个完整波峰 + 中间一个完整波谷 + 左右各半个
   波谷，即 `y = L - 0.1L*cos(2πx/L)`。旧代码用了 `+cos`（中间波峰/两侧
   波谷），已改为文献相位。

4. **左右周期缝 ghost 同步（边界伪影修复）**
   `MF_Per.Apply()` 只交换 POP，不交换 HX/HY/HMAG 的卷绕 ghost 列，
   缝上 `∇|H|` 会产生虚假尖峰（左右边界伪影，Rosen 同款根因）。已移植
   rosenMag3d/khMag3d 的 `MFGhostSyncPlans`，在每次 H 计算后同步
   HX/HY/HMAG 缝列；PSI 缝仍由扭曲周期 BC 单独处理。

5. **时间标定可配置 + 双时间输出**
   `KH_TSTAR_SCALE` 可在 ini 中修改，日志同时输出：
   ```
   [KH2D_DBG] step=... t*_eq=... t*_plot=... amp=... ncross=...
   ```
   `t*_eq=step*sqrt(g/L)` 是 Eq.(84) 原定义；`t*_plot=KH_TSTAR_SCALE*t*_eq`
   是与论文图像对齐的时间。

6. **全局磁力诊断**
   `[KH2D_MAG]` 输出全域 MPI 归约后的 Hmag 和
   `F_kelvin / F_press / F_total` 的最大值与平均值。

## 重跑 Bom 扫描（集群）

```bash
sbatch dlf.sh            # mpiexec -n 128, 512x512, 30000 步
```

每个 Bom 先改 ini 再提交：
```bash
sed -i 's/^Bom = .*/Bom = 2547/' kh2dStd.ini
```

## 日志分析

```bash
python3 analyze_logs.py .            # 旧日志 (out0.log 等)
python3 analyze_logs.py . --k 3.0    # 指定 t* 标定系数
python3 analyze_logs.py . --anchor 2 8000 9000
```

时间标定（2026-08-20 更新）：
- 用户复核：Bom=0 时 step=25600（旧代码自报 t*=6）实际对应文献 **t*=8**，
  故 `KH_TSTAR_SCALE=4.0`；t*={2,4,6,8} ↔ step={6400,12800,19200,25600}。
- 旧 `k=3` 是相位没对齐时的视觉假象；修正相位后再按 Bom=0 复核。
磁力：
- 旧 `k=3` 下 Bom=3980 step=10000 看起来像文献 Bom=1434 t*=4，
  等效缺 `3980/1434≈2.78` 倍力系数，这就是新增磁压力项的原因。
- **磁势子迭代 + 磁压力系数（已校准）**：2D 校准确认
  `PsiSolver_Iter = 100`、`MagPressure_Factor = 1.5` 时图像与文献基本对齐。
  3D 已同步为相同参数。
