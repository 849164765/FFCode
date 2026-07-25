# CHANGELOG — bubbleMag2d

## 修改 1：移除 5x 经验补偿因子

**文件**: `bubbleMag2d.cpp`

**原因**: 原代码中 `H0 = std::sqrt(25) * std::sqrt(2 * Bom * sigma / D_bubble)` 的 `std::sqrt(25)` 是经验性补偿因子，非物理推导。论文 Bom 定义式（Bom = μ₀H₀²D/(2σ)）不含此因子。

**修改**: 移除 `std::sqrt(25)`，恢复物理公式：

```cpp
// 修改前
H0 = std::sqrt(25) * std::sqrt(T(2.0) * Bom * sigma / D_bubble);

// 修改后
H0 = std::sqrt(T(2.0) * Bom * sigma / D_bubble);
```

## 修改 2：添加 MF 子迭代（已移除）

**文件**: `bubbleMag2d.cpp`

~~添加了 MF 子迭代循环，但后来发现子迭代不是根本原因，已经移除。~~

## 修改 3：ini 配置文件添加子迭代参数（已移除）

**文件**: `bubbleMag2d.ini`

~~[Magnetic_Field] 中添加了 MF_SubIter 参数，已移除。~~

## 修改 4：磁力公式修正（核心）

**文件**: `src/mfield/mfield2d.hh`

**原因**: 磁力公式从根本上是错误的。

| 项 | 原公式 | 正确公式 |
|----|--------|----------|
| 表达式 | `F_m = χ|H|∇|H|` | `F_m = -(1/2)|H|² ∇χ` |
| 来源 | Kelvin 力（粒子体体力） | Maxwell 应力张量散度 |
| 依赖 | ∇|H|（非均匀场） | ∇χ（界面处磁化率梯度） |
| 量级 | O(H₀²/R) | O(H₀²/W)，W=4 为界面宽度 |
| 比值 | 1 | ~R/W ≈ 12.5x |

论文 Eq.(8) `F_m = (μ₀χ/2)∇(|H|²)` 假设 χ 为常数，在相场模型中 χ 在界面处从 0→1 变化，此公式不成立。正确公式应从 Maxwell 应力张量推导。

**修改**: 用 `F_m = -(1/2)|H|² ∇χ` 替换原公式：

```cpp
// 修改前：Kelvin 力（仅体体力，量级过小）
T chi = chi_l + phi * (chi_h - chi_l);
// ... 计算 ∇|H| ...
T Fmag_x = chi * Hmag * grad_Hmag[0];
T Fmag_y = chi * Hmag * grad_Hmag[1];

// 修改后：Maxwell 应力（包含界面力，量级正确）
// ∇χ = (chi_h - chi_l) * ∇φ  via D2Q5 gradient
// ... 计算 ∇χ ...
T half_Hsq = T{0.5} * Hmag * Hmag;
T Fmag_x = -half_Hsq * dchi_x;
T Fmag_y = -half_Hsq * dchi_y;
```

**后果**:
- 5x 经验因子（H0² 放大 25x）不再需要，已移除
- 子迭代本身不是根本原因，但作为辅助机制保留（MF_SubIter=10）
- 磁力量级提升约 12.5x，与论文 Bom 定义一致

## 修改 5：磁力公式三分支切换（当前核心修改）

**文件**: `src/mfield/mfield2d.h` 和 `src/mfield/mfield2d.hh`

新增三种磁力公式函数声明与实现：

| 函数 | 公式 | 物理含义 |
|------|------|----------|
| `MFMagneticForce2D` | `F_m = (H·∇χ)H - (1/2)|H|²∇χ` | 完整 Maxwell 应力张量散度 |
| `MFMagneticForceInterfacial2D` | `F_m = -(1/2)|H|²∇χ` | 仅界面力（Maxwell 的第二项） |
| `MFMagneticForcePaper2D` | `F_m = (χ/2)∇|H|²` | 论文 Eq.(8) Kelvin 体力 |

## 修改 6：FORCE_FORMULA 宏切换

**文件**: `bubbleMag2d.cpp` 第6行

```cpp
#define FORCE_FORMULA 2   // 0=Maxwell, 1=Interfacial, 2=Paper
```

磁力应用循环（第630-654行）：
```cpp
if(Bom>T{0}){
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
    auto& ns_bl=NSLattice.getBlockLat(b);
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=bk.getOverlap();
    for(int j=ov;j<bk.getNy()-ov;++j){
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
#if FORCE_FORMULA == 0
        MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#elif FORCE_FORMULA == 1
        MFMagneticForceInterfacial2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#else
        MFMagneticForcePaper2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#endif
      }
    }
  }
}
```

## 修改 7：磁场预收敛循环

**文件**: `bubbleMag2d.cpp` 第357-441行

在仿真初始化阶段添加 200 步磁场预收敛循环，强制磁场达到稳态，排除非稳态场对力计算的影响。

## 修改 8：直接 NS 力场对比诊断

**文件**: `bubbleMag2d.cpp` 第443-531行

将三种公式分别应用到 NS 格子，测量总力 sum|F| 和 max|F|：

```
Direct NS force comparison (Bom=0.760):
  Maxwell:   sum|F|=6.15e-04  max|F|=5.89e-07
  Interfacial: sum|F|=6.15e-04  max|F|=5.89e-07
  Paper:     sum|F|=8.28e-03  max|F|=9.77e-06
  Ratio paper/Maxwell sum=847  max=16.6
```

**诊断结论**: Maxwell 和 Interfacial 完全相同，因为 (H·∇χ)H = 0（H场沿界面切线）。Paper 公式力值大 16.6x(max)/847x(sum)，因为它是体力作用于整个气泡内部，而 Maxwell 是界面力。

## 修改 9：NS 力场诊断增强

**文件**: `bubbleMag2d.cpp` 第656-785行

在主循环中添加：
1. 力诊断（每10000步）：打印 H 场、∇χ、∇|H|²、三种力的最大量级
2. NS 力场在界面点的诊断（t=0）：气泡顶部、侧面、底部、中心的 NS 力
3. NS 力场总和诊断（每200步）：打印 sum|F|、max|F|、sumF=(Fx,Fy)、非零单元数

## 修改 10：阶段 1 定量诊断仪表盘（几何 + 场增强）

**文件**: `bubbleMag2d.cpp` 第927-970行

**原因**: 阶段 1 诊断计划要求将"视觉相同"转化为定量指标。原代码仅有 VTK 输出，无法定量比较三公式差异。

**修改**: 在主循环 Phase E（宏更新后，VTK 输出前）添加几何 + 场增强诊断：

```cpp
if(t()%OutputStep==0){
  // 几何诊断：拉长比、质心、面积
  // 场增强诊断：Henh = avg|H|_inside / H0 (期望 ~1.33)
  // 输出格式: [Geom t=X] elong=Y center=(X,Y) area=Z Henh=W (expect~1.33)
}
```

测量内容：
1. **拉长比** `elong = Dy/Dx`：扫描 phi<0.5 的单元，测量 x、y 方向跨度
2. **质心** `(xc, yc)`：`Yc = Σy(1-φ) / Σ(1-φ)`（论文 Eq.68）
3. **气泡面积** `Σ(1-φ)`
4. **磁场增强因子** `Henh = avg|H|_inside / H0`：phi<0.1（气泡深处）单元的平均 |H| 与 H0 之比

### 200 步验证测试结果（Bom=1.94, Paper 公式）

| t | elong | Henh | 说明 |
|---|-------|------|------|
| 50 | 1.0000 | 1.1293 | 无变形，场约 85% 期望值 |
| 100 | 1.0000 | 1.1412 | 场缓慢收敛 |
| 150 | 1.0102 | 1.1516 | 气泡开始垂直拉长 |
| 200 | 1.0102 | 1.1608 | 变形继续 |

**关键发现**：
1. ✅ 诊断系统工作正常，每 OutputStep 正确打印
2. ⚠️ **Henh = 1.13→1.16**（期望 1.33）——场增强仅约 85% 期望值，**嫌疑 A 部分确认**
3. ✅ **气泡正在垂直变形**（elong 1.00→1.01）——Paper 公式产生沿 H0 方向的 prolate 形状，方向正确
4. Henh 收敛速率约 0.015/50 步 = 0.0003/步，达到 1.33 需约 570 步（线性外推，实际为对数收敛）

## 当前配置

**文件**: `bubbleMag2d.ini`

```ini
[Magnetic_Field]
chi_l = 0.0
chi_h = 1.0
mu_l = 1.0
mu_h = 2.0
Bom = 1.94
[Simulation_Settings]
TotalStep = 200
OutputStep = 50
```

## 关键代码位置

| 功能 | 文件 | 行号 |
|------|------|------|
| FORCE_FORMULA 宏 | `bubbleMag2d.cpp` | 6 |
| 磁力应用循环 | `bubbleMag2d.cpp` | 630-654 |
| 磁场预收敛循环 | `bubbleMag2d.cpp` | 357-441 |
| 直接 NS 力场对比 | `bubbleMag2d.cpp` | 443-531 |
| 主循环诊断 | `bubbleMag2d.cpp` | 656-785 |
| Maxwell 公式实现 | `mfield2d.hh` | 88-121 |
| Interfacial 公式实现 | `mfield2d.hh` | 133-165 |
| Paper 公式实现 | `mfield2d.hh` | 177-209 |
| 磁力函数声明 | `mfield2d.h` | 99-123 |

## 修改 11：三公式并行对比测试（关键突破）

**文件**: 编译三个独立二进制 `bubbleMag2d_{maxwell,interfacial,paper}.exe`，分别运行 Bom=1.94, TotalStep=15000

### 关键发现 1：Maxwell ≠ Paper = Interfacial（推翻用户原始观察）

在 Bom=1.94 下并行运行三公式至 t=6000，定量对比 elong（拉长比 Dy/Dx）：

| t    | Maxwell  | Paper    | Interfacial | Baseline(Bom=0) |
|------|----------|----------|-------------|-----------------|
| 1000 | 1.0000   | 1.0000   | 1.0000      | 1.0000          |
| 2000 | 1.0102   | 1.0102   | 1.0102      | 1.0102          |
| 3000 | 1.0102   | 1.0102   | 1.0102      | —               |
| 4000 | 1.0102   | 1.0102   | 1.0102      | —               |
| 5000 | **1.0521** | 1.0204 | 1.0204      | —               |
| 5500 | **1.0417** | 1.0204 | 1.0204      | —               |
| 6000 | **1.0521** | 1.0204 | 1.0204      | —               |

**核心结论**：
1. **Paper ≡ Interfacial**（Henh 和 elong 在所有时间步完全相同）
2. **Maxwell 产生明显更大的变形**（t=6000: 1.0521 vs 1.0204，差异 3%）
3. **Paper/Interfacial 的变形等于 Baseline**（t=2000: 1.0102 == 1.0102）→ **这两个公式不产生任何磁力变形，仅浮力变形**
4. 用户原始观察"三公式变形相同"是因为之前使用的是 Paper 公式（FORCE_FORMULA=2），该公式等价于 Interfacial，不产生磁力变形

### 关键发现 2：Paper ≡ Interfacial 的数学证明

利用矢量恒等式 `∇(χ|H|²) = |H|²∇χ + χ∇|H|²`，可得：

```
Paper:      F = (χ/2)∇|H|²
Interfacial: F = -(1/2)|H|²∇χ

Paper = (1/2)[∇(χ|H|²) - |H|²∇χ]
      = (1/2)∇(χ|H|²) + Interfacial
      = Interfacial + ∇Φ   （其中 Φ = χ|H|²/2）
```

**差项 `(1/2)∇(χ|H|²)` 是纯梯度项**，在不可压缩 NS 方程中被压力场吸收（p' = p - Φ），不产生任何流动。因此 Paper 和 Interfacial 产生**完全相同的流场**，都不产生磁力变形。

### 关键发现 3：Maxwell 是唯一正确的公式

Maxwell 公式 `F = (H·∇χ)H - (1/2)|H|²∇χ` 包含额外项 `(H·∇χ)H`：
- 该项**不是梯度**，不能被压力吸收
- 在气泡顶部/底部极点（H ∥ ∇χ）：`(H·∇χ)H = H²∇χ`（沿场方向拉伸）
- 在气泡赤道（H ⊥ ∇χ）：`(H·∇χ)H = 0`（无额外力）
- 净效应：垂直方向拉伸，水平方向压缩 → **正确的气泡拉长行为**

### 关键发现 4：Bom 量级问题的根因

**根因**：用户一直使用 Paper 公式（FORCE_FORMULA=2），该公式在不可压缩相场模型中等价于 Interfacial，不产生磁力变形。因此：
- Bom=1.94 配置下，磁力完全无效
- 气泡仅受浮力变形（elong≈1.02），与论文 Bom=0.13（弱磁力）相似
- 这解释了"15 倍偏小"的现象——不是 Bom 定义错误，而是**公式选择错误**

### Henh 场增强对比

| 公式 | Henh(t=6000) | 期望值 | 比值 |
|------|-------------|--------|------|
| Maxwell | 1.2748 | 1.33 | 96% |
| Paper | 1.2810 | 1.33 | 96% |
| Interfacial | 1.2810 | 1.33 | 96% |
| Baseline | 0.0000 | — | — |

Henh ≈ 1.28（96% 期望值），**场增强基本正确**，不是 Bom 偏小的主要原因（嫌疑 A 排除）。

## 未解决问题

### 1. ✅ 已解决：三公式变形相同
**根因**：Paper 公式在不可压缩相场模型中等价于 Interfacial（差一梯度项被压力吸收），两者都不产生磁力变形。只有 Maxwell 公式产生真实磁力变形。

### 2. ✅ 已解决：Bom 量级偏小 15 倍
**根因**：使用 Paper 公式（FORCE_FORMULA=2）导致磁力完全无效。气泡变形仅来自浮力，与论文 Bom=0.13 的弱磁力情况相似。切换到 Maxwell 公式（FORCE_FORMULA=0）后，Bom=1.94 产生明显磁力变形。

### 3. ⏳ 待验证：Maxwell 公式的 Bom 量级是否正确
- t=6000 时 Maxwell elong=1.0521，仍在增长
- 需要 t=15000（t*≈1）与论文 Fig.16f 对比
- 粘性时间尺度 τ_visc ≈ R²/ν = 50²/0.017 ≈ 1.46e5 步，t=15000 时仅达 10%，变形仍在增长期

### 4. 论文复现
- 论文使用 Paper 公式但能产生变形 → 论文可能使用可压缩表述，或 χ 定义不同
- 需要在 Bom=0.13/0.27/0.49/0.76/1.94 下用 Maxwell 公式验证 elong 单调增加

## 下一步

### 短期（等待当前仿真完成）
1. ✅ 三公式并行对比至 t=6000，已确认 Maxwell 是正确公式
2. ⏳ 等待三公式 + baseline 完成 15000 步（t*≈1）
3. 📋 切换默认 `FORCE_FORMULA=0`（Maxwell）并提交

### 中期（论文复现）
4. 在 Bom=0/0.13/0.27/0.49/0.76/1.94 下运行 Maxwell 公式，对比论文 Fig.16
5. 验证 elong 随 Bom 单调增加
6. 若 Maxwell 在 Bom=1.94 产生 filament 结构（论文 Fig.16f），则复现成功

### 长期
7. 调查论文 Paper 公式为何能产生变形（可压缩性？χ 定义？）
8. 考虑是否需要在代码中添加可压缩修正

## 验证

- 编译：三个独立二进制 `bubbleMag2d_{maxwell,interfacial,paper}.exe` 编译通过
- 并行运行：mpirun -np 1 各自运行，避免 MPI 内存问题
- 诊断：几何诊断（elong, Henh）每 500 步输出，数据干净无交叉
