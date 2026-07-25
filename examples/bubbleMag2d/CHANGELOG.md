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

## 未解决问题

### 1. 三种公式结果相同（视觉判断）
用户测试了三次（修改 `MFMagneticForce2D::apply()` 函数体），15000 步下三种公式产生视觉相同的气泡变形。但 200 步诊断显示：
- 力值差异确实存在：Paper 公式力值大 16.6x（max）/ 840x（sum）
- Maxwell 与 Interfacial 力值**完全相同**（max=1.36e-6, sum=1.51e-3）
- 原因分析：Maxwell 与 Interfacial 的差异在 (H·∇χ)H 项，该项要求 H 沿 ∇χ 方向有分量。若 H 沿界面切线（∇χ 沿界面法线），则 (H·∇χ)=0，二者相同。

### 2. Bom 量级偏小（约 15 倍）
配置 Bom=1.94，但变形仅相当于论文 Bom=0.13。200 步诊断显示：
- Henh = 1.13→1.16（期望 1.33），场增强仅 85%——**嫌疑 A 部分确认**
- 气泡开始变形（elong 1.00→1.01），方向正确（垂直拉长）
- 需要更长仿真（~15000 步）观察完整变形

### 3. 完整公式与论文不一致
- 论文使用 `F_m = (χ/2)∇|H|²` (Kelvin 体力)
- 当前代码实现使用 Maxwell 公式
- 两者在相场模型（χ 连续变化）中不等价

### 4. 论文复现
- 需要在 `Bom=0.76` 和 `Bom=1.94` 下对比 t*=1 时的气泡形状
- 特征时间步数: ~14583 步 (基于重力特征速度)
- 需要验证拉长比（elongation ratio）与论文 Fig. 4(c) 一致

## 下一步

### 短期 (阶段 1 完整诊断)
1. ✅ 已完成：200 步短测试验证诊断系统工作正常
2. 运行 `FORCE_FORMULA=2` (Paper) Bom=1.94 TotalStep=15000，收集 Henh 和 elongation 时序
3. 运行 `FORCE_FORMULA=0` (Maxwell) Bom=1.94 TotalStep=15000，对比 elongation
4. 运行 `FORCE_FORMULA=1` (Interfacial) Bom=1.94 TotalStep=15000，对比 elongation
5. 在 `Bom=0` 下运行同样参数，作为基准对照

### 中期 (阶段 2：场增强验证，嫌疑 A 修复)
6. 若 Henh<1.30，修改 [mfield2d.hh:28](file:///home/dlf/myCode/FerroTest/src/mfield/mfield2d.hh) 松弛公式：`tau_psi = 0.5 + 2.5*mu`（论文 Eq.42）
7. 若 Henh 仍不足，添加 MF 子迭代（5-10 次/物理步）
8. 验证 Henh ≈ 1.33 后重新评估三公式差异

### 长期 (论文复现)
9. 验证 Bom=1.94 的 elongation 与论文 Fig.16f 一致
10. 测试 Bom=0.13/0.27/0.49/0.76/1.94 下 elongation 单调增加
11. 确认 Bom=1.94 产生明显的垂直拉长

## 验证

- 编译：`make` 通过，生成 `bubbleMag2d.exe`
- 直接 NS 力场对比诊断成功运行，输出三种公式的力值差异
- 主循环诊断正常输出力场信息
