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

## 当前配置

**文件**: `bubbleMag2d.ini`

```ini
[Magnetic_Field]
chi_l = 0.0
chi_h = 1.0
mu_l = 1.0
mu_h = 2.0
Bom = 0.76
[Simulation_Settings]
TotalStep = 200
OutputStep = 200
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

### 1. 三种公式结果相同
用户测试了三次（修改 `MFMagneticForce2D::apply()` 函数体），结果图像完全相同。但诊断显示 Paper 公式的力值大 16.6x。可能原因：
- 仿真步数不足以产生可见差异
- 磁力被其他力（表面张力、重力）掩盖
- 修改 `apply()` 时实现有 bug

### 2. 完整公式与论文不一致
- 论文使用 `F_m = (χ/2)∇|H|²` (Kelvin 体力)
- 当前代码实现使用 Maxwell 公式
- 两者在相场模型（χ 连续变化）中不等价

### 3. 论文复现
- 需要在 `Bom=0.76` 和 `Bom=1.94` 下对比 t*=1 时的气泡形状
- 特征时间步数: ~14583 步 (基于重力特征速度)
- 需要验证拉长比（elongation ratio）与论文 Fig. 4(c) 一致

## 下一步

### 短期 (诊断验证)
1. 运行 `FORCE_FORMULA=0` (Maxwell) 和 `FORCE_FORMULA=2` (Paper) 各跑一次 `TotalStep=15000`，对比 VTK 输出中 PHI 场的差异
2. 在 `Bom=0` (无磁场) 下运行同样参数，作为基准对照
3. 比较三个案例的 `t*=1` 时气泡形状

### 中期 (问题定位)
4. 如果三种公式仍无差异，在 NS 碰撞步骤后添加力场诊断，确认力是否被 NS 求解器正确使用
5. 检查 MRTForce 碰撞是否正确处理 FORCE 场

### 长期 (论文复现)
6. 确认使用哪种公式（Maxwell 或 Paper）更符合论文物理模型
7. 调整参数使 t*=1 时的气泡形状与论文 Fig. 4(c) 对齐
8. 测试 `Bom=1.94` 和 `Bom=0.13` 等不同磁场强度

## 验证

- 编译：`make` 通过，生成 `bubbleMag2d.exe`
- 直接 NS 力场对比诊断成功运行，输出三种公式的力值差异
- 主循环诊断正常输出力场信息
