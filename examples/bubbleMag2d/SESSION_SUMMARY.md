# Session Summary — bubbleMag2d 磁力公式验证与调试

## 客观目标

**理解并验证三种磁力公式在仿真中产生相同结果的原因，并使仿真结果与论文 (Guo et al., 2025, Phys. Fluids 37, 022148) 一致。**

## 当前代码状态

### 编译通过
```bash
cd /home/dlf/myCode/FerroTest/examples/bubbleMag2d && make -j4
```

### 运行
```bash
mpirun -np 1 ./bubbleMag2d.exe
```

### 当前配置文件 (`bubbleMag2d.ini`)
- `TotalStep = 5000`, `OutputStep = 1000`
- `Bom = 0.76`, `chi_l = 0.0`, `chi_h = 1.0`, `mu_l = 1.0`, `mu_h = 2.0`
- 无子迭代

### 当前 FORCE_FORMULA
在 `bubbleMag2d.cpp` 第6行: `#define FORCE_FORMULA 2` (Paper formula)

## 修改的文件

### 1. `src/mfield/mfield2d.h`
**声明部分** (第99-123行):
- `MFMagneticForce2D` — 完整 Maxwell 公式: `F_m = (H·∇χ)H - (1/2)|H|²∇χ`
- `MFMagneticForceInterfacial2D` — 仅界面力: `F_m = -(1/2)|H|²∇χ`
- `MFMagneticForcePaper2D` — 论文 Eq.(8) 体力: `F_m = (χ/2)∇|H|²`

### 2. `src/mfield/mfield2d.hh`
**实现部分** (第69-209行):
- `MFMagneticForce2D::apply()` — 完整 Maxwell 公式实现
- `MFMagneticForceInterfacial2D::apply()` — 仅界面力实现
- `MFMagneticForcePaper2D::apply()` — 论文公式实现

### 3. `examples/bubbleMag2d/bubbleMag2d.cpp`
**关键修改**:
- 第6行: `#define FORCE_FORMULA` 宏，0=Maxwell, 1=Interfacial, 2=Paper
- 第630-654行: 磁力应用循环，使用 `FORCE_FORMULA` 选择公式
- 第370-441行: 磁场预收敛后诊断 (打印三种公式的力值)
- 第443-531行: 直接 NS 力场对比 (将三种公式分别应用到 NS 格子并测量总力)
- 第656-730行: 主循环中每步的诊断输出
- 移除: 子迭代循环、5x 经验补偿因子

### 4. `bubbleMag2d.ini`
- 移除 `MF_SubIter` 参数
- `TotalStep` 从 20000 → 200 → 100 → 5000

## 关键诊断发现

### 发现 1: Maxwell 和 Interfacial 公式完全相同
```
NS Force comparison:
  Maxwell:     sum|F|=6.15e-04  max|F|=5.89e-07
  Interfacial: sum|F|=6.15e-04  max|F|=5.89e-07  ← 完全相同
  Paper:       sum|F|=8.28e-03  max|F|=9.77e-06  ← 16.6x/847x 更大
```

**原因**: `(H·∇χ)H = 0` 处处成立，因为 H 场与界面法线垂直（2D 柱体在横向磁场中，H 场沿界面切线方向）。

### 发现 2: Maxwell 公式在界面处的物理行为
- 顶部极点 (n ∥ H): `F_m = (1/2)H₀²∇χ` (向外，拉伸)
- 侧面赤道 (n ⟂ H): `F_m = -(1/2)H₀²∇χ` (向内，压缩)
- 净效应: 气泡沿磁场方向拉长

### 发现 3: 论文公式 vs Maxwell 公式
- 论文公式 `F_m = (χ/2)∇|H|²` 是 Kelvin 体力，作用于整个气泡内部
- Maxwell 公式 `F_m = (H·∇χ)H - (1/2)|H|²∇χ` 是界面力，仅作用于界面
- 当 χ 在界面处连续变化时，两者不等价（论文公式假设 χ 为常数）
- 论文公式的力值比 Maxwell 公式大 16.6x (max) / 847x (sum)

## 未解决的问题

### 1. 三种公式结果相同的原因
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

## 下一步建议

### 短期 (诊断验证)
1. **运行 `FORCE_FORMULA=0` (Maxwell) 和 `FORCE_FORMULA=2` (Paper) 各跑一次 `TotalStep=15000`**，对比 VTK 输出中 PHI 场的差异
2. 在 `Bom=0` (无磁场) 下运行同样参数，作为基准对照
3. 比较三个案例的 `t*=1` 时气泡形状

### 中期 (问题定位)
4. 如果三种公式仍无差异，在 NS 碰撞步骤后添加力场诊断，确认力是否被 NS 求解器正确使用
5. 检查 MRTForce 碰撞是否正确处理 FORCE 场

### 长期 (论文复现)
6. 确认使用哪种公式（Maxwell 或 Paper）更符合论文物理模型
7. 调整参数使 t*=1 时的气泡形状与论文 Fig. 4(c) 对齐
8. 测试 `Bom=1.94` 和 `Bom=0.13` 等不同磁场强度

## 文件结构
```
/home/dlf/myCode/FerroTest/
├── src/mfield/
│   ├── mfield2d.h      # 磁力函数声明 (第99-123行)
│   └── mfield2d.hh     # 磁力函数实现 (第69-209行)
├── examples/bubbleMag2d/
│   ├── bubbleMag2d.cpp  # 主仿真代码 (FORCE_FORMULA 第6行)
│   ├── bubbleMag2d.ini  # 配置参数
│   ├── CHANGELOG.md     # 修改历史
│   └── SESSION_SUMMARY.md  # 本文件
└── paper/
    └── Phase-field lattice Boltzmann model with adaptive mesh refinement for ferrofluid interfacial dynamics.md
```

## 关键代码位置

| 功能 | 文件 | 行号 |
|------|------|------|
| FORCE_FORMULA 宏 | `bubbleMag2d.cpp` | 6 |
| 磁力应用循环 | `bubbleMag2d.cpp` | 630-654 |
| 磁场预收敛循环 | `bubbleMag2d.cpp` | 357-441 |
| 直接 NS 力场对比 | `bubbleMag2d.cpp` | 443-531 |
| 主循环诊断 | `bubbleMag2d.cpp` | 656-730 |
| Maxwell 公式实现 | `mfield2d.hh` | 88-121 |
| Interfacial 公式实现 | `mfield2d.hh` | 133-165 |
| Paper 公式实现 | `mfield2d.hh` | 177-209 |
| 磁力函数声明 | `mfield2d.h` | 99-123 |