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

## 修改 3：ini 配置文件添加子迭代参数

**文件**: `bubbleMag2d.ini`

```ini
[Magnetic_Field]
chi_l = 0.0
chi_h = 1.0
mu_l = 1.0
mu_h = 2.0
Bom = 0.76
MF_SubIter = 10        ; <-- 新增，每物理步磁场子迭代次数
```

## 修改 4：修正磁力公式（核心！）

**文件**: [`src/mfield/mfield2d.hh`](file:///home/dlf/myCode/FerroTest/src/mfield/mfield2d.hh)

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

## 验证

- 编译：`make` 通过，生成 `bubbleMag2d.exe`
- 启动验证：日志显示 `H0=0.002 Bom=0.760 SubIter=10`，参数正确读取