# Bubble2D 开发进度记录

## 当前状态 (2026-06-14)

### 已完成的修改

#### 1. 界面厚度修正 (最初修改)
- `bubble2d.ini`: `Interface_Width = 4.0 → 11.0`
- `bubble2d.cpp:274`: 初始化公式修正为 `tanh(2√2*(r-R)/W)`，匹配理论平衡解
- 原因: W=4 时界面过渡区仅 ~4 个网格点，对 10:1 变物性太薄。W=11 时界面有效宽度 ~11.5 lu

#### 2. U_g 速度参数设计 (加速模拟)
- `bubble2d.ini`: 新增 `U_g = 0.02` 参数
- `bubble2d.cpp:87-113`: 二选一逻辑 — 若 U_g>0，反推 η、g、σ (保持 Re、Eo 不变)
- 若 U_g=0，走原始 Re→g→U_g 流程 (向后兼容)
- 速度提升约 64x

#### 3. Phi 约定反转：气泡内 φ=1，气泡外 φ=0
- `bubble2d.ini`: `rho_l=10.0, rho_h=1.0, eta_h=0.02, eta_l=0.4`
- `bubble2d.cpp:274`: `phi = 0.5 - 0.5*tanh(2√2*(r-R)/W)`
- `_l` 和 `_h` 仅表示"φ 取该值时的属性"，不表示物理轻重

#### 4. 变粘度接入
- `bubble2d.cpp:221-222`: NS 场包加入 `OMEGA<T>`
- `bubble2d.cpp:390-395`: 注册 `FFRhoOmegaUpdate2D` 耦合任务
- `bubble2d.cpp:476`: 主循环 Step 5.5 调用 RhoOmegaUpdate
- `collision.h:59-65`: `BGKForce` 用 `if constexpr hasField<OMEGA<T>>()` 读 per-cell omega
- 向后兼容：无 OMEGA 字段时自动回退到块级 `getOmega()`

#### 5. 插值公式升级
- **密度**: `g(φ)=φ²(3-2φ)`, `ρ(φ)=ρ_l + g(φ)(ρ_h-ρ_l)` — 平滑多项式
- **粘度**: 调和平均 `1/(τ-0.5) = (1-φ)/(τ_L-0.5) + φ/(τ_G-0.5)` — 保持界面剪切应力连续性

#### 6. 化学势梯度表面张力算子
- `ff2d.h:135-139`: 新增 `FFSurfaceTensionChemPot2D` 声明
- `ff2d.hh:122-153`: 实现 `F_s = -φ·∇λ / ρ(φ)` (与 `F_s = λ·∇φ` 等价，差 ∇(φλ) 被压力吸收)
- `bubble2d.cpp:381`: 两行互换注释即可切换

#### 7. 变密度 LBM 尝试 (已回滚)
- `moment.h:461-510`: 新增 `variableDensityRhoU` 矩方案（通用版）
- `moment.ur.h:554-568`: 新增 D2Q9 展开特化
- 已在 bubble2d.cpp 中回滚到 `forcerhoU`（恢复 ρ≈1.0 标称密度）
- 代码保留未删除，`moment::variableDensityRhoU` 可随时切换回去

#### 8. 破坏性测试 (已清理)
- 在主循环 NS Stream 后、PF 碰撞前将速度场强制清零
- 结果：速度置零 → 气泡不变形 → 证明**变形是流场驱动，非相场 bug**
- 测试代码已删除

#### 9. Boussinesq 参数校准 (关键修复)
- **问题**: LBM solver 中惯性密度 ρ₀=1.0，密度比 10:1 导致有效加速度被放大 10 倍
  - `a_LBM = F/ρ₀ = Δρ·g/1.0 = 9g`，而 `a_phys = Δρ·g/ρ_L = 9g/10 = 0.9g`
  - 终端速度被放大 √10≈3.16 倍，底部动压暴涨 10 倍
- **修复 1** (`bubble2d.cpp:80`): `DeltaRho = abs(rho_h - rho_l)` — 修复反号 bug
- **修复 2** (`bubble2d.cpp:102,115`): `sigma *= rho_l` — σ 从 0.0288 放大到 0.288 (10x)
- **修复 3** (`ff2d.hh:196-197`): `ν = η` 代替 `ν = η/ρ_phys` — ν 基于 ρ₀=1 计算
  - ν_L = 5.132, ν_G = 0.2566 (20:1)，τ_L=15.9, τ_G=1.27
  - 界面调和平均 τ(φ=0.5)≈1.97 → 合理范围

### 修改涉及的文件清单

| 文件 | 改动内容 |
|------|---------|
| `bubble2d.ini` | W=11, U_g=0.02, rho_l=10, rho_h=1, eta_h=0.02, eta_l=0.4 |
| `bubble2d.cpp` | U_g 推导, φ 初始化反转, OMEGA 字段, RhoOmega 任务, Boussinesq σ 校准, 化学势 ST 选项 |
| `src/lbm/collision.h` | BGKForce per-cell OMEGA (if constexpr) |
| `src/lbm/moment.h` | 新增 `variableDensityRhoU` (保留未用) |
| `src/lbm/moment.ur.h` | 新增 `variableDensityRhoU` D2Q9 特化 (保留未用) |
| `src/ff/ff2d.h` | 新增 `FFSurfaceTensionChemPot2D` 声明 |
| `src/ff/ff2d.hh` | `FFSurfaceTension2D` 加 /ρ, `FFSurfaceTensionChemPot2D` 实现, `FFRhoOmegaUpdate2D` 平滑多项式+调和平均+Boussinesq ν |

### 当前模拟参数

```
D = 80, U_g = 0.02, Re = 6.236, Eo = 10
ρ: φ=0(外)液体=10, φ=1(内)气体=1
η: η_l=5.132, η_h=0.2566
ν: ν_L=5.132, ν_G=0.2566 (Boussinesq: ρ₀=1)
g = -5e-6, σ = 0.288
W = 11, tau_phi = 0.8
τ_L = 15.9, τ_G = 1.27, τ_interface ≈ 1.97
CFL = 0.597
```

### 待测试/待定项
- [ ] 运行 Boussinesq 校准后的模拟，观察气泡底部变形是否改善
- [ ] 若仍有内凹，考虑切换到 `FFSurfaceTensionChemPot2D`
- [ ] `variableDensityRhoU` 矩方案已实现但未启用，必要时可切回
- [ ] 其他示例 (bubble2dMRT, bubblecoalescence2d) 需同步修改
