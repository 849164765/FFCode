# Bubble2D — 全代码库重构变更记录

## 项目目标

基于 FreeLB 代码库，重新设计并实现正确的 FDLBM 铁磁流体求解器。从 bubble2d 算例开始逐步推进。

---

## Phase 1: Converter 系统重构 **[已完成]**

### 1.1 PhaseFieldConverter 改造 (`src/lbm/unit_converter.h`)
- 新增成员: `cs2`, `Lattice_RT_phi`, `OMEGA_phi`, `Lattice_beta`, `Lattice_kappa`
- 构造函数增加 `cs2` 参数: `PhaseFieldConverter(BaseConverter<T>&, T cs2)`
- 新增 `fromLattice(W_lat, M_lat, sigma_lat, tau_phi)` — 格子单位输入（绕过物理→格子转换）
- 原 `Converter(W_phys, M_phys, sigma_phys)` 增强：自动计算 lattice beta/kappa/tau/omega
- 实现完整 `AbstractConverter` 虚方法（不再返回 `T(0)`）
- 新增 getter: `getLatticeBeta()`, `getLatticeKappa()`

### 1.2 BaseConverter Bug 修复
- `Lattice_Re = Lattice_charU / Lattice_VisKine` → `Lattice_Re = Lattice_charU * Lattice_charL / Lattice_VisKine`
- 硬编码: `Lattice_g = T(9810) / Conv_Acc`（bubble2d 不受影响）

### 1.3 删除未使用转换器
- 删除: `TempConverter`, `ConcConverter`, `PhaseDiagramConverter`, `ZSConverter`, `GandinConverter`
- 简化 `UnitConvManager`：仅 BaseConv + PhaseFieldConv
- 新增 `PhaseFieldConverter::check()` 实现

### 1.4 示例文件重构（12个文件）
废除所有手动物理→格子转换，统一使用 Converter。批改文件：
`bubble2d`, `bubble2dMRT`, `bubblecoalescence2d`, `simpledrop2d/2dMRT/2dV2/3d/3dMRT/3dV2`, `rti2d/rti2dMRT`, `shearflow2dMRT`

### 1.5 无量纲数守恒验证
Re, Eo, We, Cn, βW/σ, κ/(σW) 全部验证通过。

---

## Phase 2: 主循环算法修复 **[已完成]**

参考标准代码 `sheardrop2(1).f90`（MRT LBM 多相流），重排主循环为标准顺序：

```
Phase 1 (力):    Clear FORCE → ST_force → Gravity_force
Phase 2 (碰撞):  PF_BGKSource → NS_BGKForce(WriteToField=false)
Phase 3 (迁移):  PF_BC+Stream+Comm → NS_BC+Stream+Comm
Phase 4 (宏观):  NS_macro(ρ=Σf, u=(Σcf+0.5F)/ρ) → PF_φ←Σg → PF_∇φ,λ
```

关键改动：
- 碰撞/宏观分离：NS碰撞用 `forcerhoU<false>`（纯BGK公式），宏观更新用独立 `forcerhoU<true>`（半力修正）
- 主循环前增加 PF 宏观量初始化（∇φ, λ）

---

## Phase 3: 参数系统优化 **[已完成]**

### 3.1 双模式输入
- 无量纲设计模式（默认）：Re, Eo, U_g → 自动推导 ν, g, σ
- 直接指定模式（`use_direct=1`）：nu_l, nu_h, gravity, sigma

### 3.2 新增 INI 参数
- `viscosity_ratio` — 动力粘度比 η_h/η_l（替代硬编码 /10）
- `nu_l`, `nu_h` — 直接指定模式下的轻/重流体运动粘度

### 3.3 TC2 基准算例
创建 `examples/bubble2d_direct/`，配置如下：
```ini
Ni=256, Nj=512, Radius=64 (D=128), W=6
rho_l=1, rho_h=1000
Re=35, Eo=125, U_g=0.1, viscosity_ratio=100
```
推导参数：ν=0.366, g=-7.81e-5, σ=10.23, η_h=366, η_l=3.66

---

## 当前状态

### 关键文件
| 路径 | 说明 |
|------|------|
| `examples/bubble2d/bubble2d.cpp` | 主算例代码（Phase 1-3 完整） |
| `examples/bubble2d/bubble2d.ini` | 原始测试参数 (rho_h=2, Eo=125) |
| `examples/bubble2d_direct/bubble2d.cpp` | TC2 副本（与上一致） |
| `examples/bubble2d_direct/bubble2d.ini` | TC2 基准参数 (rho_h=1000, Re=35, Eo=125) |
| `src/lbm/unit_converter.h` | BaseConverter + PhaseFieldConverter + UnitConvManager + RefineConverter |
| `src/lbm/unit_converter.hh` | check() 实现 |
| `src/ff/ff2d.h/hh` | PF field 类型、BroadcastAllParams、functor 实现 |

### 已知问题
- TC2 参数下 BGK 碰撞稳定性待验证（tau_ns≈1.6, 密度比1000:1）
- 推荐使用 Interface_Width=6（原版=4）以提升高Eo下稳定性
- 恢复原设定：rho_h 回 2，Eo=125，viscosity_ratio=10

---

## 待处理

- [ ] TC2 基准复刻运行验证
- [ ] 铁磁流体模块开发 (ff/ 目录扩展)
- [ ] 磁场相关转换器设计
- [ ] shearflow 算例基于参考代码重写
- [ ] MRT 碰撞算子适应高密度比场景
