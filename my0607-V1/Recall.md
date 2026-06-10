# Code Change Recall

## 2026-06-09 — Phase 1: Converter系统重构

### 修改文件

| 文件 | 位置 | 操作 | 摘要 |
|------|------|------|------|
| `src/lbm/unit_converter.h` | :497-538 | 修改 | PhaseFieldConverter: 新增 cs2/Lattice_RT_phi/OMEGA_phi/Lattice_beta/Lattice_kappa 成员，构造函数增加 cs2 参数，新增 fromLattice() 方法，实现完整 AbstractConverter 虚方法 |
| `src/lbm/unit_converter.h` | :112 | 修改 | Bug修复: `Lattice_Re = Lattice_charU / Lattice_VisKine` → `Lattice_Re = Lattice_charU * Lattice_charL / Lattice_VisKine` |
| `src/lbm/unit_converter.h` | :142-493 | 删除 | 删除 TempConverter, ConcConverter, PhaseDiagramConverter, ZSConverter, GandinConverter 五个类 |
| `src/lbm/unit_converter.h` | :583-625 | 修改 | UnitConvManager 简化为 BaseConv + PhaseFieldConv，构造函数从3个减为2个 |
| `src/lbm/unit_converter.h` | :163 | 修改 | 初始化顺序修复: `cs2(cs2_), BaseConv(baseconv)` 匹配声明顺序 |
| `src/lbm/unit_converter.hh` | :45-151 | 删除 | 删除 TempConverter/ConcConverter/ZSConverter/GandinConverter::check() |
| `src/lbm/unit_converter.hh` | :153-175 | 新增 | PhaseFieldConverter::check() 实现 |
| `src/lbm/unit_converter.hh` | :177-220 | 修改 | UnitConvManager::Check_and_Print() 简化，增加 PhaseFieldConv 输出 |
| `examples/bubble2d/bubble2d.cpp` | :167-174 | 修改 | Converter: SimplifiedConverterFromRT → ConvertFromRT(1, Tau_ns, rho_h, D, U_g, nu_lat); PhaseConv.fromLattice() |
| `examples/bubble2d/bubble2d.cpp` | :224 | 修改 | PFLattice 用 PhaseConv 替代 PFBaseConv |
| `examples/bubble2d/bubble2d.cpp` | :231-235 | 修改 | BroadcastAllParams 用 PhaseConv.getLatticeBeta/Kappa() |
| `examples/simpledrop2d/simpledrop2d.cpp` | :34-38, :66-72 | 删除 | 删除 Kappa/Beta/Tau_phi/Omega_phi 全局变量和手工计算 |
| `examples/simpledrop2d/simpledrop2d.cpp` | :125-141 | 修改 | dummyConv → BaseConv+PhaseConv.fromLattice(), BlockLatticeManager 用 PhaseConv |
| `examples/simpledrop2dMRT/simpledrop2d.cpp` | 同上模式 | 修改 | 与 simpledrop2d 相同改动 |
| `examples/simpledrop3d/simpledrop3d.cpp` | 同上模式 | 修改 | 3D 版本，相同改动 |
| `examples/simpledrop3dMRT/simpledrop3d.cpp` | 同上模式 | 修改 | 同上 |
| `examples/simpledrop2dV2/simpledrop2d.cpp` | 同上模式 | 修改 | 同上 |
| `examples/simpledrop3dV2/simpledrop3d.cpp` | 同上模式 | 修改 | 同上 |
| `examples/rti2d/rti2d.cpp` | Tau_ns/Kappa/Beta/Tau_phi | 删除 | 删除手动计算，NSConv→BaseConv, PFConv→PhaseConv |
| `examples/rti2dMRT/rti2d.cpp` | 同上模式 | 修改 | 与 rti2d 相同改动 |
| `examples/shearflow2dMRT/shearflow2d.cpp` | Tau_ns/Tau_phi/Omega/Kappa/Beta | 删除 | 删除6个手动计算变量，NS/PF用新Converter |
| `examples/bubble2dMRT/bubble2d.cpp` | 同bubble2d模式 | 修改 | ConvertFromTimeStep + PhaseConv.fromLattice() |
| `examples/bubblecoalescence2d/bubblecoalescence2d.cpp` | 同bubble2d模式 | 修改 | 同上 |

**原因**: 废除所有手动物理→格子转换，统一由 Converter 负责；无量纲数守恒验证全部通过 (Re/Eo/We/Cn)

---

## 2026-06-09 — Phase 2: 主循环算法修复

### 修改文件

| 文件 | 位置 | 操作 | 摘要 |
|------|------|------|------|
| `examples/bubble2d/bubble2d.cpp` | :329-336 | 修改 | NS任务: BGKForce 用 `forcerhoU<false>` (WriteToField=false, 纯BGK) |
| `examples/bubble2d/bubble2d.cpp` | :339-341 | 新增 | NSMacroTask: 独立 `forcerhoU<true>` 半力修正任务 |
| `examples/bubble2d/bubble2d.cpp` | :394-402 | 新增 | 主循环前 PF 宏观量初始化 (∇φ, λ) |
| `examples/bubble2d/bubble2d.cpp` | :410-488 | 修改 | 主循环重排: Phase1(力) → Phase2(碰撞:PF→NS) → Phase3(迁移:PF→NS) → Phase4(宏观:NS_macro→φ→∇φ,λ) |

**原因**: 参考 sheardrop2(1).f90 标准算法，修正为力→碰撞→迁移→宏观顺序；碰撞/宏观分离

---

## 2026-06-09 — Phase 3: 参数系统与TC2基准

### 新增文件

| 文件 | 行号 | 摘要 |
|------|------|------|
| `examples/bubble2d_direct/bubble2d.ini` | L1-L28 | TC2 基准参数: Radius=64, rho_h=1000, Re=35, Eo=125, viscosity_ratio=100, W=6 |
| `.claude/rules/changelog.md` | L1-L7 | 项目规则: 每次对话自动读取 CHANGELOG.md |
| `examples/bubble2d/CHANGELOG.md` | L1-L102 | 三阶段完整重构记录 |

### 修改文件

| 文件 | 位置 | 操作 | 摘要 |
|------|------|------|------|
| `examples/bubble2d/bubble2d.cpp` | :66 | 新增 | 读入 viscosity_ratio 参数 |
| `examples/bubble2d/bubble2d.cpp` | :81-82 | 修改 | `eta_l = eta_h / T(10)` → `eta_l = eta_h / viscosity_ratio` |
| `examples/bubble2d/bubble2d.ini` | :22 | 新增 | `viscosity_ratio = 10` |
| `examples/bubble2d/bubble2d.cpp` | :169-179 | 修改 | MPI_RANK → IF_MPI_RANK (main返回int时return无值bug) |
| `examples/bubble2d/bubble2d.cpp` | :234-237 | 修改 | BroadcastAllParams 用局部变量 Beta/Kappa (getLatticeBeta()返回临时量不能绑T&) |
| `examples/bubble2d/bubble2d.cpp` | :248 | 删除 | 移除未使用的 blockLat 变量 |

**原因**: TC2 粘度比 100:1 (μ₁/μ₂=10/0.1) vs 原硬编码 10:1; 编译bug修复

---

## 2026-06-10 — FF2D 重构：CA 模式设计

### 新增文件

| 文件 | 行号 | 摘要 |
|------|------|------|
| `src/ff/ff2d.h` | L1-L76 | FF2D 模块声明: GRADBase/NORMALBase/LAPLACIANBase 场基, GRAD/NORMAL/LAPLACIAN 类型别名, FFFIELDS/REFLBMFIELDS/FIELDPACK/ALLFIELDS 场包, FF2DBlock 单块类, FF2DManager 多块管理器 |
| `src/ff/ff2d.hh` | L1-L143 | FF2DBlock 实现 (构造/computeGradient/computeLaplacian/apply), FF2DManager 实现 (构造/Init/apply/Communicate) |

### FF2D 架构 (仿 CA::BlockZhuStefanescu2D 模式)

| CA 组件 | FF2D 对应 | 说明 |
|---------|----------|------|
| `BlockZhuStefanescu2D` | `FF2DBlock` | 单块算子, 继承 `BlockLatticeBase`, 持有 `Converter&` |
| `BlockZhuStefanescu2DManager` | `FF2DManager` | 多块管理器, 继承 `BlockLatticeManagerBase`, 持有 `vector<FF2DBlock>` |
| `CAFIELDS` (自管) | `FFFIELDS` | GRAD, NORMAL, LAPLACIAN — FF2D 计算的场 |
| `REFFIELDS` (引用) | `REFLBMFIELDS` | PHI (from PFLattice), VELOCITY (from NSLattice) |
| 维度 | 硬编码 D=2 | 2D 专用 |

### 修改文件

| 文件 | 位置 | 操作 | 摘要 |
|------|------|------|------|
| `src/utils/alias.h` | :155-156 | 删除 | GRADBase, NORMALBase 移入 ff2d.h |
| `src/utils/alias.h` | :190-193 | 删除 | 非 CUDA 版 GRAD<T,D>, NORMAL<T,D> 别名 |
| `src/utils/alias.h` | :263-266 | 删除 | CUDA 版 GRAD<T,D>, NORMAL<T,D> 别名 |
| `examples/bubble2dImprove/bubble2d.cpp` | :202-218 | 修改 | PFFIELDS 精简为 TypePack<PHI, POP>; PFFIELDREFS 加入 ff::GRAD/NORMAL/LAPLACIAN; PFLattice 构造传 nullptr 占位 + addField 填入 |
| `examples/bubble2dImprove/bubble2d.cpp` | :268-270 | 新增 | NS 宏场显式初始化: VELOCITY置零, RHO置rho_init |
| `examples/bubble2dImprove/bubble2d.cpp` | :272-286 | 新增 | FF2DMgr 创建(自管GRAD/NORMAL/LAPLACIAN, 引用PHI+VELOCITY) + Init + addField |
| `examples/bubble2dImprove/bubble2d.cpp` | :122-400 | 删除 | 移除所有旧 per-cell 任务 (FF2D<CELL>/FFLaplacian2D/FFChemPotential2D 等), 主循环清空待重写 |

**原因**: 原 FF2D 是 per-cell 模板操作符, 重构为 CA 风格的三层架构 (场定义 → 单块算子 → 多块管理器)。GRAD/NORMAL/LAPLACIAN 从全局 alias.h 移入 ff 命名空间自管, 与 CA::STATE/FS 等一致。PFLattice 通过 addField 引用 FF2D 场 (同 SOLattice 引用 CA::EXCESSC)。
