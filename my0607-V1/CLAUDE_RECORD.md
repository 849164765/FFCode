# Claude Record

[2026-06-05] [开始] [全局] [] 开始实现相场+流场耦合气泡上升模拟，扩展ff/目录并创建bubble2d/

[2026-06-05] [修改] [src/utils/alias.h] [L158-163] 添加LAPLACIANBase、CHEMICALPOTENTIALBase字段基础类和新类型别名LAPLACIAN<T>、CHEMICALPOTENTIAL<T>

[2026-06-05] [修改] [src/ff/ff2d.h] [全文] 添加ff::params命名空间中的双相流全局参数变量，新增5个结构体声明：FFLaplacian2D、FFChemPotential2D、FFSurfaceTension2D、FFGravityForce2D、FFRhoOmegaUpdate2D

[2026-06-05] [修改] [src/ff/ff2d.hh] [全文] 实现5个新ff结构体：FFLaplacian2D（∇²φ标准各向同性FD）、FFChemPotential2D（化学势λ=4βφ(φ-1)(φ-0.5)-κ∇²φ）、FFSurfaceTension2D（耦合：F_s=λ∇φ→NS FORCE）、FFGravityForce2D（耦合：F_b[1]=ρg）、FFRhoOmegaUpdate2D（耦合：ρ、η插值更新）

[2026-06-05] [创建] [examples/bubble2d/bubble2d.cpp] [全文] 气泡上升主程序：双D2Q9晶格（相场+NS流场），BGKSource+BGKForce碰撞，表面张力+浮力耦合，反弹边界

[2026-06-05] [创建] [examples/bubble2d/Makefile] [全文] 编译配置

[2026-06-05] [创建] [examples/bubble2d/bubble2d.ini] [全文] 仿真参数配置文件(Ni=256,Nj=512,R=50,Eo=125,Re=40等)

[2026-06-05] [验证] [examples/bubble2d/] [] 编译成功，32x32网格3步测试通过，VTK输出正常，TotalPhi守恒监测正常

[2026-06-05] [重构] [examples/bubble2d/bubble2d.cpp] [全文] 参照simpledrop2dV2的MPI并行风格重构主程序：添加mpi().init()、BlockGeometryHelper2D域分解、block迭代PHI更新、NormalCommunicate通信模式、IF_MPI_RANK宏、分离NS/PF两个BaseConverter分别设置弛豫时间

[2026-06-06] [诊断] [examples/bubble2d/bubble2d.cpp] [全文] 分析3个问题：气泡不上升、形状异常、20000步后停止变形。根因：1)重力极小(|g|=1.25e-7，因rho_h=100,D=100使Eo公式分母过大)；2)FFRhoOmegaUpdate2D写入的逐格OMEGA<T>场完全无效（Cell::getOmega()返回BlockLattice标量而非OMEGA场，见cell.h:168）；3)表面张力松弛完成后无浮力驱动。

[2026-06-06] [修复] [examples/bubble2d/bubble2d.cpp] [多处] 修复方案：1)移除无效的OMEGA逐格更新耦合（OMEGA<T>从NS场移除）；2)eta_h通过SimplifiedConverterFromRT自动计算τ_ns来准确控制粘度；3)增添Umax诊断输出；4)改用晶格尺度参数(R=20D=40,rho_h=2,sigma=0.005,Eo=100)→|g|=1.56e-4(1000倍增大)，Re_actual≈40自动匹配

[2026-06-06] [修复] [examples/bubble2d/bubble2d.ini] [全文] 更新为晶格尺度参数：R=20,rho_h=2,sigma=0.005,Eo=100，eta_h=0.158自动匹配Re≈40,τ_ns=0.737

[2026-06-06] [注释] [src/ff/ff2d.hh] [L106-109] 添加警告注释：Cell::getOmega()返回块级标量，逐格OMEGA<T>场不会被碰撞读取

[2026-06-06] [修复] [examples/bubble2d/bubble2d.cpp] [L84-90] 支持直接指定Gravity参数（非零→直接使用，为零→从Eo公式推导），解决大R导致|g|极小的问题

[2026-06-06] [修复] [examples/bubble2d/bubble2d.cpp] [L355-357] 添加phi钳位到[0,1]，消除BGKSource源项产生的数值噪点（黑色噪声点）

[2026-06-06] [诊断] [examples/bubble2d/bubble2d.cpp] [L93-108] 增加Eo_actual/Re_actual输出和Umax诊断，帮助判断参数一致性

[2026-06-06] [诊断] [examples/bubble2d/bubble2d.ini] [L7-9] MPI并行问题：部分进程数（如8）出现segfault是BlockGeometryHelper2D库级问题，与grid无法被total_blocks整除有关。256x512,BlockCellLen=64→32blocks，应使用np=1,2,4,8,16,32等整除32的值

[2026-06-06] [重构] [examples/bubble2d/bubble2d.cpp] [L67-114] 用5步参数设计法完全重写readParam：Step1-D, Step2-U_g(Ma), Step3-g=U_g²/D, Step4-ν=U_g*D/Re(η_h=ν*ρ_h,η_l=η_h/10), Step5-σ=Δρ*g*D²/Eo。新增CFL检查(U_g+cs<1.2)。移除了手动指定eta_h/eta_l/sigma/Gravity的过期机制

[2026-06-06] [创建] [examples/bubble2d/MPI_ISSUE.md] [全文] MPI并行SegFault问题说明文档：详细的block分配机制、崩溃原因、3种解决方案（调整BlockCellLen/安全进程数/单进程+OpenMP）

[2026-06-06] [更新] [examples/bubble2d/bubble2d.ini] [全文] 更新为自洽参数：Gravity=-5e-5→sigma_auto=0.0064(保持Eo=100), eta_h=0.253→Re_actual≈40, tau_ns=0.879

[2026-06-06] [诊断] [src/ff/ff2d.hh] [L26] 根因分析：BGKSource的Allen-Cahn源项source_i=fomega*w_i*(c_i·n)*4φ(1-φ)/W中，FF2D的n=∇φ/|∇φ|截断阈值ε=1e-10过小，导致体区域噪声（φ≈0.01）的|∇φ|≈0.01>>1e-10产生虚假法向量，正反馈→噪声放大。∇λ方案不可行（|∇λ|<<|n|使源项过弱，无法抵消BGK数值扩散）

[2026-06-06] [修复] [src/ff/ff2d.hh] [L26-35] FF2D中法向量梯度截断阈值：ε=1e-10→ε=0.005。体区|∇φ|<0.005→n=0→source=0，噪声源切断。界面中心|∇φ|≈1/W≈0.33>>0.005正常。同时修正分母：grad_mag+epsilon→grad_mag（已知>ε直接除）

[2026-06-06] [新增] [src/ff/ff2d.h] [L68-78] 添加FFChemPotentialGradient2D<CELL>结构体声明（∇λ计算，存储到NORMAL）。当前方案未启用（注释留着备查）。

[2026-06-06] [新增] [src/ff/ff2d.hh] [L74-101] 实现FFChemPotentialGradient2D::apply()：按D2Q9权重求λ邻居值的各向同性FD梯度∇λ，写入NORMAL字段。以注释块形式保留，未在主循环中激活。

[2026-06-06] [新增] [examples/bubble2d/bubble2d.cpp] [L324-328] 添加FFChemPotGradTask和FFChemPotGradSel模板（∇λ计算任务），以注释保留备用。

[2026-06-06] [新增] [examples/bubble2d/bubble2d.cpp] [L461-466] 主循环中FFChemPotGradSel调用注释块，保留∇λ→NORMAL覆盖方案的接入点，方便后续切换CH/AC源。

[2026-06-06] [修复] [examples/bubble2d/bubble2d.ini] [L16] Mobility: 0.1→0.2，τ_phi: 0.8→1.1，ω_phi: 1.25→0.909（过松弛→欠松弛，减少振荡，phi_min品质略有改善）

[2026-06-06] [重构] [src/ff/ff2d.h] [L26-44] 将ff::params全局模板变量（gravity/Beta/Kappa/twoPhase_rho_l/h/twoPhase_eta_l/h）替换为字段式设计：新增TwoPhaseParamsBase基类、TwoPhaseParams<T>参数结构体、TWOPHASEPARAMS<T>字段别名。移除ff::params命名空间。根因：全局模板变量不在MPI通信路径上，若仅在Rank0计算会导致其他Rank使用零初始化默认值（gravity=0,Beta=0），气泡不上升。

[2026-06-06] [重构] [src/ff/ff2d.hh] [L67-69,L121-124,L148-152] 更新FFChemPotential2D::apply()、FFGravityForce2D::apply()、FFRhoOmegaUpdate2D::apply()三个函子：从cell.template get<TWOPHASEPARAMS<T>>()读取参数结构体替代原来的params::Beta<T>/params::gravity<T>等全局模板变量访问。

[2026-06-06] [重构] [examples/bubble2d/bubble2d.cpp] [L52,L157-163,L222-224,L228-231,L399] 添加ff::TWOPHASEPARAMS<T>到PFFIELDS列表和PFInitValues初始值；readParam()中将参数写入tp_params结构体（替代ff::params::*赋值）；调试输出改用tp_params访问。编译通过无错误。

[2026-06-06] [删除] [src/ff/ff2d.h] [L107-128] 移除ff::CommunicatePhi<T>(lattice)和ff::CommunicateDerivedFields<T,D>(lattice)两个通信辅助函数。根因：与CA模块对比发现，ca/zhu_stefanescu2d.h的通信函数只通信自有字段（STATEBase/FSBase定义在ca/内），而ff/通信函数引用的PHI/GRAD/NORMAL/LAPLACIAN/CHEMICALPOTENTIAL字段全部定义在utils/alias.h中，不属于ff模块自有。ff2d本质是纯函子模块，不应越权拥有非自有字段的通信逻辑。

[2026-06-06] [恢复] [examples/bubble2d/bubble2d.cpp] [L426,L429-432,L437,L471-472] 恢复为逐字段手动Communicate()调用，与ff2d模块边界一致。编译通过无错误。

[2026-06-07] [重构] [src/utils/alias.h] [L158-159,L198-201,L270-272] 将LAPLACIANBase/CHEMICALPOTENTIALBase和LAPLACIAN<T>/CHEMICALPOTENTIAL<T>别名从alias.h移除。根因：这两字段仅被ff模块内部函子使用，应按CA模块模式定义为ff自有字段而非全局通用字段。

[2026-06-07] [重构] [src/ff/ff2d.h] [全文] 参照src/ca/zhu_stefanescu2d.h的完整设计模式重构ff2d.h：
Section 1: LaplacianBase/ChemPotentialBase/TwoPhaseParamsBase 自有字段基类
Section 2: LAPLACIAN<T>/CHEMICALPOTENTIAL<T>/TwoPhaseParams<T>/TWOPHASEPARAMS<T> 自有字段别名
Section 3: FFFIELDS<T>（自有）/ FFEXTERNALFIELDS<T,D>（外部：PHI/GRAD/NORMAL/INTERFACEWIDTH）/ FFFIELDPACK（自有+外部）/ ALLFF_FIELDS（展平）
Section 4: 函子声明（不变）
Section 5: 通信函数：CommunicateLAPLACIAN/CommunicateCHEMICALPOTENTIAL（逐cell通信）/ BroadcastAndSetParams（逐成员MPI_Bcast+InitValue）/ CommunicateAllSelfFields（批量：仅CHEMICALPOTENTIAL，LAPLACIAN无需ghost同步）

[2026-06-07] [重构] [examples/bubble2d/bubble2d.cpp] [L222-225,L432] PFFIELDS中LAPLACIAN/CHEMICALPOTENTIAL改为ff::LAPLACIAN/ff::CHEMICALPOTENTIAL；通信调用改为ff::CommunicateAllSelfFields<T>(PFLattice)。编译通过无错误。

[2026-06-07] [重构] [src/ff/ff2d.h] [全文] 取消TWOPHASEPARAMS<T>聚合字段，将7个参数拆为独立字段：
Type A (逐cell变化需ghost sync): LAPLACIAN, CHEMICALPOTENTIAL (不变)
Type B (均匀字段需MPI_Bcast): GRAVITY, BETA, KAPPA, RHO_L, RHO_H, ETA_L, ETA_H (新增7个自有Base+别名)
TwoPhaseParams<T>保留为纯数据传递结构体(非字段类型)。BroadcastAndSetParams→BroadcastAllParams(对7个TypeB字段逐成员广播+InitValue)。

[2026-06-07] [重构] [src/ff/ff2d.hh] [L67-68,L121-124,L148-152] 三个函子改为独立字段访问：cell.get<BETA<T>>() / cell.get<KAPPA<T>>() / pf_cell.get<GRAVITY<T>>() / pf_cell.get<RHO_L<T>>() 等，替代原来的cell.get<TWOPHASEPARAMS<T>>()结构体访问。

[2026-06-07] [重构] [src/ff/ff2d.h] [L39-58] 将7个Type B字段(GRAVITY/BETA/KAPPA/RHO_L/RHO_H/ETA_L/ETA_H)从GenericField<GenericArray<T>>改为Data<T>。根因：这些是块级常量（同CONSTFORCE机制），无需逐cell数组+ghost同步。Data<T>的isField=false意味着：无ghost层、无Communicate开销、每个block只存一个T值。cell.get<FIELD<T>>()接口不变。

[2026-06-07] [重构] [examples/bubble2d/bubble2d.cpp] [L49-52,L156-157,L228-231,L392] 移除tp_params全局变量；PFInitValues改用独立变量(gravity,Beta,Kappa,rho_l,rho_h,eta_l,eta_h)；BroadcastAllParams调用改为传递7个独立变量；调试打印改为使用独立变量。编译通过无错误。
