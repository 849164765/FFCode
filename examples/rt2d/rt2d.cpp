// rt2d.cpp — 2D Rayleigh-Taylor instability in ferrofluid (Phase field + NS + Magnetic)
// 依据论文 Section E (行917–1029): 均匀磁场下 Rayleigh–Taylor 不稳定性算例,
// 重相(铁磁流体, phi=1)在上, 轻相(非磁介质, phi=0)在下, 界面 y=2L+0.1L*cos(2πx/L).
// 复用 bubbleMag2d 的三场格子与耦合框架, 仅以下 4 处不同:
//   ① readParam 参数读取与推导 (RT 无量纲数 -> LBM 物理参数)
//   ② 计算域尺寸 256×1024 (由 ini 读入)
//   ③ 初始相场: 余弦扰动平界面 (重相在上)
//   ④ 壁面相场边界: 顶=1 底=0 (区别于 bubble 顶底均=1)
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

// ---- Simulation Parameters ----
int Ni, Nj;
T Cell_Len;
int BlockCellLen, Thread_Num;

// RT dimensionless params (论文 Section E)
T L, sqrt_gL, At, Re_num, Ca, W_rt, ViscRatio;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (推导)
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, U0, Tau_ns, DeltaRho;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bom, MF_Epsilon;
int MF_SubSteps;  // MF子迭代次数: 补偿标准LBM流动无法像论文Lax-Wendroff(式39)那样加速
int MF_SOR_Iter;   // SOR求解器迭代次数
T MF_SOR_Omega;    // SOR松弛因子(1-2, 1.8-1.9最优)

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("rt2d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  // RT 无量纲参数 (论文行962)
  L = r.getValue<T>("RT_Params","L");
  sqrt_gL = r.getValue<T>("RT_Params","sqrt_gL");
  At = r.getValue<T>("RT_Params","At");
  Re_num = r.getValue<T>("RT_Params","Re");
  Ca = r.getValue<T>("RT_Params","Ca");
  W_rt = r.getValue<T>("RT_Params","W");
  ViscRatio = r.getValue<T>("RT_Params","ViscRatio");
  // 相场
  Mobility = r.getValue<T>("Phase_Field","Mobility");
  Interface_Width = r.getValue<T>("Phase_Field","Interface_Width");
  // 两相
  rho_l = r.getValue<T>("Two_Phase","rho_l");
  // 磁场
  chi_l = r.getValue<T>("Magnetic_Field","chi_l");
  chi_h = r.getValue<T>("Magnetic_Field","chi_h");
  mu_l = r.getValue<T>("Magnetic_Field","mu_l");
  mu_h = r.getValue<T>("Magnetic_Field","mu_h");
  Bom = r.getValue<T>("Magnetic_Field","Bom");
  // MF_Epsilon: 伪时间加速参数(论文式37: ε=1/μ₀=1, 论文行273)
  // 论文使用Lax-Wendroff流动(式39, γ>1)加速收敛, 本代码用标准LBM+子迭代替代
  // ε=1 → τ=0.5+3μ (轻相3.5, 重相6.5), 碰撞有效, 数值稳定
  MF_Epsilon = r.getValue<T>("Magnetic_Field","MF_Epsilon");
  // MF_SubSteps: 每个时间步的MF子迭代次数, 补偿标准LBM流动的慢收敛
  // 有效伪时间 = ε × SubSteps, 扩散长度 = √(2D×SubSteps) 需 > 界面宽度W
  MF_SubSteps = r.getValue<int>("Magnetic_Field","MF_SubSteps");
  // SOR求解器参数: 直接有限差分求解∇·(μ∇ψ)=0, 替代LBM碰撞-流动
  // 优点: 正确处理变磁导率μ(调和平均), 无条件稳定, 收敛快
  MF_SOR_Iter = r.getValue<int>("Magnetic_Field","MF_SOR_Iter");
  MF_SOR_Omega = r.getValue<T>("Magnetic_Field","MF_SOR_Omega");
  // 模拟设置
  MaxStep = r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep = r.getValue<int>("Simulation_Settings","OutputStep");

  // ===== 推导 LBM 物理参数 (按论文公式) =====
  // U0 = sqrt(g*L) (论文行960)
  U0 = sqrt_gL;
  // g = (sqrt(gL))^2 / L (由 U0=sqrt(gL) 反推)
  gravity = sqrt_gL * sqrt_gL / L;
  // rho_h = rho_l*(1+At)/(1-At) (论文式73反解: At=(rho_h-rho_l)/(rho_h+rho_l))
  rho_h = rho_l * (T(1.0)+At) / (T(1.0)-At);
  // eta_l = rho_l*U0*L/Re (论文式73后 Re定义: Re=rho_l*U0*L/eta_l)
  eta_l = rho_l * U0 * L / Re_num;
  // eta_h = ViscRatio * eta_l (用户确认粘度比100, 同bubble)
  eta_h = ViscRatio * eta_l;
  // sigma = eta_l*U0/Ca (论文式74: Ca=mu_l*U0/sigma, mu_l即eta_l动力粘度)
  sigma = eta_l * U0 / Ca;
  // H0 = sqrt(2*sigma*Bom/L) (论文式76反解: Bom=mu0*H0^2*L/(2*sigma), mu0=1)
  H0 = std::sqrt(T(2.0) * sigma * Bom / L);
  // Beta = 12*sigma/W (论文式14)
  Beta = T(12.0) * sigma / Interface_Width;
  // Kappa = 3*W*sigma/2 (论文式14)
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  // Tau_phi = 3*Mobility + 0.5 (论文式22: 1/s3^f = M_phi/(cs^2*dt)+0.5, dt=1, cs^2=1/3)
  Tau_phi = T(3.0) * Mobility + T(0.5); Omega_phi = T(1.0) / Tau_phi;
  // Tau_ns = 0.5 + eta_h/rho_h/cs^2 (论文式33, 取重相为参考)
  Tau_ns = T(0.5) + eta_h / rho_h / LatSet::cs2;
  // DeltaRho = rho_h - rho_l (密度差, 用于F_p/F_v)
  DeltaRho = rho_h - rho_l;

  // Mach 数检查
  T cs = std::sqrt(LatSet::cs2), Ma = U0 / cs;
  MPI_RANK(0){
    printf("==== Rayleigh-Taylor Instability in Ferrofluid ====\n");
    printf("Mesh: %dx%d (Lx4L)  BlockCellLen=%d\n", Ni, Nj, BlockCellLen);
    printf("RT params: L=%.0f sqrt(gL)=%.4f At=%.2f Re=%.0f Ca=%.2f W=%.1f ViscRatio=%.0f\n",
           L, sqrt_gL, At, Re_num, Ca, W_rt, ViscRatio);
    printf("Derived: g=%.3e U0=%.4f rho_l=%.4f rho_h=%.4f eta_l=%.4e eta_h=%.4e\n",
           gravity, U0, rho_l, rho_h, eta_l, eta_h);
    printf("         sigma=%.3e H0=%.5f (Bom=%.2f) Beta=%.3e Kappa=%.3e\n",
           sigma, H0, Bom, Beta, Kappa);
    printf("         tau_phi=%.4f tau_ns=%.4f DeltaRho=%.4f Ma=%.4f\n",
           Tau_phi, Tau_ns, DeltaRho, Ma);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) SOR: iter=%d omega=%.2f\n",
           chi_l, chi_h, mu_l, mu_h, MF_SOR_Iter, MF_SOR_Omega);
    printf("====================================================\n");
  }
  if(Ma > T(0.2)) MPI_RANK(0){ fprintf(stderr,"[Warn] Ma=%.4f > 0.2\n",Ma); }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Rayleigh-Taylor Instability in Ferrofluid..."));
  readParam();

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),T(1.0));
  UnitConvManager<T> ConvManager(&BaseConv); ConvManager.Check_and_Print();

  // -- geometry --
  AABB<T,2> domain({0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len)});
  AABB<T,2> left({T(-Cell_Len),0},{0,T(Nj*Cell_Len)});
  AABB<T,2> right({T(Ni*Cell_Len),0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(Ni,Nj,domain,Cell_Len,BlockCellLen);
  GeoHelper.CreateBlocks(8,16);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // -- flag --
  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);

  // -- NS lattice --
  using NSFIELDS=TypePack<DENSITY<T>,VELOCITY<T,2>,POP<T,LatSet::q>,FORCE<T,2>,OMEGA<T>,PRESSURE<T>>;
  T omega_ns=T{1}/Tau_ns;
  ValuePack NSI(T{1},Vector<T,2>{0,0},T{},Vector<T,2>{0,0},omega_ns,T{});
  using NSCELL=Cell<T,LatSet,NSFIELDS>;
  BlockLatticeManager<T,LatSet,NSFIELDS> NSLattice(Geo,NSI,BaseConv);

  // -- PF lattice --
  using PFFIELDS=TypePack<PHI<T>,POP<T,LatSet::q>,GRAD<T,2>,NORMAL<T,2>,INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>,ff::CHEMICALPOTENTIAL<T>,
    ff::GRAVITY<T>,ff::BETA<T>,ff::KAPPA<T>,
    ff::RHO_L<T>,ff::RHO_H<T>,ff::ETA_L<T>,ff::ETA_H<T>,ff::DELTARHO<T>>;
  using PFREF=TypePack<VELOCITY<T,2>>;
  using PFPACK=TypePack<PFFIELDS,PFREF>;
  ValuePack PFI(T{},T{},Vector<T,2>{0,0},Vector<T,2>{0,0},Interface_Width,
    T{},T{},gravity,Beta,Kappa,rho_l,rho_h,eta_l,eta_h,DeltaRho);
  using PFCELL=Cell<T,LatSet,ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T,LatSet,PFPACK> PFLattice(Geo,PFI,PFBaseConv,
    &NSLattice.getField<VELOCITY<T,2>>());

  ff::BroadcastAllParams<T>(PFLattice,rho_l,rho_h,eta_l,eta_h,gravity,Beta,Kappa);
  PFLattice.template getField<ff::DELTARHO<T>>().InitValue(DeltaRho);

  // -- MF lattice (D2Q5) --
  using MFFIELDS=TypePack<PSI<T>,OMEGA_PSI<T>,MU_PERCELL<T>,CHI_PERCELL<T>,
    HX<T>,HY<T>,HMAG<T>,POP<T,MFLatSet::q>,
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,MF_EPSILON<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0,MF_Epsilon);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0,MF_Epsilon);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi (RT: 重相在上 phi=1, 轻相在下 phi=0, 界面 y=2L+0.1L*cos(2pi*x/L)) --
  // 论文式(72): y = 2L + 0.1L*cos(2*pi*x/L)
  // 论文式(11): phi = 0.5 + 0.5*tanh(2*(y-y_int)/W)
  T W_phys = Interface_Width * Cell_Len;
  auto& phiField = PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=my+T(j)*vs;
      for(int i=ov;i<bk.getNx()-ov;++i){
        T x=mx+T(i)*vs;
        // 界面位置 (式72): y_int = 2L + 0.1L*cos(2*pi*x/L)
        T y_int = T(2.0)*L + T(0.1)*L * std::cos(T(2.0)*T(3.14159265358979323846)*x/L);
        // phi (式11): 重相(y>y_int) phi=1, 轻相(y<y_int) phi=0
        T phi = T(0.5) + T(0.5)*std::tanh(T(2.0)*(y - y_int)/W_phys);
        bPhi.get(j*pr[1]+i) = phi;
      }
    }
  }
  // init PF pops
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=PFLattice.getBlockLat(b); auto& bPhi=phiField.getBlockField(b);
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j)
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; PFCELL c(id,bl); T phi=bPhi.get(id);
        for(unsigned k=0;k<LatSet::q;++k) c[k]=latset::w<LatSet>(k)*phi*(T{1}+LatSet::InvCs2*T{0});
      }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // init NS pops (p=0, u=0)
  Vector<T,2> uz{0,0}; T pz=0;
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j)
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL c(id,bl);
        for(unsigned k=0;k<LatSet::q;++k){
          T uc=uz*latset::c<LatSet>(k);
          c[k]=latset::w<LatSet>(k)*(pz+LatSet::InvCs2*uc+uc*uc*T{0.5}*LatSet::InvCs4-LatSet::InvCs2*T{0});
        }
      }
  }

  // init MF: psi = 分段线性 (正确的变μ稳态解, 使H场在界面处产生跳变)
  // 论文式(5): ∇·(μ∇ψ)=0 的1D稳态解 (水平界面 y_int):
  //   轻相 (y<y_int, μ=μ_l): ψ = -H0*y,  → H = H0
  //   重相 (y≥y_int, μ=μ_h): ψ = -H0*y_int - (H0/μ_h)*(y-y_int), → H = H0/μ_h
  //   界面处 B_n=μH 连续: μ_l*H0 = μ_h*(H0/μ_h) = H0 ✓
  // 初始界面位置 y_int(x) = 2L + 0.1L*cos(2πx/L) (论文式72)
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs;
        for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs;
          T y_int = T(2.0)*L + T(0.1)*L * std::cos(T(2.0)*T(3.14159265358979323846)*x/L);
          T psi;
          if(y < y_int){
            psi = -H0 * y;  // 轻相: H = H0
          } else {
            psi = -H0 * y_int - (H0 / mu_h) * (y - y_int);  // 重相: H = H0/μ_h
          }
          std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
          for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi;
          bPsi.get(id)=psi;
        }
      }
    }
  }

  // -- BCs --
  using LM_NS=BlockLatticeManager<T,LatSet,NSFIELDS>;
  using LM_PF=BlockLatticeManager<T,LatSet,PFPACK>;
  using LM_MF=BlockLatticeManager<T,MFLatSet,MFPACK>;
  using FM=BlockFieldManager<FLAG,T,2>;

  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,LM_NS,FM>
    NS_BB("NS_BB",NSLattice,FlagFM,BouncebackFlag,VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>,LM_PF,FM>
    PF_BB("PF_BB",PFLattice,FlagFM,BouncebackFlag,VoidFlag);

  FixedPeriodicBoundaryManager<LM_NS,FM> NS_Per("NS_Per",NSLattice,FlagFM,PeriodicFlag,VoidFlag);
  NS_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  NS_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  FixedPeriodicBoundaryManager<LM_PF,FM> PF_Per("PF_Per",PFLattice,FlagFM,PeriodicFlag,VoidFlag);
  PF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  PF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
#ifdef MPI_ENABLED
  NS_Per.SetupMPI(GeoHelper); PF_Per.SetupMPI(GeoHelper); MF_Per.SetupMPI(GeoHelper);
#endif

  // -- PF tasks --
  using PFNT=tmp::Key_TypePair<BulkFlag,ff::FF2D<PFCELL>>;
  using PFLT=tmp::Key_TypePair<BulkFlag,ff::FFLaplacian2D<PFCELL>>;
  using PFCT=tmp::Key_TypePair<BulkFlag,ff::FFChemPotential2D<PFCELL>>;
  using PFSelN=TaskSelector<std::uint8_t,PFCELL,PFNT>;
  using PFSelL=TaskSelector<std::uint8_t,PFCELL,PFLT>;
  using PFSelC=TaskSelector<std::uint8_t,PFCELL,PFCT>;
  using PFColT=tmp::Key_TypePair<BulkFlag,
    collision::MRTSource<equilibrium::FirstOrder<PFCELL>,NORMAL<T,2>,true,true>>;
  using PFPerT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<PFCELL>>;
  using PFAll=tmp::TupleWrapper<PFColT,PFPerT>;
  using PFSel=tmp::TaskSelector<PFAll,std::uint8_t,PFCELL>;

  // -- NS tasks --
  using NSMT=tmp::Key_TypePair<BulkFlag,collision::MRTForce<NSCELL,FORCE<T,2>>>;
  using NSPT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<NSCELL>>;
  using NSAll=tmp::TupleWrapper<NSMT,NSPT>;
  using NSSel=tmp::TaskSelector<NSAll,std::uint8_t,NSCELL>;

  // -- Coupling tasks --
  using STT=tmp::Key_TypePair<BulkFlag,ff::FFSurfaceTension2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling STC(PFLattice,NSLattice);
  using GrT=tmp::Key_TypePair<BulkFlag,ff::FFGravityForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling GrC(PFLattice,NSLattice);
  using PrT=tmp::Key_TypePair<BulkFlag,ff::FFPreForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling PrC(PFLattice,NSLattice);
  using ViT=tmp::Key_TypePair<BulkFlag,ff::FFViscoForce2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling ViC(PFLattice,NSLattice);
  using RoT=tmp::Key_TypePair<BulkFlag,ff::FFRhoOmegaUpdate2D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling RoC(PFLattice,NSLattice);

  // MF coupling: PF→MF (coeff update)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs2D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // Writers
  vtmo::ScalarWriter PW("PHI",PFLattice.getField<PHI<T>>());
  vtmo::ScalarWriter PS("PSI",MFLattice.getField<PSI<T>>());
  vtmo::VectorWriter VW("Velocity",NSLattice.getField<VELOCITY<T,2>>());
  vtmo::ScalarWriter Dw("Density",NSLattice.getField<DENSITY<T>>());
  vtmo::VectorWriter Fw("Force",NSLattice.getField<FORCE<T,2>>());
  vtmo::ScalarWriter Hxw("HX",MFLattice.getField<HX<T>>());
  vtmo::ScalarWriter Hyw("HY",MFLattice.getField<HY<T>>());
  vtmo::ScalarWriter Hmw("HMAG",MFLattice.getField<HMAG<T>>());
  vtmo::vtmWriter<T,2> MW("rt2d",Geo);
  MW.addWriterSet(PW,PS,VW,Dw,Fw,Hxw,Hyw,Hmw);

  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();
  NS_Per.Apply(); PF_Per.Apply();

  PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
  PFLattice.getField<NORMAL<T,2>>().Communicate();
  PFLattice.getField<GRAD<T,2>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);
  RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(0,FlagFM);
  MW.WriteBinary(0);

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  Timer t; Timer ot;
  T H_global=T(Nj)*Cell_Len;

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b-0e: SOR直接求解器 for ∇·(μ∇ψ) = 0 (替代LBM碰撞-流动)
    // 优点: ①正确处理变磁导率μ(调和平均) ②无条件稳定 ③不受τ限制
    // 离散: ψ_new = (μ_E·ψ_E + μ_W·ψ_W + μ_N·ψ_N + μ_S·ψ_S) / (μ_E+μ_W+μ_N+μ_S)
    // μ_face = 2·μ_c·μ_n/(μ_c+μ_n) (调和平均, 正确处理界面跳变)
    // Red-Black SOR: 先更新红格(i+j偶), 再更新黑格(i+j奇), 可并行
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      auto& muF =MFLattice.getField<MU_PERCELL<T>>();
      // 壁面BC函数 (Dirichlet: 底壁ψ=-H0*y, 顶壁ψ=-H0·2L-(H0/μ_h)·(y-2L))
      auto setWallBC=[&](){
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bPsi=psiF.getBlockField(b);
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
          T minY=bk.getMin()[1],maxY=bk.getMax()[1],vs=bk.getVoxelSize();
          if(minY<Cell_Len*T{1.5})
            for(int jj=0;jj<=ov;++jj){ T y=minY+T(jj)*vs;
              for(int ii=0;ii<nx;++ii) bPsi.get(jj*pr[1]+ii)=-H0*y; }
          if(maxY>H_global-Cell_Len*T{1.5}){ T yi=T(2.0)*L;
            for(int jj=ny-ov-1;jj<ny;++jj){ T y=minY+T(jj)*vs;
              for(int ii=0;ii<nx;++ii) bPsi.get(jj*pr[1]+ii)=-H0*yi-(H0/mu_h)*(y-yi); }}
        }
      };
      setWallBC();
      CommunicatePSI<T>(MFLattice);
      // SOR迭代
      for(int iter=0;iter<MF_SOR_Iter;++iter){
        for(int color=0;color<2;++color){ // 0=红(i+j偶) 1=黑(i+j奇)
          for(int b=0;b<Geo.getBlockNum();++b){
            auto& bPsi=psiF.getBlockField(b); auto& bMu=muF.getBlockField(b);
            const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
            int ov=bk.getOverlap(),nx=bk.getNx(),ny=bk.getNy();
            T minY=bk.getMin()[1],maxY=bk.getMax()[1];
            bool is_bottom=(minY<Cell_Len*T{1.5}), is_top=(maxY>H_global-Cell_Len*T{1.5});
            for(int j=ov;j<ny-ov;++j){
              if(is_bottom&&j==ov) continue;      // 跳过底壁(Dirichlet)
              if(is_top&&j==ny-ov-1) continue;     // 跳过顶壁(Dirichlet)
              for(int i=ov;i<nx-ov;++i){
                if(((i+j)&1)!=color) continue;
                std::size_t id=j*pr[1]+i;
                T muc=bMu.get(id);
                T muE=bMu.get(id+1),muW=bMu.get(id-1);
                T muN=bMu.get((j+1)*pr[1]+i),muS=bMu.get((j-1)*pr[1]+i);
                // 调和平均 (正确处理界面μ跳变)
                T hE=2*muc*muE/(muc+muE),hW=2*muc*muW/(muc+muW);
                T hN=2*muc*muN/(muc+muN),hS=2*muc*muS/(muc+muS);
                T psi_new=(hE*bPsi.get(id+1)+hW*bPsi.get(id-1)
                          +hN*bPsi.get((j+1)*pr[1]+i)+hS*bPsi.get((j-1)*pr[1]+i))
                          /(hE+hW+hN+hS);
                bPsi.get(id)=(T{1}-MF_SOR_Omega)*bPsi.get(id)+MF_SOR_Omega*psi_new;
              }
            }
          }
          setWallBC(); // 每次颜色pass后重设壁面BC
        }
      }
      CommunicatePSI<T>(MFLattice);
    }

    // 0f: Compute H = -∇ψ
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          MFComputeH2D<MFCELL>::apply(c);
        }
      }
    }
    CommunicateAllMFFields<T>(MFLattice);

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic force (if chi>0 and H0>0)
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
            MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
          }
        }
      }
    }

    GrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,GrT>>(t(),FlagFM);
    PrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,PrT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().Communicate();

    // ===== Phase B-C: PF + NS collision =====
    PFLattice.template ApplyInnerCellDynamics<PFSel>(FlagFM);
    PF_Per.Apply(); PFLattice.NormalFullCommunicate();

    ViC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,ViT>>(t(),FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSSel>(FlagFM);
    NS_Per.Apply(); NSLattice.NormalFullCommunicate();

    // ===== Phase D: Streaming =====
    PF_BB.Apply(t());
    { // PF Y ghost pop fix
      T Hg=T(Nj)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=PFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){PFCELL g(0*pr[1]+i,bl),w(ov*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
        if(My>Hg-Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){PFCELL g((ny-1)*pr[1]+i,bl),w((ny-1-ov)*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
      }
    }
    PFLattice.Stream(); PFLattice.NormalFullCommunicate();

    NS_BB.Apply(t());
    { // NS Y ghost pop fix
      T Hg=T(Nj)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){NSCELL g(0*pr[1]+i,bl),w(ov*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
        if(My>Hg-Cell_Len*T{1.5})
          for(int i=0;i<nx;++i){NSCELL g((ny-1)*pr[1]+i,bl),w((ny-1-ov)*pr[1]+i,bl); for(unsigned k=0;k<LatSet::q;++k)g[k]=w[k];}
      }
    }
    NSLattice.Stream(); NSLattice.NormalFullCommunicate();

    // ===== Phase E: Macro update =====
    // phi = sum(g_i)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=PFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j)
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; PFCELL c(id,bl);
            T pn=0; for(unsigned k=0;k<LatSet::q;++k)pn+=c[k];
            if(pn<T{0})pn=T{0}; if(pn>T{1})pn=T{1}; bP.get(id)=pn;
          }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi BC: 顶壁 phi=1 (重相/铁磁流体), 底壁 phi=0 (轻相/非磁介质)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        // 底壁 phi=0
        if(my<Cell_Len*T{1.5}) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=T{0};
        // 顶壁 phi=1
        if(My>H_global-Cell_Len*T{1.5}){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=T{1};}
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();
    // wall phi ghost fix (底ghost=0, 顶ghost=1)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{0};
        if(My>H_global-Cell_Len*T{1.5}) for(int j=ny-ov;j<ny;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{1};
      }
    }

    // gradients
    PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
    PFLattice.getField<NORMAL<T,2>>().Communicate();
    PFLattice.getField<GRAD<T,2>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // wall grad/chempot fix (PF)
    {
      auto& gF=PFLattice.getField<GRAD<T,2>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bG=gF.getBlockField(b); int nx=bk.getNx(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=ov;i<nx-ov;++i) bG.get(ov*pr[1]+i)[1]=bG.get((ov+1)*pr[1]+i)[1];
        if(My>H_global-Cell_Len*T{1.5}){int jj=bk.getNy()-1-ov;for(int i=ov;i<nx-ov;++i)bG.get(jj*pr[1]+i)[1]=bG.get((jj-1)*pr[1]+i)[1];}
      }
    }
    {
      auto& cF=PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bC=cF.getBlockField(b); int nx=bk.getNx(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=ov;i<nx-ov;++i) bC.get(ov*pr[1]+i)=(T{4}*bC.get((ov+1)*pr[1]+i)-bC.get((ov+2)*pr[1]+i))/T{3};
        if(My>H_global-Cell_Len*T{1.5}){int jj=bk.getNy()-1-ov;for(int i=ov;i<nx-ov;++i)bC.get(jj*pr[1]+i)=(T{4}*bC.get((jj-1)*pr[1]+i)-bC.get((jj-2)*pr[1]+i))/T{3};}
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // NS macro
    {
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& rF=NSLattice.getField<DENSITY<T>>(); auto& pF=NSLattice.getField<PRESSURE<T>>();
        auto& vF=NSLattice.getField<VELOCITY<T,2>>(); auto& fF=NSLattice.getField<FORCE<T,2>>();
        auto& bR=rF.getBlockField(b); auto& bP=pF.getBlockField(b);
        auto& bV=vF.getBlockField(b); auto& bF=fF.getBlockField(b);
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j)
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; NSCELL c(id,bl);
            T p=0,ux=0,uy=0;
            for(unsigned k=0;k<LatSet::q;++k){p+=c[k];ux+=latset::c<LatSet>(k)[0]*c[k];uy+=latset::c<LatSet>(k)[1]*c[k];}
            T rho=bR.get(id); auto F=bF.get(id);
            bP.get(id)=p; bV.get(id)=Vector<T,2>{ux+T{0.5}*F[0]/rho,uy+T{0.5}*F[1]/rho};
          }
      }
    }
    NSLattice.getField<VELOCITY<T,2>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    ++t; ++ot;
    if(t()%OutputStep==0){
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
    }
    // 诊断: 每500步打印H场统计, 监控MF求解器是否维持界面跳变
    if(t()%500==0){
      IF_MPI_RANK(0){
        T hhmax=0,hhmin=1e10,hhsum=0; int hcnt=0;
        T lhmax=0,lhmin=1e10,lhsum=0; int lcnt=0;
        auto& hF=MFLattice.getField<HMAG<T>>();
        auto& pF=PFLattice.getField<PHI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bH=hF.getBlockField(b); auto& bP=pF.getBlockField(b);
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          int ov=bk.getOverlap();
          for(int j=ov;j<bk.getNy()-ov;++j)
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=j*pr[1]+i;
              T h=bH.get(id), phi=bP.get(id);
              if(phi>T{0.9}){ hhsum+=h; hcnt++; if(h>hhmax)hhmax=h; if(h<hhmin)hhmin=h; }
              else if(phi<T{0.1}){ lhsum+=h; lcnt++; if(h>lhmax)lhmax=h; if(h<lhmin)lhmin=h; }
            }
        }
        T hh_mean=hhsum/hcnt, lh_mean=lhsum/lcnt;
        T ratio=hh_mean/lh_mean;
        printf("[Diag T=%d] H_heavy: mean=%.4e min=%.4e max=%.4e | H_light: mean=%.4e | ratio=%.4f (theory 0.500) max_ratio=%.1f\n",
               t(), hh_mean, hhmin, hhmax, lh_mean, ratio, hhmax/H0);
        fflush(stdout);
      }
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
