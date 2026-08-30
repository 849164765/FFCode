// kh2dStd.cpp — 2D Kelvin-Helmholtz instability in horizontal magnetic field (Phase field + NS + Magnetic)
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

// KH 论文图像时间校准系数。由 kh2dStd 与文献 Fig.31 对比确定:
// 文献绘制的 t* 约等于本代码 Eq.(84) 的 3 倍; 现在改为 ini 可配置,
// 输出同时打印 t*_eq 与 t*_plot, 便于用 Bom=0 的 Fig.31 图像重新精确拟合。
T KH_TSTAR_SCALE;

// ---- Simulation Parameters ----
int Ni, Nj;                    // 网格数 (文献Sec.F: 计算域[0,L]x[0,L], L=256)
T Cell_Len;                    // 格子物理长度
int BlockCellLen, Thread_Num;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (文献 Sec.F 第1098行, Eq.81-83)
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Re, We, Fr, U, Tau_ns;
T DeltaRho;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bom;
T PsiSolver_K;
int PsiSolver_Iter;   // 每流体步磁势子迭代次数 (MRT 伪时间求解; 3D 用 20)
// KH 专用磁压力项 F_p = -0.5*Hx^2*grad(mu) 的开关与系数
// (仅 kh2dStd/khMag3d/khMagDeepSeek3D 案例内启用; src 的 Kelvin 力保持不动)
int MagPressureEnabled;
T MagPressureFactor;
int MagForceFilterEnabled;   // 壁面/x 缝磁力尖峰过滤开关 (1=过滤, 与 khMag3d 一致)

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("kh2dStd.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  // 两相参数 (文献 Sec.F 第1098行)
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  Re=r.getValue<T>("Two_Phase","Re");      // 雷诺数 (Eq.81: Re=rho_h*U*L/eta_h)
  We=r.getValue<T>("Two_Phase","We");      // Weber数 (Eq.82: We=rho_h*L*U^2/sigma)
  Fr=r.getValue<T>("Two_Phase","Fr");      // Froude数 (Eq.83: Fr=U^2/(g*L))
  U=r.getValue<T>("Two_Phase","U");        // 参考速度(剪切速度) (Sec.F 第1098行: U=0.02)
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  KH_TSTAR_SCALE=r.getValue<T>("Simulation_Settings","KH_TSTAR_SCALE",T{3.0});
  MagPressureEnabled=r.getValue<int>("Simulation_Settings","MagPressure",1);
  MagPressureFactor=r.getValue<T>("Simulation_Settings","MagPressure_Factor",T{1.0});
  MagForceFilterEnabled=r.getValue<int>("Simulation_Settings","MagForce_Filter",1);
  // Magnetic (文献 Sec.F 第1107行, Eq.84)
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  Bom=r.getValue<T>("Magnetic_Field","Bom");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K",T{0.5});
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter",1);

  // 文献 Sec.F 第1039行: 计算域[0,2L]x[0,2L], L为半域长度
  // 文献 L=256 对应 512x512 网格; 此处 L=半域=Ni*Cell_Len/2
  // 若 Ni=512 则 L=256 (完全匹配文献); 若 Ni=256 则 L=128 (半尺度, 参数正确)
  T L = T(Ni) * Cell_Len * T(0.5);  // 参考长度 L = 半域 (文献 Sec.F: 域=2Lx2L)

  // 派生参数 (文献公式推导, LBM单位 mu0=1)
  // sigma = rho_h*L*U^2/We  (文献 Eq.82)
  sigma = rho_h * L * U * U / We;
  // gravity = U^2/(Fr*L)  (文献 Eq.83)
  gravity = U * U / (Fr * L);
  // eta_h = rho_h*U*L/Re  (文献 Eq.81)
  eta_h = rho_h * U * L / Re;
  // eta_l = eta_h  (假设: 文献未明确KH的eta_l, 因rho_l/rho_h=0.99近等密度)
  eta_l = eta_h;
  // H0 = sqrt(2*sigma*Bom/L)  (文献 Eq.84, mu0=1)
  H0 = std::sqrt(T(2.0) * sigma * Bom / L);
  // Beta = 12*sigma/W, Kappa = 3*W*sigma/2  (文献 Eq.14)
  Beta = T(12.0) * sigma / Interface_Width;
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  // Tau_phi = 3*Mobility + 0.5  (文献 Eq.22)
  Tau_phi = T(3.0) * Mobility + T(0.5); Omega_phi = T(1.0) / Tau_phi;
  // Tau_ns = 0.5 + eta_h/(rho_h*cs2)  (文献 Eq.33)
  Tau_ns = T(0.5) + eta_h / rho_h / LatSet::cs2;
  DeltaRho = rho_h - rho_l;

  T cs = std::sqrt(LatSet::cs2), Ma = U / cs;
  MPI_RANK(0) {
    printf("---- Kelvin-Helmholtz Instability in Horizontal Magnetic Field ----\n");
    printf("Mesh: %dx%d  BlockCellLen=%d  L=%.0f\n", Ni, Nj, BlockCellLen, L);
    printf("rho: l=%.4f h=%.4f  eta: l=%.6e h=%.6e  sigma=%.6e g=%.6e\n", rho_l, rho_h, eta_l, eta_h, sigma, gravity);
    printf("Re=%.0f We=%.0f Fr=%.1f U=%.4f W=%.1f M=%.3f tau_phi=%.3f tau_ns=%.5f\n", Re, We, Fr, U, Interface_Width, Mobility, Tau_phi, Tau_ns);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) H0=%.6f Bom=%.0f\n", chi_l, chi_h, mu_l, mu_h, H0, Bom);
    printf("KH magnetic-pressure term: enabled=%d factor=%.3f\n", MagPressureEnabled, MagPressureFactor);
    printf("MF solver: MRT pseudo-time sub-iterations/step = %d\n", PsiSolver_Iter);
    printf("Time scale: t*_eq = step*sqrt(g/L); t*_plot = %.3f * t*_eq\n", KH_TSTAR_SCALE);
    printf("Ma=%.4f\n", Ma);
    printf("---------------------------------------------------------------\n");
  }
  // 修复 MPI_RANK 宏悬挂问题: 需用大括号将 MPI_RANK 与语句包裹,
  // 否则宏展开为 if(...) if(rank!=0){return;} {fprintf(...);} 导致 fprintf 无条件执行
  if(Ma > T(0.1)){
    MPI_RANK(0)
    fprintf(stderr,"[Warn] Ma=%.4f > 0.1\n", Ma);
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Kelvin-Helmholtz Instability in Horizontal Magnetic Field..."));
  readParam();

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.001),Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni,T(0.001),Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.001),T(1.0));
  UnitConvManager<T> ConvManager(&BaseConv); ConvManager.Check_and_Print();

  // -- geometry (方形域 [0,L]x[0,L]) --
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
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,PSI_K<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi (KH interface, flipped tanh) --
  // 初始相场 phi (文献 Eq.78-79, 翻转tanh使底层phi=1铁磁流体)
  // y(x) = L + 0.1*L*cos(2*pi*x/L)  (文献 Eq.78, 适配用户域[0,L]x[0,L])
  // phi = 0.5 + 0.5*tanh(2*(y(x)-y)/W)  (文献 Eq.79翻转版, 底层phi=1)
  T L_phys = T(Ni) * Cell_Len * T(0.5);
  T W_phys = Interface_Width * Cell_Len;
  auto& phiField = PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b); auto& bl=PFLattice.getBlockLat(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=my+T(j)*vs;
      for(int i=ov;i<bk.getNx()-ov;++i){
        T x=mx+T(i)*vs;
        // 文献 Eq.78 的 OCR 相位易误读。用户对照 Fig.30/31 确认:
        //   y(x) = L - 0.1*L*cos(2*pi*x/L)
        // 即在 [0,2L] 内出现 2 个完整波峰、中间 1 个完整波谷、左右两侧
        // 各半个波谷(而不是旧代码的中间波峰+两侧波谷)。相位反转。
        T y_interface = L_phys*T(1.0) - T(0.1)*L_phys*std::cos(T(2.0)*M_PI*x/L_phys);
        // phi = 0.5 + 0.5*tanh(2*(y_interface-y)/W) (文献 Eq.79翻转, 底层phi=1)
        T phi = T(0.5) + T(0.5)*std::tanh(T(2.0)*(y_interface - y)/W_phys);
        bPhi.get(j*pr[1]+i)=phi;
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

  // init NS pops (p=0, u=shear flow) 文献 Eq.80: u=(0.5-phi)*U
  // phi=1(下层铁磁): u=-0.5*U (向左); phi=0(上层有机): u=+0.5*U (向右)
  T pz=0;
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=NSLattice.getBlockLat(b); auto& bPhi=phiField.getBlockField(b);
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j)
      for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL c(id,bl);
        T phi=bPhi.get(id);
        T ux=(T(0.5)-phi)*U;  // 文献 Eq.80 剪切流
        Vector<T,2> uz{ux,T(0)};
        for(unsigned k=0;k<LatSet::q;++k){
          T uc=uz*latset::c<LatSet>(k);
          c[k]=latset::w<LatSet>(k)*(pz+LatSet::InvCs2*uc+uc*uc*T(0.5)*LatSet::InvCs4-LatSet::InvCs2*T(0));
        }
      }
  }

  // init MF: psi = -H0*x (水平均匀场, H=(H0,0))
  // 注意: bubbleMag2d用psi=-H0*y(垂直场), 此处改为psi=-H0*x(水平场)
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0];
      for(int j=0;j<bk.getNy();++j){
        for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs, psi=-H0*x;  // 水平场: psi=-H0*x
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
  // MF 周期边界 (左右壁): 用于扰动 δψ=ψ+H0*x 的周期传播
  // 注意: 水平场 ψ=-H0*x 非周期, 需在 Apply 后施加扭曲校正 (见 step 0c)
  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
#ifdef MPI_ENABLED
  NS_Per.SetupMPI(GeoHelper); PF_Per.SetupMPI(GeoHelper); MF_Per.SetupMPI(GeoHelper);
#endif

  // -- MF wrapped-ghost field sync plans (MPI-aware) --
  // x 周期缝的 per-cell MF 场 (HX/HY/HMAG) 必须镜像对侧物理缝列。
  // MF_Per.Apply() 只复制 POP + GenericRho(对 MFCELL 解析为 PHI), 所以
  // PSI/HX/HY/HMAG 的卷绕 ghost 列不会被框架交换; 旧代码让这些 ghost 停在
  // 初始值 0, MFComputeH2D / MFMagneticForce2D 在缝上读取陈旧的邻居,
  // 产生虚假 ∇|H| (~0.5*H0^2) 尖峰 -> 左右边界伪影。这里预计算跨 rank/
  // 同 rank 的缝伙伴, 在每次 H 计算后重同步 HX/HY/HMAG。
  // (PSI 因水平场 ψ=-H0*x 有斜坡, 缝 ghost 仍由 step 0c 的扭曲周期 BC 单独
  //  重建, 不走本同步——H 场本身在 x 方向是周期的, 直接镜像即可。)
  struct MFGhostSyncPlan {
    std::size_t ghostBlockIdx; // local edge block holding the wrapped ghost
    std::size_t srcBlockIdx;   // local partner (same rank)
    int ghostCol;              // wrapped-ghost column (0 or nx-1)
    int seamCol;               // MY physical seam column sent to the partner (1 or nx-2)
    int srcCol;                // partner's opposite-edge column copied into ghostCol (nx-2 or 1)
    int ny;
    bool crossRank;
    int partnerRank, sendTag, recvTag; // crossRank only
  };
  std::vector<MFGhostSyncPlan> MFGhostSyncPlans;
  {
    const T xLeft = T{0};
    const T xRight = T(Ni) * Cell_Len;
    auto edgeOf = [&](const auto& bk, bool& atL, bool& atR) {
      atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
      atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
    };
    for (std::size_t b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& bk = Geo.getBlock(b);
      bool atL, atR; edgeOf(bk, atL, atR);
      if (!atL && !atR) continue;
      MFGhostSyncPlan pl;
      pl.ghostBlockIdx = b;
      pl.ghostCol = atL ? 0 : bk.getNx() - 1;
      pl.seamCol = atL ? 1 : bk.getNx() - 2;   // my column that wraps to the partner
      pl.srcCol = atL ? bk.getNx() - 2 : 1;    // partner's column that wraps to me
      pl.ny = bk.getNy();
      pl.crossRank = false;
      std::size_t partner = Geo.getBlockNum();
      for (std::size_t bp = 0; bp < Geo.getBlockNum(); ++bp) {
        const auto& bk2 = Geo.getBlock(bp);
        if (bk2.getMin()[1] != bk.getMin()[1]) continue;
        bool pL, pR; edgeOf(bk2, pL, pR);
        if (!((atL && pR) || (atR && pL) || (atL && atR && bp == b))) continue;
        partner = bp; break;
      }
      if (partner < Geo.getBlockNum()) {
        pl.srcBlockIdx = partner;
        MFGhostSyncPlans.push_back(pl);
      } else {
#ifdef MPI_ENABLED
        const auto& globalGeo = GeoHelper.getBlockGeometry();
        int partnerId = -1;
        for (std::size_t g = 0; g < globalGeo.getBlockNum(); ++g) {
          const auto& gb = globalGeo.getBlock(g);
          if (gb.getMin()[1] != bk.getMin()[1]) continue;
          bool gL, gR; edgeOf(gb, gL, gR);
          if (!((atL && gR) || (atR && gL))) continue;
          partnerId = gb.getBlockId(); break;
        }
        if (partnerId < 0) continue; // geometry invariant broken — skip
        pl.crossRank = true;
        pl.partnerRank = GeoHelper.whichRank(partnerId);
        pl.sendTag = bk.getBlockId();
        pl.recvTag = partnerId;
        MFGhostSyncPlans.push_back(pl);
#endif
      }
    }
  }

  // Execute the MFGhostSyncPlans for one per-cell field: direct column copies
  // for same-rank partners, non-blocking column exchange for cross-rank ones.
  auto SyncMFPeriodicGhosts = [&](auto& field) {
#ifdef MPI_ENABLED
    constexpr int MF_SYNC_TAG_BASE = 9500; // distinct from PERIODIC_TAG_BASE
    std::vector<std::vector<T>> sendBufs, recvBufs;
    std::vector<MPI_Request> sendReqs, recvReqs;
    std::vector<std::size_t> recvBlock;
    std::vector<int> recvCol, recvNy;
#endif
    for (const auto& pl : MFGhostSyncPlans) {
      auto& fg = field.getBlockField(pl.ghostBlockIdx);
      const auto& pr = Geo.getBlock(pl.ghostBlockIdx).getProjection();
      if (!pl.crossRank) {
        auto& fs = field.getBlockField(pl.srcBlockIdx);
        for (int j = 0; j < pl.ny; ++j)
          fg.get(j * pr[1] + pl.ghostCol) = fs.get(j * pr[1] + pl.srcCol);
      }
#ifdef MPI_ENABLED
      else {
        sendBufs.emplace_back(pl.ny);
        auto& sb = sendBufs.back();
        for (int j = 0; j < pl.ny; ++j) sb[j] = fg.get(j * pr[1] + pl.seamCol);
        MPI_Request rq;
        mpi().iSend(sb.data(), pl.ny, pl.partnerRank, &rq, MF_SYNC_TAG_BASE + pl.sendTag);
        sendReqs.push_back(rq);
        recvBufs.emplace_back(pl.ny);
        recvBlock.push_back(pl.ghostBlockIdx);
        recvCol.push_back(pl.ghostCol);
        recvNy.push_back(pl.ny);
        mpi().iRecv(recvBufs.back().data(), pl.ny, pl.partnerRank, &rq, MF_SYNC_TAG_BASE + pl.recvTag);
        recvReqs.push_back(rq);
      }
#endif
    }
#ifdef MPI_ENABLED
    if (!sendReqs.empty()) {
      MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
      for (std::size_t i = 0; i < recvReqs.size(); ++i) {
        MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
        auto& fg = field.getBlockField(recvBlock[i]);
        const auto& pr = Geo.getBlock(recvBlock[i]).getProjection();
        const auto& rb = recvBufs[i];
        for (int j = 0; j < recvNy[i]; ++j) fg.get(j * pr[1] + recvCol[i]) = rb[j];
      }
    }
#endif
  };

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
  vtmo::vtmWriter<T,2> MW("kh2dStd",Geo);
  MW.addWriterSet(PW,PS,VW,Dw,Fw,Hxw,Hyw);

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
    if(Bom>T{0}){
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: 设置磁场上下壁 Dirichlet 边界 (ψ=-H0*x, 维持外加水平场)
    // 左右壁改用扭曲周期 BC (见 step 0c), 允许扰动 δψ 周期传播.
    // 物理依据: 上下壁远离界面(y=0,2L), 场扰动衰减, ψ≈-H0*x 为良好近似;
    //           左右壁存在周期扰动 cos(2πx/L), 必须允许 δψ 周期穿越.
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minX=bk.getMin()[0],vs=bk.getVoxelSize();
        T minY=bk.getMin()[1],maxY=bk.getMax()[1];
        T H_global=T(Nj)*Cell_Len;
        // 下壁: ψ=-H0*x (精确解析值, 维持水平场)
        if(minY<Cell_Len*T(1.5)){
          for(int jj=0;jj<ov;++jj){
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              T psi_ghost=-H0*x_ghost;
              MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_ghost;
              bPsi.get(id)=psi_ghost;
            }
          }
        }
        // 上壁: ψ=-H0*x (精确解析值, 维持水平场)
        if(maxY>H_global-Cell_Len*T(1.5)){
          for(int jj=ny-ov;jj<ny;++jj){
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              T psi_ghost=-H0*x_ghost;
              MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_ghost;
              bPsi.get(id)=psi_ghost;
            }
          }
        }
      }
    }

    // ===== 磁势子迭代 (PsiSolver_Iter 次 MRT 伪时间松弛) =====
    // 旧代码每流体步只做 1 次 MRT 伪时间松弛, 界面快速演化时磁场常常
    // “追不上”界面 -> 磁力偏小。这里把 collision→周期BC→stream→PSI
    // 整段重复 PsiSolver_Iter 次, 让每步磁势更接近 div(mu*grad psi)=0。
    for(int sub=0; sub<PsiSolver_Iter; ++sub){
    // 0c: MF collision (direct iteration)
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
        }
      }
    }

    // 0c2: MF 扭曲周期边界 (Twisted Periodic BC for horizontal field)
    // 标准周期 BC 复制 g_k: ghost ← source (MF_Per.Apply())
    // 但水平场 ψ=-H0*x 非周期: ψ(0)=0, ψ(L_domain)=-H0*L_domain
    // 扭曲校正: 对扰动 δψ=ψ+H0*x 施加标准周期
    //   左壁 ghost (x<0): g_k += w_k*H0*L_domain,  ψ 重算
    //   右壁 ghost (x>=L_domain): g_k -= w_k*H0*L_domain,  ψ 重算
    MF_Per.Apply();
    {
      T L_domain=T(Ni)*Cell_Len;   // 全域宽度 = 2L
      T shift=H0*L_domain;          // 扭曲位移 = H0 * 2L
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minX=bk.getMin()[0],maxX=bk.getMax()[0],vs=bk.getVoxelSize();
        T W_global=T(Ni)*Cell_Len;
        T H_global=T(Nj)*Cell_Len;
        T minY=bk.getMin()[1],maxY=bk.getMax()[1];
        // 左壁 ghost: x<0, 加 shift (源来自右壁, ψ 应增加 H0*L_domain)
        // 仅处理非上下壁的 ghost (避免覆盖 Dirichlet BC)
        if(minX<T(0)+Cell_Len*T(1.5)){
          for(int jj=ov;jj<ny-ov;++jj){  // 排除上下壁 ghost
            for(int ii=0;ii<ov;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              if(x_ghost<T(0)){
                MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k){
                  c[k]+=latset::w<MFLatSet>(k)*shift;
                }
                // 从 g_k 重算 PSI (MF_Per 不更新 PSI 字段)
                T psi_new=T(0);
                for(unsigned k=0;k<MFLatSet::q;++k) psi_new+=c[k];
                bPsi.get(id)=psi_new;
              }
            }
          }
        }
        // 右壁 ghost: x>=L_domain, 减 shift (源来自左壁, ψ 应减少 H0*L_domain)
        if(maxX>W_global-Cell_Len*T(1.5)){
          for(int jj=ov;jj<ny-ov;++jj){  // 排除上下壁 ghost
            for(int ii=nx-ov;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              if(x_ghost>=W_global){
                MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k){
                  c[k]-=latset::w<MFLatSet>(k)*shift;
                }
                T psi_new=T(0);
                for(unsigned k=0;k<MFLatSet::q;++k) psi_new+=c[k];
                bPsi.get(id)=psi_new;
              }
            }
          }
        }
      }
    }

    // 0c3: 重新设置上下壁 Dirichlet BC (修复 MF_Per 覆盖的角落 ghost)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minX=bk.getMin()[0],vs=bk.getVoxelSize();
        T minY=bk.getMin()[1],maxY=bk.getMax()[1];
        T H_global=T(Nj)*Cell_Len;
        if(minY<Cell_Len*T(1.5)){
          for(int jj=0;jj<ov;++jj){
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              T psi_ghost=-H0*x_ghost;
              MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_ghost;
              bPsi.get(id)=psi_ghost;
            }
          }
        }
        if(maxY>H_global-Cell_Len*T(1.5)){
          for(int jj=ny-ov;jj<ny;++jj){
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii;
              T x_ghost=minX+T(ii)*vs;
              T psi_ghost=-H0*x_ghost;
              MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_ghost;
              bPsi.get(id)=psi_ghost;
            }
          }
        }
      }
    }

    // 0d: Stream
    MFLattice.NormalFullCommunicate();
    MFLattice.Stream();
    MFLattice.NormalFullCommunicate();

    // 0e: PSI = Σg_i (interior only, walls keep psi_bc)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j){
          for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
            T psi=0; for(unsigned k=0;k<MFLatSet::q;++k)psi+=c[k];
            bPsi.get(id)=psi;
          }
        }
      }
    }
    CommunicatePSI<T>(MFLattice);
    }  // end PsiSolver_Iter sub-iterations

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
    // 左右周期缝的 H 场 ghost 同步: 修复缝上陈旧 H ghost 造成的
    // 虚假 ∇|H| 与左右边界伪影 (与 rosenMag3d/bubbleMag3d 的
    // MFGhostSyncPlans 同源)。PSI 缝不回同步, 由 0c 扭曲周期 BC 单独处理。
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>());
    }  // end if(Bom>T{0})

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic force (if chi>0 and H0>0)
    // Kelvin term (src): F_k = chi*|H|*grad|H|   (paper Eq.8, mu0=1)
    // KH-only magnetic-pressure supplement: F_p = -0.5*Hx^2*grad(mu)
    //   grad(mu) = (mu_h-mu_l)*grad(phi). 与 khMag3d/khMagDeepSeek3d 完全一致;
    //   src 不动 (bubbleMag/rosenMag 仍只用 Kelvin 力).
    if(Bom>T{0}){
      // 上下壁 Dirichlet 磁势会在壁面 ghost 处留下 Hmag 伪梯度, 使 Kelvin 力
      // 在最靠壁的 1~2 行出现数值尖峰(与界面无关)。与 khMag3d 一样, 在壁面
      // 16 个格子宽(且不超过 10% 域高)的带内不施加磁力, 诊断也不统计该带。
      const T H_global=T(Nj)*Cell_Len;
      const T wall_band=std::min(T{16.0}, H_global*T{0.1});
      T Hmag_min=T{1e30}, Hmag_max=T{0}, Hmag_sum=T{0}; long long Hmag_cnt=0;
      T Fk_max=T{0}, Fk_sum=T{0}, Fp_max=T{0}, Fp_sum=T{0}, Ft_max=T{0}, Ft_sum=T{0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap(); T minY=bk.getMin()[1], vs=bk.getVoxelSize();
        for(int j=ov;j<bk.getNy()-ov;++j){
          T y=minY+T(j-ov)*vs;
          if(MagForceFilterEnabled && (y<=wall_band || y>=H_global-wall_band)) continue;
          for(int i=ov;i<bk.getNx()-ov;++i){
            // 注: 左右缝的 H ghost 已在 Phase0 末尾由 SyncMFPeriodicGhosts
            // 更新, 这里不再跳过缝列 (避免损失缝上真实的磁力)。
            std::size_t id=j*pr[1]+i;
            PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
            // Record force before magnetic force
            T fx_before=ns.template get<FORCE<T,2>>()[0];
            T fy_before=ns.template get<FORCE<T,2>>()[1];
            MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
            T fkx=ns.template get<FORCE<T,2>>()[0]-fx_before;
            T fky=ns.template get<FORCE<T,2>>()[1]-fy_before;

            // KH 专用磁压力项 (切向场下抑制波浪的主项)
            T fpx=T{0}, fpy=T{0};
            if(MagPressureEnabled){
              T Hx=mf.template get<HX<T>>();
              const auto& grad_phi=pf.template get<GRAD<T,2>>();
              T dmu=mf.template get<MU_H<T>>()-mf.template get<MU_L<T>>();
              T pref=MagPressureFactor*(-T{0.5})*Hx*Hx*dmu;
              fpx=pref*grad_phi[0];
              fpy=pref*grad_phi[1];
              auto& F=ns.template get<FORCE<T,2>>();
              F[0]+=fpx; F[1]+=fpy;
            }

            // Diagnostics
            T Hmag=mf.template get<HMAG<T>>();
            Hmag_min=std::min(Hmag_min,Hmag);
            Hmag_max=std::max(Hmag_max,Hmag);
            Hmag_sum+=Hmag; Hmag_cnt++;
            T fk=std::sqrt(fkx*fkx+fky*fky);
            T fp=std::sqrt(fpx*fpx+fpy*fpy);
            T ft=std::sqrt((fkx+fpx)*(fkx+fpx)+(fky+fpy)*(fky+fpy));
            Fk_max=std::max(Fk_max,fk); Fk_sum+=fk;
            Fp_max=std::max(Fp_max,fp); Fp_sum+=fp;
            Ft_max=std::max(Ft_max,ft); Ft_sum+=ft;
          }
        }
      }
#ifdef MPI_ENABLED
      if(mpi().getSize()>1){
        double locmax[5]={Hmag_max,-Hmag_min,Fk_max,Fp_max,Ft_max}, globmax[5]={0,0,0,0,0};
        MPI_Allreduce(locmax,globmax,5,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        Hmag_max=globmax[0]; Hmag_min=-globmax[1]; Fk_max=globmax[2]; Fp_max=globmax[3]; Ft_max=globmax[4];
        double locsum[5]={Hmag_sum,Fk_sum,Fp_sum,Ft_sum,double(Hmag_cnt)}, globsum[5]={0,0,0,0,0};
        MPI_Allreduce(locsum,globsum,5,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        Hmag_sum=globsum[0]; Fk_sum=globsum[1]; Fp_sum=globsum[2]; Ft_sum=globsum[3]; Hmag_cnt=static_cast<long long>(globsum[4]);
      }
#endif
      if(mpi().getRank()==0){
        if(t()%OutputStep==0 || t()<3){
          T nrm=T{1}/T(Hmag_cnt>0 ? Hmag_cnt : 1);
          printf("[KH2D_MAG] step=%d Hmag: min=%.6e max=%.6e mean=%.6e range=%.6e (H0=%.6e)\n",
                 t(),Hmag_min,Hmag_max,Hmag_sum*nrm,Hmag_max-Hmag_min,H0);
          printf("[KH2D_MAG] step=%d F_kelvin: max=%.6e mean=%.6e | F_press: max=%.6e mean=%.6e | F_total: max=%.6e mean=%.6e\n",
                 t(),Fk_max,Fk_sum*nrm,Fp_max,Fp_sum*nrm,Ft_max,Ft_sum*nrm);
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

    // 壁面phi约束: 下壁(y=0)phi=1铁磁, 上壁(y=L)phi=0有机
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T(1.5)) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=T(1);  // 下壁phi=1
        if(My>H_global-Cell_Len*T(1.5)){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=T(0);}  // 上壁phi=0
      }
      PFLattice.getField<PHI<T>>().Communicate();
    }

    // wall phi ghost fix (下壁ghost=1, 上壁ghost=0)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T(1.5)) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T(1);   // 下壁ghost=1
        if(My>H_global-Cell_Len*T(1.5)) for(int j=ny-ov;j<ny;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T(0);  // 上壁ghost=0
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
      // ---- KH 2D debug: interface height/amplitude + effective Re ----
      T hmin=T{1e30}, hmax=T{-1e30}; long ncross=0;
      {
        auto& phiF=PFLattice.getField<PHI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPhi=phiF.getBlockField(b);
          int ov=bk.getOverlap(); T minY=bk.getMin()[1],vs=bk.getVoxelSize();
          for(int i=ov;i<bk.getNx()-ov;++i){
            T yprev=minY+T(ov)*vs, pprev=bPhi.get(ov*pr[1]+i);
            for(int j=ov+1;j<bk.getNy()-ov;++j){
              T ycur=minY+T(j)*vs, pcur=bPhi.get(j*pr[1]+i);
              if((pprev-T{0.5})*(pcur-T{0.5})<=T{0}){
                T h=yprev+(T{0.5}-pprev)*(ycur-yprev)/(pcur-pprev);
                if(h<hmin)hmin=h; if(h>hmax)hmax=h; ++ncross;
              }
              yprev=ycur; pprev=pcur;
            }
          }
        }
      }
#ifdef MPI_ENABLED
      if(mpi().getSize()>1){
        double loc[3]={hmax,-hmin,double(ncross)}, glob[3];
        MPI_Allreduce(loc,glob,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        hmax=glob[0]; hmin=-glob[1];
        long long gcross=0; MPI_Allreduce(&ncross,&gcross,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD); ncross=long(gcross);
      }
#endif
      T L_phys=T(Ni)*Cell_Len*T{0.5};
      T tstar_eq=T(t())*std::sqrt(gravity/L_phys);              // Eq.(84) 原定义
      T tstar=T(t())*KH_TSTAR_SCALE*std::sqrt(gravity/L_phys);  // 论文图像校准: t*_plot = k*t*_eq
      if(mpi().getRank()==0){
        printf("[KH2D_DBG] step=%d t*_eq=%.4f t*_plot=%.4f  hmin=%.3f hmax=%.3f amp=%.3f  ncross=%ld\n",
               t(),tstar_eq,tstar,hmin,hmax,(hmax-hmin),ncross);
        printf("[KH2D_DBG] L_ref=%.1f U0=%.4f eta_l=%.6e eta_h=%.6e sigma=%.6e g=%.6e omega_ns=%.5f tstar_scale=%.3f\n",
               L_phys,U,eta_l,eta_h,sigma,gravity,T{1}/Tau_ns,KH_TSTAR_SCALE);
      }
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
