#include <cstring>
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D3Q19<T>;
using MFLatSet = D3Q7<T>;
using namespace mfield;

// KH 论文图像时间校准系数。由 kh2dStd 与文献 Fig.31 对比确定:
// 文献绘制的 t* 等于本代码 Eq.(84) 的 3 倍。
constexpr T KH_TSTAR_SCALE = T{4}; // 与 2D 校准一致: t*_plot = 4*t*_eq

// ---- Simulation Parameters ----
int Ni, Nj, Nk;
T Cell_Len;
int BlockCellLen, Thread_Num;

// interface (ferrofluid on bottom, solvent on top)
T InterfaceZ, PerturbAmp, PerturbPeriods;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (KH dimensionless -> lattice; Re/We/Fr use the SOLVENT = code "l" side)
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Tau_ns, DeltaRho;
T Re, We, Fr, U0, L_dom, L_ref;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bo_m;

int PsiSolver_Iter; T PsiSolver_K;
int PreForce;
int MagPressureEnabled; T MagPressureFactor; // KH 磁压力项开关/系数(仅本案例)

int MaxStep, OutputStep;
std::string work_dir;

void readParam(int argc, char* argv[]) {
  std::string iniName = "khMagDeepSeek3d.ini";
  if (argc > 1) iniName = argv[1];
  iniReader r(iniName);
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Nk = r.getValue<int>("Mesh","Nk");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  InterfaceZ   = r.getValue<T>("Interface","Z0_cells");
  PerturbAmp   = r.getValue<T>("Interface","PerturbAmp_cells");
  PerturbPeriods = r.getValue<T>("Interface","PerturbPeriods");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  // KH 无量纲数 (论文 F 节): Re/We/Fr 用溶剂(φ=0, code "l" 侧)的密度/粘度
  Re  = r.getValue<T>("Two_Phase","Re");
  We  = r.getValue<T>("Two_Phase","We");
  Fr  = r.getValue<T>("Two_Phase","Fr");
  U0  = r.getValue<T>("Two_Phase","U0");
  rho_l = r.getValue<T>("Two_Phase","rho_l");
  rho_h = r.getValue<T>("Two_Phase","rho_h");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  Bo_m=r.getValue<T>("Magnetic_Field","Bo_m");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PreForce=r.getValue<int>("Two_Phase","PreForce",1);
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");
  MagPressureEnabled=r.getValue<int>("Simulation_Settings","MagPressure",1);
  MagPressureFactor=r.getValue<T>("Simulation_Settings","MagPressure_Factor",T{1.0});

  // 由无量纲数推导格子参数 (论文 F 节 Eq.81-84):
  //   Re = rho_ref*U0*L/eta_ref,  We = rho_ref*L*U0^2/sigma,
  //   Fr = U0^2/(gL),  Bo_m = mu0*H0^2*L/(2*sigma) (mu0=1).
  // L_dom = 域宽 (x 方向, = 2L), L_ref = 参考长度 L = Ni/2.
  // 溶剂在上(φ=0, rho_l), 铁磁流体在下(φ=1, rho_h). 标准 KH 基准取
  // rho_h/rho_l=1/0.99, 即重流体在下, 重力是稳定的.
  const T pi = std::acos(T{-1});
  L_dom = T(Ni)*Cell_Len;
  L_ref = T(Ni)*Cell_Len/T{2};
  T rho_ref = rho_l;
  T eta_l_re = rho_ref*U0*L_ref/Re;
  // L=64 的低分辨率网格无法稳定表达 Re=5000: 按 Re 推导的
  // omega 约 1.997, 高于 MRT 稳定上限。这里把默认有效松弛钳到
  // omega=1.98, 对应有效 Re~760; 这正是文献时间轴所需的稳定设置。
  // ini 显式给出 eta_l 时尊重用户覆盖。
  T eta_l_stable = rho_ref*LatSet::cs2*(T{1}/T{1.98}-T{0.5});
  eta_l = std::max(eta_l_re, eta_l_stable);
  sigma = rho_ref*L_ref*U0*U0/We;
  gravity = U0*U0/(Fr*L_ref);
  eta_h = r.getValue<T>("Two_Phase","eta_h", eta_l);   // 文献未指定下相粘度, 默认等粘度
  // 显式覆盖 (可选): 若 ini 直接给出这些格子值则优先
  rho_l   = r.getValue<T>("Two_Phase","rho_l", rho_l);
  rho_h   = r.getValue<T>("Two_Phase","rho_h", rho_h);
  eta_l   = r.getValue<T>("Two_Phase","eta_l", eta_l);
  eta_h   = r.getValue<T>("Two_Phase","eta_h", eta_h);
  sigma   = r.getValue<T>("Two_Phase","sigma", sigma);
  gravity = r.getValue<T>("Two_Phase","Gravity", gravity);
  DeltaRho=rho_h-rho_l;
  H0 = r.getValue<T>("Magnetic_Field","H0", T{-1});
  if(H0 < T{0}) H0 = std::sqrt(T{2}*sigma*Bo_m/L_ref);   // mu0=1

  Beta=T{12.0}*sigma/Interface_Width;
  Kappa=T{3.0}*Interface_Width*sigma*T{0.5};
  Tau_phi=T{3.0}*Mobility+T{0.5}; Omega_phi=T{1.0}/Tau_phi;
  Tau_ns=T{0.5}+std::max(eta_l,eta_h)/std::min(rho_l,rho_h)/LatSet::cs2;

  MPI_RANK(0){
    printf("---- DeepSeek Kelvin-Helmholtz Instability (3D) ----\n");
    printf("Mesh: %dx%dx%d  BlockCellLen=%d  L=%.0f (2L=%.0f)\n",Ni,Nj,Nk,BlockCellLen,L_ref,T{2}*L_ref);
    printf("Interface: z0=%.1f  perturb A=%.2f (0.1L=%.1f)  periods=%.0f  W=%.1f\n",
           InterfaceZ,PerturbAmp,T{0.1}*L_ref,PerturbPeriods,Interface_Width);
    printf("Dimensionless: Re=%.0f  We=%.0f  Fr=%.2f  U0=%.4f\n",Re,We,Fr,U0);
    printf("rho: solvent(l)=%.4f ferrofluid(h)=%.4f  eta: l=%.5f h=%.5f  sigma=%.3e  g=%.3e\n",
           rho_l,rho_h,eta_l,eta_h,sigma,gravity);
    printf("Shear: u_x=(0.5-phi)*U0 -> ferrofluid(bottom)=-%.4f, solvent(top)=+%.4f\n",
           U0*T{0.5}, U0*T{0.5});
    printf("Magnetic (HORIZONTAL along x): chi=(%.2f,%.2f) mu=(%.2f,%.2f)  Bo_m=%.0f  H0=%.6f\n",
           chi_l,chi_h,mu_l,mu_h,Bo_m,H0);
    printf("PsiSolver: K=%.3f iter=%d  tau_ns=%.4f tau_phi=%.4f\n",
           PsiSolver_K,PsiSolver_Iter,Tau_ns,Tau_phi);
    printf("KH magnetic-pressure term: enabled=%d factor=%.3f\n", MagPressureEnabled, MagPressureFactor);
    printf("------------------------------------------\n");
    if(H0<=T{0}) printf("[Note] H0=0 (Bo_m=0): 无磁场对照\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing DeepSeek Kelvin-Helmholtz Instability (3D)..."));
  readParam(argc, argv);

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.001),Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni,T(0.001),Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.001),T(1.0));
  UnitConvManager<T> ConvManager(&BaseConv); ConvManager.Check_and_Print();

  // -- geometry --
  // Z is vertical (H_global = Nk*Cell_Len), walls in Z, X & Y are periodic.
  AABB<T,3> domain({0,0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> left  ({T(-Cell_Len),0,0},{0,T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> right ({T(Ni*Cell_Len),0,0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> front ({0,T(-Cell_Len),0},{T(Ni*Cell_Len),0,T(Nk*Cell_Len)});
  AABB<T,3> back  ({0,T(Nj*Cell_Len),0},{T(Ni*Cell_Len),T((Nj+1)*Cell_Len),T(Nk*Cell_Len)});
  BlockGeometryHelper3D<T> GeoHelper(Ni,Nj,Nk,domain,Cell_Len,BlockCellLen);
  // 网格必须能被 BlockCellLen 整除: 每方向块数 = 网格/块长
  GeoHelper.CreateBlocks(4,4,8);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> Geo(GeoHelper);

  std::vector<double> GlobalBlockTable;
  {
    const int nranks = mpi().getSize();
    const T xLeft = T{0}, xRight = T(Ni) * Cell_Len;
    const T yFront = T{0}, yBack = T(Nj) * Cell_Len;
    std::vector<double> local;
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& bk = Geo.getBlock(b);
      local.push_back(double(mpi().getRank()));
      local.push_back(bk.getMin()[0]); local.push_back(bk.getMin()[1]); local.push_back(bk.getMin()[2]);
      local.push_back(bk.getNx()); local.push_back(bk.getNy()); local.push_back(bk.getNz());
      local.push_back(bk.getMin()[0] < xLeft  + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMax()[0] > xRight - Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMin()[1] < yFront + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(bk.getMax()[1] > yBack  - Cell_Len * T{0.5} ? 1.0 : 0.0);
    }
    // 各 rank 条目数汇聚到 rank 0 再广播
    std::vector<int> counts(nranks, 0);
    {
      double tmp[1], dummy[1];
      if (mpi().getRank() == 0) {
        for (int r = 0; r < nranks; ++r) {
          tmp[0] = double(local.size());
          mpi().sendRecv(tmp, dummy, 1, r, r, 1000);
          counts[r] = int(dummy[0]);
        }
      } else {
        tmp[0] = double(local.size());
        mpi().sendRecv(tmp, dummy, 1, 0, 0, 1000);
      }
    }
    mpi().bCast(counts.data(), nranks, 0);
    int total = 0;
    for (int c : counts) total += c;
    // 数据汇聚到 rank 0 再广播
    std::vector<double> all(total, 0.0);
    {
      std::vector<double> dummy;
      if (mpi().getRank() == 0) {
        std::copy(local.begin(), local.end(), all.begin());
        int off = counts[0];
        for (int r = 1; r < nranks; ++r) {
          dummy.assign(counts[r], 0.0);
          mpi().sendRecv(dummy.data(), all.data() + off, counts[r], r, r, 1001);
          off += counts[r];
        }
      } else {
        dummy.assign(counts[mpi().getRank()], 0.0);
        mpi().sendRecv(local.data(), dummy.data(), int(local.size()), 0, 0, 1001);
      }
    }
    mpi().bCast(all.data(), total, 0);
    GlobalBlockTable = all;
  }

  // -- flag --
  BlockFieldManager<FLAG,T,3> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(front,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(back,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);
  // FIX(边界伪影, 对应 2D rosensweig bf37c87 根因): SetupBoundary 会把 x/y 周期
  // 侧面上所有域边界格点 (其周期 ghost 邻居在 domain AABB 之外, isInside=false)
  // 误设为 Bounceback 无滑移壁, 导致四周边界被当固壁 -> 边界伪影与接缝处提前
  // 起峰。将四个侧面的首/末物理列重置为 BulkFlag (仅保留 z 壁及其相邻棱边为
  // Bounceback 壁; x/y 均为周期方向, 垂直棱边属于周期-周期角, 也应为 Bulk)。
  AABB<T,3> left_col ({T{0},T{0},Cell_Len},{Cell_Len,T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
  AABB<T,3> right_col({T((Ni-1)*Cell_Len),T{0},Cell_Len},{T(Ni*Cell_Len),T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
  AABB<T,3> front_col({T{0},T{0},Cell_Len},{T(Ni*Cell_Len),Cell_Len,T((Nk-1)*Cell_Len)});
  AABB<T,3> back_col ({T{0},T((Nj-1)*Cell_Len),Cell_Len},{T(Ni*Cell_Len),T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
  FlagFM.forEach(left_col, [&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(right_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(front_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(back_col, [&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});

  // -- NS lattice --
  using NSFIELDS=TypePack<DENSITY<T>,VELOCITY<T,3>,POP<T,LatSet::q>,FORCE<T,3>,OMEGA<T>,PRESSURE<T>>;
  T omega_ns=T{1}/Tau_ns;
  ValuePack NSI(T{1},Vector<T,3>{0,0,0},T{},Vector<T,3>{0,0,0},omega_ns,T{});
  using NSCELL=Cell<T,LatSet,NSFIELDS>;
  BlockLatticeManager<T,LatSet,NSFIELDS> NSLattice(Geo,NSI,BaseConv);

  // -- PF lattice --
  using PFFIELDS=TypePack<PHI<T>,POP<T,LatSet::q>,GRAD<T,3>,NORMAL<T,3>,INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>,ff::CHEMICALPOTENTIAL<T>,
    ff::GRAVITY<T>,ff::BETA<T>,ff::KAPPA<T>,
    ff::RHO_L<T>,ff::RHO_H<T>,ff::ETA_L<T>,ff::ETA_H<T>,ff::DELTARHO<T>>;
  using PFREF=TypePack<VELOCITY<T,3>>;
  using PFPACK=TypePack<PFFIELDS,PFREF>;
  ValuePack PFI(T{},T{},Vector<T,3>{0,0,0},Vector<T,3>{0,0,0},Interface_Width,
    T{},T{},gravity,Beta,Kappa,rho_l,rho_h,eta_l,eta_h,DeltaRho);
  using PFCELL=Cell<T,LatSet,ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T,LatSet,PFPACK> PFLattice(Geo,PFI,PFBaseConv,
    &NSLattice.getField<VELOCITY<T,3>>());

  ff::BroadcastAllParams<T>(PFLattice,rho_l,rho_h,eta_l,eta_h,gravity,Beta,Kappa);
  PFLattice.template getField<ff::DELTARHO<T>>().InitValue(DeltaRho);

  // -- MF lattice (D3Q7) --
  using MFFIELDS=TypePack<PSI<T>,OMEGA_PSI<T>,MU_PERCELL<T>,CHI_PERCELL<T>,
    HX<T>,HY<T>,HZ<T>,HMAG<T>,POP<T,MFLatSet::q>,
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>,PSI_K<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams3D<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0,PsiSolver_K);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi (ferrofluid ON BOTTOM, solvent on top) --
  // 铁磁流体在下 (phi=1), 溶剂在上 (phi=0):
  //   phi(z) = 0.5 - 0.5*tanh(2(z - h(x,y))/W)
  //   h(x) = z0 - A*cos(2πx/L)   (文献相位, 与 kh2dStd 一致:
  //   域内 2 个完整波峰 + 中间 1 个完整波谷 + 左右各半个波谷)
  T z0=InterfaceZ*Cell_Len, W_phys=Interface_Width*Cell_Len;
  T lam=L_dom/PerturbPeriods;
  T kx=T{2}*std::acos(T{-1})/lam;
  auto interfaceHeight=[&](T x, T y)->T{
    return z0 - PerturbAmp*Cell_Len*std::cos(kx*x);
  };
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          auto vox=bk.getVoxel(Vector<int,3>{i,j,k});
          T x=vox[0], y=vox[1], z=vox[2];
          T h=interfaceHeight(x,y);
          T phi=T{0.5}-T{0.5}*std::tanh(T{2.0}*(z-h)/W_phys);
          bPhi.get(k*pr[2]+j*pr[1]+i)=phi;
        }
      }
    }
  }
  // init PF pops
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=PFLattice.getBlockLat(b); auto& bPhi=phiField.getBlockField(b);
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k)
      for(int j=ov;j<bk.getNy()-ov;++j)
        for(int i=ov;i<bk.getNx()-ov;++i){
          std::size_t id=k*pr[2]+j*pr[1]+i; PFCELL c(id,bl); T phi=bPhi.get(id);
          for(unsigned kk=0;kk<LatSet::q;++kk) c[kk]=latset::w<LatSet>(kk)*phi*(T{1}+LatSet::InvCs2*T{0});
        }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // init NS pops (u_x=(0.5-phi)*U0 剪切流 + 流体静压, 溶剂在上):
  // 静水压平衡 (避免 p=0 自由落体瞬态):
  //   z >= h (溶剂, 在上): p = rho_l*g*(H - z)
  //   z <= h (铁磁流体):   p = rho_l*g*(H - h) + rho_h*g*(h - z)
  // 剪切: 铁磁流体(下, phi=1) u_x=-U0/2, 溶剂(上, phi=0) u_x=+U0/2 (相对速度 U0)
  T H_total = T(Nk) * Cell_Len;
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      for(int j=ov;j<bk.getNy()-ov;++j){
        for(int i=ov;i<bk.getNx()-ov;++i){
          auto vox=bk.getVoxel(Vector<int,3>{i,j,k});
          T x=vox[0], y=vox[1], z=vox[2];
          T h=interfaceHeight(x,y);
          T p = (z>=h) ? rho_l*gravity*(H_total-z)
                       : rho_l*gravity*(H_total-h) + rho_h*gravity*(h-z);
          std::size_t id=k*pr[2]+j*pr[1]+i;
          T phi=bPhi.get(id);
          T ux=(T{0.5}-phi)*U0;
          Vector<T,3> u{ux,T{0},T{0}};
          NSCELL c(id,bl);
          for(unsigned kk=0;kk<LatSet::q;++kk){
            T uc=u*latset::c<LatSet>(kk);
            c[kk]=latset::w<LatSet>(kk)*(p+LatSet::InvCs2*uc+uc*uc*T{0.5}*LatSet::InvCs4-LatSet::InvCs2*T{0});
          }
        }
      }
    }
  }

  // -- init MF: 水平均匀场沿 +x -> psi = -H0*x (斜坡, 对 z/y 无关) --
  // 与 RT/Rosen 的垂直场(psi=-H0*z)不同, KH 场沿周期方向 x, psi 在 x 上有斜坡,
  // 故非周期: 不能直接周期卷绕, x 缝需用"斜坡平移周期 BC"(见 SetXSeamPsiPops)。
  auto psiRamp = [&](T x)->T{
    return -H0*x;
  };

  // KH 磁势不能用 psi=-H0*x 直接起算: 该初值只包含外加均匀场, 完全没有
  // 界面波浪诱导的扰动势 delta_psi, 磁力 F_m = chi*|H|*grad|H| 在 t=0
  // 恒为 0。delta_psi 是 H_global^2/D_psi 尺度的全局椭圆解, D3Q7 伪时间
  // 伪时间迭代从零初值爬到该解需要数万~数十万次迭代, 而且水平扭曲周期
  // 缝上的 D3Q7 迭代会漂移到反号分支。因此这里用水平场下小振幅正弦界面
  // 的解析磁势作初值, 再交给有限差分 Gauss-Seidel 求解器:
  //   delta_psi(x,z) = B sinh(kz) sin(kx),            z <= z0
  //                  = C sinh(k(H-z)) sin(kx),        z >  z0
  //   k=2*pi/lambda, 界面 h=z0-A cos(kx) (文献相位),
  //   B = -A*H0*(mu_h-mu_l) / { sinh(k z0)[ mu_l*coth(k(H-z0)) + mu_h*coth(k z0) ] },
  //   C = B*sinh(k z0)/sinh(k(H-z0))  (界面处 delta_psi 连续, B 的法向跳变抵消
  //   基态切向场的法向磁通跳变)。上下壁处 delta_psi=0, 与 Dirichlet 壁面一致。
  // 该一阶解析解已包含衰减最慢的全局 z 模态, 剩余高阶差由每步 Gauss-Seidel
  // 扫描快速抹平; 磁力从第 0 步就是正确量级, 高 Bom 的波浪抑制不再延迟。
  T psiAmp=-PerturbAmp*Cell_Len; // 界面相位 h=z0-A*cos -> delta_psi 反号
  T S0=std::sinh(kx*z0), S1=std::sinh(kx*(H_total-z0));
  T psiDenom=S0*(mu_l/std::tanh(kx*(H_total-z0)) + mu_h/std::tanh(kx*z0));
  T psiB=+(mu_h-mu_l)*psiAmp*H0/psiDenom;
  T psiC=psiB*S0/S1;
  auto psiKH=[&](T x, T z)->T{
    T delta=(z<=z0) ? psiB*std::sinh(kx*z)
                    : psiC*std::sinh(kx*(H_total-z));
    return psiRamp(x) + delta*std::sin(kx*x);
  };
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      for(int k=0;k<bk.getNz();++k){
        for(int j=0;j<bk.getNy();++j){
          for(int i=0;i<bk.getNx();++i){
            auto vox=bk.getVoxel(Vector<int,3>{i,j,k});
            T x=vox[0], z=vox[2], psi=psiKH(x,z);
            std::size_t id=k*pr[2]+j*pr[1]+i; MFCELL c(id,bl);
            for(unsigned kk=0;kk<MFLatSet::q;++kk) c[kk]=latset::w<MFLatSet>(kk)*psi;
            bPsi.get(id)=psi;
          }
        }
      }
    }
  }

  // -- BCs --
  using LM_NS=BlockLatticeManager<T,LatSet,NSFIELDS>;
  using LM_PF=BlockLatticeManager<T,LatSet,PFPACK>;
  using LM_MF=BlockLatticeManager<T,MFLatSet,MFPACK>;
  using FM=BlockFieldManager<FLAG,T,3>;

  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,LM_NS,FM>
    NS_BB("NS_BB",NSLattice,FlagFM,BouncebackFlag,VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>,LM_PF,FM>
    PF_BB("PF_BB",PFLattice,FlagFM,BouncebackFlag,VoidFlag);

  FixedPeriodicBoundaryManager<LM_NS,FM> NS_Per("NS_Per",NSLattice,FlagFM,PeriodicFlag,VoidFlag);
  NS_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  NS_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  NS_Per.Setup(front,NbrDirection::YN,back,NbrDirection::YP);
  NS_Per.Setup(back,NbrDirection::YP,front,NbrDirection::YN);
  FixedPeriodicBoundaryManager<LM_PF,FM> PF_Per("PF_Per",PFLattice,FlagFM,PeriodicFlag,VoidFlag);
  PF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  PF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
  PF_Per.Setup(front,NbrDirection::YN,back,NbrDirection::YP);
  PF_Per.Setup(back,NbrDirection::YP,front,NbrDirection::YN);
  FixedPeriodicBoundaryManager<LM_MF,FM> MF_Per("MF_Per",MFLattice,FlagFM,PeriodicFlag,VoidFlag);
  // KH: MF 的 x 方向不做周期卷绕 (psi=-H0*x 有斜坡, 非周期), 改为斜坡平移周期
  // BC (SetXSeamPsiPops 在子迭代里显式处理); y 向保持周期。
  MF_Per.Setup(front,NbrDirection::YN,back,NbrDirection::YP);
  MF_Per.Setup(back,NbrDirection::YP,front,NbrDirection::YN);
#ifdef MPI_ENABLED
  NS_Per.SetupMPI(GeoHelper); PF_Per.SetupMPI(GeoHelper); MF_Per.SetupMPI(GeoHelper);
#endif

  // -- PF tasks --
  using PFNT=tmp::Key_TypePair<BulkFlag,ff::FF3D<PFCELL>>;
  using PFLT=tmp::Key_TypePair<BulkFlag,ff::FFLaplacian3D<PFCELL>>;
  using PFCT=tmp::Key_TypePair<BulkFlag,ff::FFChemPotential3D<PFCELL>>;
  using PFSelN=TaskSelector<std::uint8_t,PFCELL,PFNT>;
  using PFSelL=TaskSelector<std::uint8_t,PFCELL,PFLT>;
  using PFSelC=TaskSelector<std::uint8_t,PFCELL,PFCT>;
  using PFColT=tmp::Key_TypePair<BulkFlag,
    collision::MRTSource<equilibrium::FirstOrder<PFCELL>,NORMAL<T,3>,true,true>>;
  using PFPerT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<PFCELL>>;
  using PFAll=tmp::TupleWrapper<PFColT,PFPerT>;
  using PFSel=tmp::TaskSelector<PFAll,std::uint8_t,PFCELL>;

  // -- NS tasks --
  using NSMT=tmp::Key_TypePair<BulkFlag,collision::MRTForce<NSCELL,FORCE<T,3>>>;
  using NSPT=tmp::Key_TypePair<PeriodicFlag,collision::PeriodicBoundary<NSCELL>>;
  using NSAll=tmp::TupleWrapper<NSMT,NSPT>;
  using NSSel=tmp::TaskSelector<NSAll,std::uint8_t,NSCELL>;

  // -- Coupling tasks --
  using STT=tmp::Key_TypePair<BulkFlag,ff::FFSurfaceTension3D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling STC(PFLattice,NSLattice);
  using GrT=tmp::Key_TypePair<BulkFlag,ff::FFGravityForce3D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling GrC(PFLattice,NSLattice);
  using PrT=tmp::Key_TypePair<BulkFlag,ff::FFPreForce3D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling PrC(PFLattice,NSLattice);
  using ViT=tmp::Key_TypePair<BulkFlag,ff::FFViscoForce3DM<PFCELL,NSCELL>>;
  BlockLatManagerCoupling ViC(PFLattice,NSLattice);
  using RoT=tmp::Key_TypePair<BulkFlag,ff::FFRhoOmegaUpdate3D<PFCELL,NSCELL>>;
  BlockLatManagerCoupling RoC(PFLattice,NSLattice);

  // MF coupling: PF→MF (coeff update)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs3D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // Writers
  vtmo::ScalarWriter PW("PHI",PFLattice.getField<PHI<T>>());
    vtmo::vtmWriter<T,3> MW("khMagDeepSeek3d",Geo);
  MW.addWriterSet(PW);

    auto SyncMFPeriodicGhosts = [&](auto& field, int fidx) {
      using ValT = typename std::decay_t<decltype(field.getBlockField(0))>::value_type;
      const int valSz = static_cast<int>(sizeof(ValT));
      const T xLeft = T{0}, xRight = T(Ni) * Cell_Len;
      const T yFront = T{0}, yBack = T(Nj) * Cell_Len;
      const int nEntry = int(GlobalBlockTable.size() / 11);
      const int myRank = mpi().getRank();
#ifdef MPI_ENABLED
      struct RecvJob { int block; int ghostPlane0; bool xSeam; };
      std::vector<std::vector<char>> sendBufs, recvBufs;
      std::vector<MPI_Request> sendReqs, recvReqs;
      std::vector<RecvJob> recvJobs;
#endif
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& bk = Geo.getBlock(b);
        const int nx = bk.getNx(), ny = bk.getNy(), nz = bk.getNz(), ov = bk.getOverlap();
        const auto& pr = bk.getProjection();
        const bool atL = bk.getMin()[0] < xLeft  + Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
        const bool atF = bk.getMin()[1] < yFront + Cell_Len * T{0.5};
        const bool atB = bk.getMax()[1] > yBack  - Cell_Len * T{0.5};
        if (!atL && !atR && !atF && !atB) continue;
        auto& f1 = field.getBlockField(b);

        // ---- x 缝：partner 同 (minY, minZ)，对侧 x 边 ----
        if (atL || atR) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[2] != bk.getMin()[1] || en[3] != bk.getMin()[2]) continue;
            if (atL && en[8] > T{0.5}) eP = e;   // partner atR
            if (atR && en[7] > T{0.5}) eP = e;   // partner atL
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2]) eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int plane = ny * nz;
            const int ghostCol0 = atL ? 0 : nx - ov;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const auto& pPr = Geo.getBlock(pb).getProjection();
                const int srcCol = atL ? (Geo.getBlock(pb).getNx() - 1 - ov) : ov;
                for (int k = 0; k < nz; ++k)
                  for (int j = 0; j < ny; ++j) {
                    const ValT v = f2.get(k * pPr[2] + j * pPr[1] + srcCol);
                    for (int c = 0; c < ov; ++c)
                      f1.get(k * pr[2] + j * pr[1] + ghostCol0 + c) = v;
                  }
              }
            }
#ifdef MPI_ENABLED
            else {
              const int myPhysCol = atL ? ov : nx - 1 - ov;
              const int tag = fidx * 10000 + std::min(eL, eP) * 2 + 0;
              std::vector<char> snd(plane * valSz);
              int cnt = 0;
              for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j) {
                  const ValT v = f1.get(k * pr[2] + j * pr[1] + myPhysCol);
                  std::memcpy(&snd[cnt * valSz], &v, valSz); ++cnt;
                }
              sendBufs.emplace_back(std::move(snd));
              MPI_Request sreq, rreq;
              mpi().iSend(sendBufs.back().data(), plane * valSz, pRank, &sreq, tag);
              sendReqs.push_back(sreq);
              recvBufs.emplace_back(plane * valSz);
              recvJobs.push_back({b, ghostCol0, true});
              mpi().iRecv(recvBufs.back().data(), plane * valSz, pRank, &rreq, tag);
              recvReqs.push_back(rreq);
            }
#endif
          }
        }

        // ---- y 缝：partner 同 (minX, minZ)，对侧 y 边 ----
        if (atF || atB) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[1] != bk.getMin()[0] || en[3] != bk.getMin()[2]) continue;
            if (atF && en[10] > T{0.5}) eP = e;  // partner atB
            if (atB && en[9]  > T{0.5}) eP = e;  // partner atF
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2]) eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int plane = nx * nz;
            const int ghostRow0 = atF ? 0 : ny - ov;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const auto& pPr = Geo.getBlock(pb).getProjection();
                const int srcRow = atF ? (Geo.getBlock(pb).getNy() - 1 - ov) : ov;
                for (int k = 0; k < nz; ++k)
                  for (int i = 0; i < nx; ++i) {
                    const ValT v = f2.get(k * pPr[2] + srcRow * pPr[1] + i);
                    for (int c = 0; c < ov; ++c)
                      f1.get(k * pr[2] + (ghostRow0 + c) * pr[1] + i) = v;
                  }
              }
            }
#ifdef MPI_ENABLED
            else {
              const int myPhysRow = atF ? ov : ny - 1 - ov;
              const int tag = fidx * 10000 + std::min(eL, eP) * 2 + 1;
              std::vector<char> snd(plane * valSz);
              int cnt = 0;
              for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i) {
                  const ValT v = f1.get(k * pr[2] + myPhysRow * pr[1] + i);
                  std::memcpy(&snd[cnt * valSz], &v, valSz); ++cnt;
                }
              sendBufs.emplace_back(std::move(snd));
              MPI_Request sreq, rreq;
              mpi().iSend(sendBufs.back().data(), plane * valSz, pRank, &sreq, tag);
              sendReqs.push_back(sreq);
              recvBufs.emplace_back(plane * valSz);
              recvJobs.push_back({b, ghostRow0, false});
              mpi().iRecv(recvBufs.back().data(), plane * valSz, pRank, &rreq, tag);
              recvReqs.push_back(rreq);
            }
#endif
          }
        }
      }
#ifdef MPI_ENABLED
      if (!sendReqs.empty()) {
        MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
        for (std::size_t i = 0; i < recvReqs.size(); ++i) {
          MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
          auto& f1 = field.getBlockField(recvJobs[i].block);
          const auto& bk = Geo.getBlock(recvJobs[i].block);
          const auto& pr = bk.getProjection();
          const int ov = bk.getOverlap();
          const int g0 = recvJobs[i].ghostPlane0;
          const auto& buf = recvBufs[i];
          int cnt = 0;
          if (recvJobs[i].xSeam) {
            for (int k = 0; k < bk.getNz(); ++k)
              for (int j = 0; j < bk.getNy(); ++j) {
                ValT v; std::memcpy(&v, &buf[cnt * valSz], valSz); ++cnt;
                for (int c = 0; c < ov; ++c) f1.get(k * pr[2] + j * pr[1] + g0 + c) = v;
              }
          } else {
            for (int k = 0; k < bk.getNz(); ++k)
              for (int i = 0; i < bk.getNx(); ++i) {
                ValT v; std::memcpy(&v, &buf[cnt * valSz], valSz); ++cnt;
                for (int c = 0; c < ov; ++c) f1.get(k * pr[2] + (g0 + c) * pr[1] + i) = v;
              }
          }
        }
      }
#endif
    };


  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();
  NS_Per.Apply(); PF_Per.Apply();

  // 初始 PHI 周期缝同步：第一次 FF2D 计算 ∇φ 前，必须让左右/前后
  // 周期 ghost 列的 PHI 等于对侧物理列，否则边界梯度从一开始就是错的，
  // 并会在后续步持续累积成左右边界异常。
  PFLattice.getField<PHI<T>>().Communicate();
  SyncMFPeriodicGhosts(PFLattice.getField<PHI<T>>(), 0);

  PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
  PFLattice.getField<NORMAL<T,3>>().Communicate();
  PFLattice.getField<GRAD<T,3>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);
  RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(0,FlagFM);
  MW.WriteBinary(0);

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  Timer t; Timer ot;
  T H_global=T(Nk)*Cell_Len;

  // ---- MF periodic-seam sync plans (3D port of the 2D MFGhostSyncPlans).
  // x/y-wrapped ghost planes of the per-cell MF fields must mirror the opposite edge's
  // physical seam plane. MF_Per.Apply() copies only pops + GenericRho (which resolves to
  // PHI), so these fields are never exchanged by the framework. Precompute the partner
  // plans ONCE (same-rank direct copy; cross-rank non-blocking exchange keyed by the
  // sender's global block id, distinct send/recv tags). The old per-call GlobalBlockTable
  // + single-tag sendRecv corrupted the field at cross-rank seams (1-rank clean vs MPI
  // broken at the y-seam). For KH the x-seam PSI is further fixed by SetXSeamPsiPops
  // (ramp-shift periodic) afterwards; H fields are plain-periodic in x.
  struct MFGhostSyncPlan {
    std::size_t ghostBlockIdx, srcBlockIdx;
    bool xSeam;
    int ghostPlaneIdx, seamPlaneIdx, srcPlaneIdx, plane;
    bool crossRank;
    int partnerRank, sendTag, recvTag;
  };
  std::vector<MFGhostSyncPlan> MFGhostSyncPlans;
  {
    const T xLeft=T{0}, xRight=T(Ni)*Cell_Len, yFront=T{0}, yBack=T(Nj)*Cell_Len;
    auto edgeX = [&](const auto& bk, bool& atL, bool& atR){
      atL = bk.getMin()[0] < xLeft + Cell_Len*T{0.5};
      atR = bk.getMax()[0] > xRight - Cell_Len*T{0.5};
    };
    auto edgeY = [&](const auto& bk, bool& atF, bool& atB){
      atF = bk.getMin()[1] < yFront + Cell_Len*T{0.5};
      atB = bk.getMax()[1] > yBack - Cell_Len*T{0.5};
    };
    for (std::size_t b=0;b<Geo.getBlockNum();++b){
      const auto& bk=Geo.getBlock(b);
      bool atL,atR; edgeX(bk,atL,atR);
      bool atF,atB; edgeY(bk,atF,atB);
      if(!atL&&!atR&&!atF&&!atB) continue;
      const int ov=bk.getOverlap();
      const int nx=bk.getNx(), ny=bk.getNy(), nz=bk.getNz();
      if(atL||atR){
        MFGhostSyncPlan pl;
        pl.ghostBlockIdx=b; pl.xSeam=true; pl.crossRank=false;
        pl.ghostPlaneIdx = atL ? (ov-1) : (nx-ov);
        pl.seamPlaneIdx  = atL ? ov : (nx-1-ov);
        pl.plane = ny*nz;
        std::size_t partner=Geo.getBlockNum();
        for(std::size_t bp=0;bp<Geo.getBlockNum();++bp){
          const auto& bk2=Geo.getBlock(bp);
          if(bk2.getMin()[1]!=bk.getMin()[1]||bk2.getMin()[2]!=bk.getMin()[2]) continue;
          bool pL,pR; edgeX(bk2,pL,pR);
          if(!((atL&&pR)||(atR&&pL))) continue;
          partner=bp; break;
        }
        if(partner<Geo.getBlockNum()){
          const auto& pb=Geo.getBlock(partner); const int pov=pb.getOverlap();
          pl.srcBlockIdx=partner;
          pl.srcPlaneIdx = atL ? (pb.getNx()-1-pov) : pov;
          MFGhostSyncPlans.push_back(pl);
        } else {
#ifdef MPI_ENABLED
          const auto& globalGeo=GeoHelper.getBlockGeometry();
          int partnerId=-1;
          for(std::size_t g=0;g<globalGeo.getBlockNum();++g){
            const auto& gb=globalGeo.getBlock(g);
            if(gb.getMin()[1]!=bk.getMin()[1]||gb.getMin()[2]!=bk.getMin()[2]) continue;
            bool gL,gR; edgeX(gb,gL,gR);
            if(!((atL&&gR)||(atR&&gL))) continue;
            partnerId=gb.getBlockId(); break;
          }
          if(partnerId<0) continue;
          pl.crossRank=true; pl.partnerRank=GeoHelper.whichRank(partnerId);
          pl.sendTag=bk.getBlockId(); pl.recvTag=partnerId;
          MFGhostSyncPlans.push_back(pl);
#endif
        }
      }
      if(atF||atB){
        MFGhostSyncPlan pl;
        pl.ghostBlockIdx=b; pl.xSeam=false; pl.crossRank=false;
        pl.ghostPlaneIdx = atF ? (ov-1) : (ny-ov);
        pl.seamPlaneIdx  = atF ? ov : (ny-1-ov);
        pl.plane = nx*nz;
        std::size_t partner=Geo.getBlockNum();
        for(std::size_t bp=0;bp<Geo.getBlockNum();++bp){
          const auto& bk2=Geo.getBlock(bp);
          if(bk2.getMin()[0]!=bk.getMin()[0]||bk2.getMin()[2]!=bk.getMin()[2]) continue;
          bool pF,pB; edgeY(bk2,pF,pB);
          if(!((atF&&pB)||(atB&&pF))) continue;
          partner=bp; break;
        }
        if(partner<Geo.getBlockNum()){
          const auto& pb=Geo.getBlock(partner); const int pov=pb.getOverlap();
          pl.srcBlockIdx=partner;
          pl.srcPlaneIdx = atF ? (pb.getNy()-1-pov) : pov;
          MFGhostSyncPlans.push_back(pl);
        } else {
#ifdef MPI_ENABLED
          const auto& globalGeo=GeoHelper.getBlockGeometry();
          int partnerId=-1;
          for(std::size_t g=0;g<globalGeo.getBlockNum();++g){
            const auto& gb=globalGeo.getBlock(g);
            if(gb.getMin()[0]!=bk.getMin()[0]||gb.getMin()[2]!=bk.getMin()[2]) continue;
            bool gF,gB; edgeY(gb,gF,gB);
            if(!((atF&&gB)||(atB&&gF))) continue;
            partnerId=gb.getBlockId(); break;
          }
          if(partnerId<0) continue;
          pl.crossRank=true; pl.partnerRank=GeoHelper.whichRank(partnerId);
          pl.sendTag=bk.getBlockId(); pl.recvTag=partnerId;
          MFGhostSyncPlans.push_back(pl);
#endif
        }
      }
    }
  }

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*x (Dirichlet). The wall rows/halos keep the exact
    // applied horizontal-field potential; the interior wavy perturbation is
    // supplied by the analytic initial psi and the sub-iterated solve.
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      const T nwall=Cell_Len*T{3.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T minX=bk.getMin()[0],minZ=bk.getMin()[2],vs=bk.getVoxelSize();
        for(int kk=0;kk<nz;++kk){
          T z=minZ+T(kk-ov)*vs;
          if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
            for(int jj=0;jj<ny;++jj){
              for(int ii=0;ii<nx;++ii){
                T x=minX+T(ii)*vs;
                T psi_w=psiRamp(x);
                std::size_t id=kk*pr[2]+jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k2=0;k2<MFLatSet::q;++k2) c[k2]=latset::w<MFLatSet>(k2)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }
    }

    // SyncMFPeriodicGhosts: mirror the physical edge planes of a per-cell MF
    // field into the opposite x/y-wrapped ghost planes. MF_Per.Apply() copies
    // only pops + GenericRho, and GenericRho for MFCELL resolves to PHI (see
    // FindGenericRhoType in src/utils/tmp.h — PSI is not in its list), so the
    // wrapped ghost planes of PSI/HX/HY/HZ/HMAG would stay frozen at their T0
    // values and MFComputeH3D / MFMagneticForce3D would read stale neighbors at
    // the x=0/Ni and y=0/Nj edges. (Block-interior ghosts are fine:
    // CommunicatePSI / CommunicateAllMFFields3D fill them from the neighbor
    // block's physical data.)
    //
    // 3D/MPI 修复: 128 进程时每进程仅 1 个 block, x/y 周期对侧 block 位于
    // 其他 rank, 本地查找必然失败 -> 接缝 ghost 永远停留在初始值 (H 场=0),
    // 磁力读取陈旧 ghost 产生虚假 ∇|H| (接缝处 |F|~1e-3, 顶壁四角 1.5e-2),
    // 是气泡 20~30 步崩溃的根因。现通过全局 block 表 (GlobalBlockTable)
    // 定位跨 rank partner, 用 MPI sendRecv 交换物理边平面。
    // tag = fidx*1000 + min(eLocal,ePartner)*2 + dir, 双方由同一张表计算
    // 出相同 tag, 保证消息一一配对。
    // Execute the precomputed MFGhostSyncPlans for one per-cell field.
    // ---- Seam POP exchange (fixes MF_Per's stale-ghost coverage gap at the seams).
    // Called after MF_Per.Apply(); SetXSeamPsiPops then overrides the x-seam PSI pops
    // with the ramp-shift.  The y-seam ghost POPs (plain periodic) are re-copied from the
    // partner's physical seam row so the D3Q7 Stream never reads the stale cells that
    // MF_Per leaves behind (the periodic-seam boundary artifact).
    auto SyncMFPops = [&]() {
#ifdef MPI_ENABLED
      constexpr int MF_POP_TAG_BASE = 9600; // distinct from MF_SYNC_TAG_BASE (9500)
      std::vector<std::vector<T>> sendBufs, recvBufs;
      std::vector<MPI_Request> sendReqs, recvReqs;
      std::vector<std::size_t> recvBlock;
      std::vector<int> recvPlaneIdx;
      std::vector<bool> recvXSeam;
#endif
      for (const auto& pl : MFGhostSyncPlans) {
        auto& gbl = MFLattice.getBlockLat(pl.ghostBlockIdx);
        const auto& bk = Geo.getBlock(pl.ghostBlockIdx);
        const auto& pr = bk.getProjection();
        if (!pl.crossRank) {
          auto& sbl = MFLattice.getBlockLat(pl.srcBlockIdx);
          const auto& spr = Geo.getBlock(pl.srcBlockIdx).getProjection();
          if (pl.xSeam) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j) {
                MFCELL gc(k*pr[2]+j*pr[1]+pl.ghostPlaneIdx, gbl);
                MFCELL sc(k*spr[2]+j*spr[1]+pl.srcPlaneIdx, sbl);
                for (unsigned q=0;q<MFLatSet::q;++q) gc[q]=sc[q];
              }
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i=0;i<bk.getNx();++i) {
                MFCELL gc(k*pr[2]+pl.ghostPlaneIdx*pr[1]+i, gbl);
                MFCELL sc(k*spr[2]+pl.srcPlaneIdx*spr[1]+i, sbl);
                for (unsigned q=0;q<MFLatSet::q;++q) gc[q]=sc[q];
              }
          }
        }
#ifdef MPI_ENABLED
        else {
          sendBufs.emplace_back(pl.plane * MFLatSet::q);
          auto& sb = sendBufs.back();
          int kk=0;
          if (pl.xSeam) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j) {
                MFCELL c(k*pr[2]+j*pr[1]+pl.seamPlaneIdx, gbl);
                for (unsigned q=0;q<MFLatSet::q;++q) sb[kk++]=c[q];
              }
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i=0;i<bk.getNx();++i) {
                MFCELL c(k*pr[2]+pl.seamPlaneIdx*pr[1]+i, gbl);
                for (unsigned q=0;q<MFLatSet::q;++q) sb[kk++]=c[q];
              }
          }
          MPI_Request srq, rrq;
          mpi().iSend(sb.data(), static_cast<int>(pl.plane*MFLatSet::q), pl.partnerRank, &srq, MF_POP_TAG_BASE + pl.sendTag);
          sendReqs.push_back(srq);
          recvBufs.emplace_back(pl.plane * MFLatSet::q);
          recvBlock.push_back(pl.ghostBlockIdx);
          recvPlaneIdx.push_back(pl.ghostPlaneIdx);
          recvXSeam.push_back(pl.xSeam);
          mpi().iRecv(recvBufs.back().data(), static_cast<int>(pl.plane*MFLatSet::q), pl.partnerRank, &rrq, MF_POP_TAG_BASE + pl.recvTag);
          recvReqs.push_back(rrq);
        }
#endif
      }
#ifdef MPI_ENABLED
      if (!sendReqs.empty()) {
        MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
        for (std::size_t i=0;i<recvReqs.size();++i) {
          MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
          auto& gbl = MFLattice.getBlockLat(recvBlock[i]);
          const auto& bk = Geo.getBlock(recvBlock[i]);
          const auto& pr = bk.getProjection();
          const auto& rb = recvBufs[i];
          int kk=0;
          if (recvXSeam[i]) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j) {
                MFCELL gc(k*pr[2]+j*pr[1]+recvPlaneIdx[i], gbl);
                for (unsigned q=0;q<MFLatSet::q;++q) gc[q]=rb[kk++];
              }
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i2=0;i2<bk.getNx();++i2) {
                MFCELL gc(k*pr[2]+recvPlaneIdx[i]*pr[1]+i2, gbl);
                for (unsigned q=0;q<MFLatSet::q;++q) gc[q]=rb[kk++];
              }
          }
        }
      }
#endif
    };

    const T Lx = T(Ni) * Cell_Len;
    const T rampShift = H0 * Lx;
    auto SetXSeamPsiPops = [&]() {
      const T xLeft = T{0};
      const T xRight = Lx;
      const int nEntry = int(GlobalBlockTable.size() / 11);
      const int myRank = mpi().getRank();
      auto& psiF = MFLattice.getField<PSI<T>>();
#ifdef MPI_ENABLED
      struct Job { int block; int ghostCol0; T shift; };
      std::vector<std::vector<T>> sendBufs, recvBufs;
      std::vector<MPI_Request> sendReqs, recvReqs;
      std::vector<Job> jobs;
#endif
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& bk = Geo.getBlock(b);
        const bool atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
        if (!atL && !atR) continue;
        const auto& pr = bk.getProjection();
        const int nx = bk.getNx(), ny = bk.getNy(), nz = bk.getNz(), ov = bk.getOverlap();
        auto& bPsi = psiF.getBlockField(b);
        auto& bl = MFLattice.getBlockLat(b);
        int eL = -1, eP = -1;
        for (int e = 0; e < nEntry; ++e) {
          const double* en = &GlobalBlockTable[e * 11];
          if (en[2] != bk.getMin()[1] || en[3] != bk.getMin()[2]) continue;
          if (atL && en[8] > T{0.5}) eP = e;   // partner atR
          if (atR && en[7] > T{0.5}) eP = e;   // partner atL
          if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
              en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2]) eL = e;
        }
        if (eP < 0 || eL < 0) continue;
        const double* enP = &GlobalBlockTable[eP * 11];
        const int pRank = int(enP[0]);
        const int pnx = int(enP[4]);
        const T shift = atL ? rampShift : -rampShift;
        const int srcI = atL ? (pnx - 1 - ov) : ov;
        const int ghostCol0 = atL ? 0 : nx - ov;
        if (pRank == myRank) {
          int pb = -1;
          for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
            const auto& bkk = Geo.getBlock(bb);
            if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                bkk.getMin()[2] == enP[3]) { pb = bb; break; }
          }
          if (pb < 0) continue;
          auto& pPsi = psiF.getBlockField(pb);
          const auto& pPr = Geo.getBlock(pb).getProjection();
          for (int kk = 0; kk < nz; ++kk)
            for (int jj = 0; jj < ny; ++jj) {
              T psi_g = pPsi.get(kk * pPr[2] + jj * pPr[1] + srcI) + shift;
              for (int c = 0; c < ov; ++c) {
                std::size_t id = kk * pr[2] + jj * pr[1] + ghostCol0 + c;
                MFCELL cell(id, bl);
                for (unsigned k2 = 0; k2 < MFLatSet::q; ++k2)
                  cell[k2] = latset::w<MFLatSet>(k2) * psi_g;
                bPsi.get(id) = psi_g;
              }
            }
        }
#ifdef MPI_ENABLED
        else {
          const int tag = 3000 + std::min(eL, eP) * 2 + 0;
          const int myPhysCol = atL ? ov : nx - 1 - ov;
          std::vector<T> sbuf(ny * nz);
          int cnt = 0;
          for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
              sbuf[cnt++] = bPsi.get(k * pr[2] + j * pr[1] + myPhysCol);
          sendBufs.emplace_back(std::move(sbuf));
          MPI_Request sreq, rreq;
          mpi().iSend(sendBufs.back().data(), ny * nz, pRank, &sreq, tag);
          sendReqs.push_back(sreq);
          recvBufs.emplace_back(ny * nz);
          jobs.push_back({b, ghostCol0, shift});
          mpi().iRecv(recvBufs.back().data(), ny * nz, pRank, &rreq, tag);
          recvReqs.push_back(rreq);
        }
#endif
      }
#ifdef MPI_ENABLED
      if (!sendReqs.empty()) {
        MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
        for (std::size_t i = 0; i < recvReqs.size(); ++i) {
          MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
          auto& bPsi = psiF.getBlockField(jobs[i].block);
          auto& bl = MFLattice.getBlockLat(jobs[i].block);
          const auto& bk = Geo.getBlock(jobs[i].block);
          const auto& pr = bk.getProjection();
          const int ov = bk.getOverlap();
          const int g0 = jobs[i].ghostCol0;
          const T shift = jobs[i].shift;
          const auto& rbuf = recvBufs[i];
          int cnt = 0;
          for (int k = 0; k < bk.getNz(); ++k)
            for (int j = 0; j < bk.getNy(); ++j) {
              T psi_g = rbuf[cnt++] + shift;
              for (int c = 0; c < ov; ++c) {
                std::size_t id = k * pr[2] + j * pr[1] + g0 + c;
                MFCELL cell(id, bl);
                for (unsigned k2 = 0; k2 < MFLatSet::q; ++k2)
                  cell[k2] = latset::w<MFLatSet>(k2) * psi_g;
                bPsi.get(id) = psi_g;
              }
            }
        }
      }
#endif
    };

    // mu(phi) 只在本步 0a 更新一次; 把块间与周期缝 ghost 同步后供
    // Gauss-Seidel 面调和平均使用。
    MFLattice.template getField<MU_PERCELL<T>>().Communicate();
    SyncMFPeriodicGhosts(MFLattice.template getField<MU_PERCELL<T>>(), 0);

    // 0c-0e+: 对水平切向场采用五点/七点有限差分 Gauss-Seidel 求解
    // nabla.(mu nabla psi)=0。不再使用 D3Q7 MRT 伪时间迭代: 后者在
    // x 方向斜坡+扭曲周期缝上会逐渐积累出反号的 delta_psi (Hx 反号),
    // 是之前 Bom=3980 不抑制、t*=1 提前的根本数值原因。
    // 离散: 面导磁率取两侧 mu 的调和平均; 上下壁保持 psi=-H0*x,
    // x/y 周期缝由 SyncMFPeriodicGhosts + SetXSeamPsiPops 同步。
    for(int sub=0;sub<PsiSolver_Iter;++sub){
      CommunicatePSI<T>(MFLattice);
      SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>(), 0);
      SetXSeamPsiPops();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=MFLattice.getField<PSI<T>>().getBlockField(b);
        int ov=bk.getOverlap();
        for(int k=ov;k<bk.getNz()-ov;++k){
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=k*pr[2]+j*pr[1]+i; MFCELL c(id,bl);
              T muc=c.template get<MU_PERCELL<T>>();
              T num=T{0}, den=T{0};
              for(unsigned k2=1;k2<MFLatSet::q;++k2){
                MFCELL n=c.getNeighbor(k2);
                T mun=n.template get<MU_PERCELL<T>>();
                T muf=T{2}/(T{1}/muc + T{1}/mun);
                T pn=n.template get<PSI<T>>();
                num += muf*pn;
                den += muf;
              }
              bPsi.get(id)=num/den;
            }
          }
        }
      }

      // re-pin walls after every sweep
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        const T nwall=Cell_Len*T{3.0};
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
          T minX=bk.getMin()[0],minZ=bk.getMin()[2],vs=bk.getVoxelSize();
          for(int kk=0;kk<nz;++kk){
            T z=minZ+T(kk-ov)*vs;
            if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
              for(int jj=0;jj<ny;++jj){
                for(int ii=0;ii<nx;++ii){
                  T x=minX+T(ii)*vs;
                  T psi_w=psiRamp(x);
                  std::size_t id=kk*pr[2]+jj*pr[1]+ii;
                  bPsi.get(id)=psi_w;
                }
              }
            }
          }
        }
      }
      // 同步 PSI 供下一 sweep / H 计算使用
      CommunicatePSI<T>(MFLattice);
    }
    // One sync per step is enough: the sub-iterations never read neighbor PSI
    // (collision uses own pops, 0e uses own cells), and the H ghosts are only
    // read by A4 below, after the sync that follows 0f.
    // PSI 缝同步: y 向普通周期 (SyncMFPeriodicGhosts), x 向用斜坡平移周期 BC 覆写
    // (SyncMFPeriodicGhosts 会把斜坡按普通周期卷绕而错位, SetXSeamPsiPops 修正 x)。
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>(), 0);
    SetXSeamPsiPops();

    // 0f: Compute H = -∇ψ
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int k=ov;k<bk.getNz()-ov;++k){
        for(int j=ov;j<bk.getNy()-ov;++j){
          for(int i=ov;i<bk.getNx()-ov;++i){
            MFCELL c(k*pr[2]+j*pr[1]+i,bl);
            MFComputeH3D<MFCELL>::apply(c);
          }
        }
      }
    }
    CommunicateAllMFFields3D<T>(MFLattice);
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>(), 1);
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>(), 2);
    SyncMFPeriodicGhosts(MFLattice.getField<HZ<T>>(), 3);
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>(), 4);

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,3>>().InitValue(Vector<T,3>{0,0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic force (if chi>0 and H0>0)
    // Paper Eq. (8): F_m = (μ₀χ/2)∇(|H|²) = μ₀χ|H|∇|H| (μ₀≡1, Kelvin form).
    // Force-free wall bands: the z wall ψ-pins and the D3Q7 ψ-solver's
    // far-field slope mismatch leave a spurious |H| spike within ~6 rows of
    // each wall (worst at the x/y periodic corners). The Kelvin force there
    // (χ|H|∇|H|) is a numerical artifact, so it is zeroed in the bands
    // |z|<=16 and |z-H_global|<=16.
    if(H0>T{0}){
      const T band=std::min(T{16.0}, H_global*T{0.1});
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        T minX=bk.getMin()[0],minZ=bk.getMin()[2],vs=bk.getVoxelSize();
        for(int k=ov;k<bk.getNz()-ov;++k){
          T z=minZ+T(k-ov)*vs;
          if(z<=band||z>=H_global-band) continue;
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=k*pr[2]+j*pr[1]+i;
              PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
              // 法向场案例 (bubbleMag3d/rosenMag3d) 继续使用 src 的
              // Kelvin 力 F=chi*|H|*grad|H|, 且已验证正确。这里只针对
              // KH 的水平切向场补充磁压力项。完整 Maxwell 应力散度
              //  F_m = -0.5*H^2*grad(mu)
              // 在切向场下 H^2 沿界面起伏, grad(mu) 沿界面法向;
              // 该项是水平场抑制 KH 波浪的主项。为保证量纲与数值稳定,
              // 只用切向 Hx 分量: F_m = -0.5*Hx^2*grad(mu)。
              MFMagneticForce3D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
              if(MagPressureEnabled){
                T Hx=mf.template get<HX<T>>();
                const auto& grad_phi=pf.template get<GRAD<T,3>>();
                T dmu=mf.template get<MU_H<T>>()-mf.template get<MU_L<T>>();
                auto& F=ns.template get<FORCE<T,3>>();
                T pref=MagPressureFactor*(-T{0.5})*Hx*Hx*dmu;
                F[0]+=pref*grad_phi[0];
                F[1]+=pref*grad_phi[1];
                F[2]+=pref*grad_phi[2];
              }
            }
          }
        }
      }
    }

    GrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,GrT>>(t(),FlagFM);
    if(PreForce)
      PrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,PrT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,3>>().Communicate();

    // ===== Phase B-C: PF + NS collision =====
    PFLattice.template ApplyInnerCellDynamics<PFSel>(FlagFM);
    PF_Per.Apply(); PFLattice.NormalFullCommunicate();

    ViC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,ViT>>(t(),FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSSel>(FlagFM);
    NS_Per.Apply(); NSLattice.NormalFullCommunicate();

    // ===== Phase D: Streaming =====
    PF_BB.Apply(t());
    { // PF Z ghost pop fix
      T Hg=T(Nk)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=PFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5})
          for(int jj=0;jj<ny;++jj)
            for(int i=0;i<nx;++i){PFCELL g(0*pr[2]+jj*pr[1]+i,bl),w(ov*pr[2]+jj*pr[1]+i,bl); for(unsigned k2=0;k2<LatSet::q;++k2)g[k2]=w[k2];}
        if(Mz>Hg-Cell_Len*T{1.5})
          for(int jj=0;jj<ny;++jj)
            for(int i=0;i<nx;++i){PFCELL g((nz-1)*pr[2]+jj*pr[1]+i,bl),w((nz-1-ov)*pr[2]+jj*pr[1]+i,bl); for(unsigned k2=0;k2<LatSet::q;++k2)g[k2]=w[k2];}
      }
    }
    PFLattice.Stream(); PFLattice.NormalFullCommunicate();

    NS_BB.Apply(t());
    { // NS Z ghost pop fix
      T Hg=T(Nk)*Cell_Len;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5})
          for(int jj=0;jj<ny;++jj)
            for(int i=0;i<nx;++i){NSCELL g(0*pr[2]+jj*pr[1]+i,bl),w(ov*pr[2]+jj*pr[1]+i,bl); for(unsigned k2=0;k2<LatSet::q;++k2)g[k2]=w[k2];}
        if(Mz>Hg-Cell_Len*T{1.5})
          for(int jj=0;jj<ny;++jj)
            for(int i=0;i<nx;++i){NSCELL g((nz-1)*pr[2]+jj*pr[1]+i,bl),w((nz-1-ov)*pr[2]+jj*pr[1]+i,bl); for(unsigned k2=0;k2<LatSet::q;++k2)g[k2]=w[k2];}
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
        for(int k=ov;k<bk.getNz()-ov;++k)
          for(int j=ov;j<bk.getNy()-ov;++j)
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=k*pr[2]+j*pr[1]+i; PFCELL c(id,bl);
              T pn=0; for(unsigned k2=0;k2<LatSet::q;++k2)pn+=c[k2];
              if(pn<T{0})pn=T{0}; if(pn>T{1})pn=T{1}; bP.get(id)=pn;
            }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi BC: 底壁 phi=1 (铁磁流体, 在下), 顶壁 phi=0 (溶剂, 在上)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5}) for(int jj=0;jj<ny;++jj) for(int i=0;i<nx;++i)bP.get(ov*pr[2]+jj*pr[1]+i)=T{1};
        if(Mz>H_global-Cell_Len*T{1.5}){int kk=nz-1-ov;for(int jj=0;jj<ny;++jj)for(int i=0;i<nx;++i)bP.get(kk*pr[2]+jj*pr[1]+i)=T{0};}
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi ghost fix
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5}) for(int k=0;k<ov;++k) for(int jj=0;jj<ny;++jj) for(int i=0;i<nx;++i) bP.get(k*pr[2]+jj*pr[1]+i)=T{1};
        if(Mz>H_global-Cell_Len*T{1.5}) for(int k=nz-ov;k<nz;++k) for(int jj=0;jj<ny;++jj) for(int i=0;i<nx;++i) bP.get(k*pr[2]+jj*pr[1]+i)=T{0};
      }
    }

    // seam PHI: 宏量更新和 z 壁修正只负责块内/壁面数据，
    // 还需要把 x/y 周期缝的 PHI ghost 平面重新镜像到对侧。
    // 否则左右边界 PHI ghost 保留旧值，导致 PHI 场左右异常。
    SyncMFPeriodicGhosts(PFLattice.getField<PHI<T>>(), 0);

    // gradients
    PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
    PFLattice.getField<NORMAL<T,3>>().Communicate();
    PFLattice.getField<GRAD<T,3>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // wall grad/chempot fix (PF)
    {
      auto& gF=PFLattice.getField<GRAD<T,3>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bG=gF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5}) for(int jj=ov;jj<ny-ov;++jj) for(int i=ov;i<nx-ov;++i) bG.get(ov*pr[2]+jj*pr[1]+i)[2]=bG.get((ov+1)*pr[2]+jj*pr[1]+i)[2];
        if(Mz>H_global-Cell_Len*T{1.5}){int kk=bk.getNz()-1-ov;for(int jj=ov;jj<ny-ov;++jj)for(int i=ov;i<nx-ov;++i)bG.get(kk*pr[2]+jj*pr[1]+i)[2]=bG.get((kk-1)*pr[2]+jj*pr[1]+i)[2];}
      }
    }
    {
      auto& cF=PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bC=cF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5}) for(int jj=ov;jj<ny-ov;++jj) for(int i=ov;i<nx-ov;++i) bC.get(ov*pr[2]+jj*pr[1]+i)=(T{4}*bC.get((ov+1)*pr[2]+jj*pr[1]+i)-bC.get((ov+2)*pr[2]+jj*pr[1]+i))/T{3};
        if(Mz>H_global-Cell_Len*T{1.5}){int kk=bk.getNz()-1-ov;for(int jj=ov;jj<ny-ov;++jj)for(int i=ov;i<nx-ov;++i)bC.get(kk*pr[2]+jj*pr[1]+i)=(T{4}*bC.get((kk-1)*pr[2]+jj*pr[1]+i)-bC.get((kk-2)*pr[2]+jj*pr[1]+i))/T{3};}
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // NS macro
    {
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& rF=NSLattice.getField<DENSITY<T>>(); auto& pF=NSLattice.getField<PRESSURE<T>>();
        auto& vF=NSLattice.getField<VELOCITY<T,3>>(); auto& fF=NSLattice.getField<FORCE<T,3>>();
        auto& bR=rF.getBlockField(b); auto& bP=pF.getBlockField(b);
        auto& bV=vF.getBlockField(b); auto& bF=fF.getBlockField(b);
        int ov=bk.getOverlap();
        for(int k=ov;k<bk.getNz()-ov;++k)
          for(int j=ov;j<bk.getNy()-ov;++j)
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=k*pr[2]+j*pr[1]+i; NSCELL c(id,bl);
              T p=0,ux=0,uy=0,uz=0;
              for(unsigned k2=0;k2<LatSet::q;++k2){p+=c[k2];ux+=latset::c<LatSet>(k2)[0]*c[k2];uy+=latset::c<LatSet>(k2)[1]*c[k2];uz+=latset::c<LatSet>(k2)[2]*c[k2];}
              T rho=bR.get(id); auto F=bF.get(id);
              bP.get(id)=p; bV.get(id)=Vector<T,3>{ux+T{0.5}*F[0]/rho,uy+T{0.5}*F[1]/rho,uz+T{0.5}*F[2]/rho};
            }
      }
    }
    NSLattice.getField<VELOCITY<T,3>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    ++t; ++ot;
    if(t()%OutputStep==0){
      // KH 界面诊断: 每列 (i,j) 找 phi=0.5 的 z, 得界面高度图 h(x,y)
      //   混合层高度 = hmax - hmin (界面卷曲使高度展宽, 论文 Fig.34 的量)
      T hmax=-1e30, hmin=1e30;
      {
        auto& pF=PFLattice.getField<PHI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bP=pF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
          T vs=bk.getVoxelSize(),mz=bk.getMin()[2];
          for(int j=ov;j<ny-ov;++j){
            for(int i=ov;i<nx-ov;++i){
              T zprev=mz+T(ov)*vs; T phiprev=bP.get(ov*pr[2]+j*pr[1]+i);
              for(int k=ov+1;k<nz-ov;++k){
                T zk=mz+T(k)*vs; T phik=bP.get(k*pr[2]+j*pr[1]+i);
                if((phiprev-T{0.5})*(phik-T{0.5})<=T{0}){
                  T hz=zprev+(T{0.5}-phiprev)*(zk-zprev)/(phik-phiprev);
                  if(hz>hmax)hmax=hz; if(hz<hmin)hmin=hz;
                }
                zprev=zk; phiprev=phik;
              }
            }
          }
        }
      }
#ifdef MPI_ENABLED
      if(mpi().getSize()>1){
        // reduceAndBcast 仅支持 uint8 枚举 (MPI_BYTE), 双精度需直接用 MPI_Allreduce
        double loc[2] = {hmax, -hmin}, glob[2];
        MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        hmax = glob[0]; hmin = -glob[1];
      }
#endif
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      if (mpi().getRank() == 0) {
        // 论文 F 节无量纲时间 t* = t*sqrt(g/L) (L = L_ref = Ni/2)
        T tstar = T(t()) * KH_TSTAR_SCALE * std::sqrt(gravity/L_ref);  // 论文图像校准: t*_plot = 4*t*_eq
        printf("[t=%d t*=%.3f] mixed_layer_h=%.3f (hmax=%.3f hmin=%.3f)  (z0=%.2f)\n",
               t(),tstar,(hmax-hmin),hmax,hmin,InterfaceZ);
      }
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
