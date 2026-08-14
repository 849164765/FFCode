#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D3Q19<T>;
using MFLatSet = D3Q7<T>;
using namespace mfield;

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

int MaxStep, OutputStep;
std::string work_dir;

void readParam(int argc, char* argv[]) {
  std::string iniName = "khMag3d.ini";
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
  eta_h = r.getValue<T>("Two_Phase","eta_h");
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

  // 由无量纲数推导格子参数 (论文 F 节):
  //   rho_h = rho_l*ratio_hl (铁磁流体在 φ=1, 略轻: rho_h<rho_l),
  //   Re = rho_l*U0*L/eta_l,  We = rho_l*L*U0^2/sigma,  Fr = U0^2/(gL),
  //   Bo_m = mu0*H0^2*L/(2*sigma) (mu0=1)
  // L_dom = 域宽 (x 方向, = 2L), L_ref = 参考长度 L = Ni/2 (论文 F 节 L=256 的 4x 缩小).
  // 无量纲数/界面波长必须都用 L_ref, 与 ini 几何 (λ=L, A=0.1L, z0=2L) 一致.
  const T pi = std::acos(T{-1});
  L_dom = T(Ni)*Cell_Len;
  L_ref = T(Ni)*Cell_Len/T{2};
  T ratio_hl = r.getValue<T>("Two_Phase","rho_hl_ratio", T{0.99});
  rho_h = rho_l*ratio_hl;
  eta_l = rho_l*U0*L_ref/Re;
  sigma = rho_l*L_ref*U0*U0/We;
  gravity = U0*U0/(Fr*L_ref);
  // 显式覆盖 (可选): 若 ini 直接给出这些格子值则优先
  rho_h   = r.getValue<T>("Two_Phase","rho_h", rho_h);
  eta_l   = r.getValue<T>("Two_Phase","eta_l", eta_l);
  sigma   = r.getValue<T>("Two_Phase","sigma", sigma);
  gravity = r.getValue<T>("Two_Phase","Gravity", gravity);
  DeltaRho=rho_h-rho_l;
  H0 = r.getValue<T>("Magnetic_Field","H0", T{-1});
  if(H0 < T{0}) H0 = std::sqrt(T{2}*sigma*Bo_m/L_ref);   // mu0=1

  Beta=T{12.0}*sigma/Interface_Width;
  Kappa=T{3.0}*Interface_Width*sigma*T{0.5};
  Tau_phi=T{3.0}*Mobility+T{0.5}; Omega_phi=T{1.0}/Tau_phi;
  Tau_ns=T{0.5}+eta_h/rho_h/LatSet::cs2;

  MPI_RANK(0){
    printf("---- Kelvin-Helmholtz Instability in Ferrofluid (3D) ----\n");
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
    printf("------------------------------------------\n");
    if(H0<=T{0}) printf("[Note] H0=0 (Bo_m=0): 无磁场对照\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Kelvin-Helmholtz Instability in 3D..."));
  readParam(argc, argv);

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni,T(0.01),T(1.0));
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
  //   h(x) = z0 + A*cos(2πx/L)   (论文 F 节, 扰动只沿剪切方向 x, y 不变;
  //   域宽 2L 故含 2 个周期; 3D 为沿 y 挤出的卷状结构)
  T z0=InterfaceZ*Cell_Len, W_phys=Interface_Width*Cell_Len;
  T lam=L_dom/PerturbPeriods;
  T kx=T{2}*std::acos(T{-1})/lam;
  auto interfaceHeight=[&](T x, T y)->T{
    return z0 + PerturbAmp*Cell_Len*std::cos(kx*x);
  };
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1],mz=bk.getMin()[2];
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      T z=mz+T(k)*vs;
      for(int j=ov;j<bk.getNy()-ov;++j){
        T y=my+T(j)*vs;
        for(int i=ov;i<bk.getNx()-ov;++i){
          T x=mx+T(i)*vs;
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
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1],mz=bk.getMin()[2];
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      T z=mz+T(k)*vs;
      for(int j=ov;j<bk.getNy()-ov;++j){
        T y=my+T(j)*vs;
        for(int i=ov;i<bk.getNx()-ov;++i){
          T x=mx+T(i)*vs;
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
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0];
      for(int k=0;k<bk.getNz();++k){
        for(int j=0;j<bk.getNy();++j){
          for(int i=0;i<bk.getNx();++i){
            T x=mx+T(i)*vs, psi=psiRamp(x);
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
    vtmo::vtmWriter<T,3> MW("khMag3d",Geo);
  MW.addWriterSet(PW);

  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();
  NS_Per.Apply(); PF_Per.Apply();

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

    // 0b: Set wall psi (layered: 轻流体 -H0*z, 重流体 -H0*(mu_l/mu_h)*z+C) and
    // wall pops = feq(psi_bc). NOTE: lattice row kk holds physical z =
    // minZ + (kk-ov)*vs (ov halo rows below the block). Pinning by absolute z
    // (halo + wall rows) gives the exact layered profile.
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
    auto SyncMFPeriodicGhosts = [&](auto& field) {
#ifdef MPI_ENABLED
      constexpr int MF_SYNC_TAG_BASE = 9500; // distinct from PERIODIC_TAG_BASE (9000)
      std::vector<std::vector<T>> sendBufs, recvBufs;
      std::vector<MPI_Request> sendReqs, recvReqs;
      std::vector<std::size_t> recvBlock;
      std::vector<int> recvPlaneIdx;
      std::vector<bool> recvXSeam;
#endif
      for (const auto& pl : MFGhostSyncPlans) {
        auto& fg = field.getBlockField(pl.ghostBlockIdx);
        const auto& bk = Geo.getBlock(pl.ghostBlockIdx);
        const auto& pr = bk.getProjection();
        if (!pl.crossRank) {
          auto& fs = field.getBlockField(pl.srcBlockIdx);
          const auto& spr = Geo.getBlock(pl.srcBlockIdx).getProjection();
          if (pl.xSeam) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j)
                fg.get(k*pr[2]+j*pr[1]+pl.ghostPlaneIdx) = fs.get(k*spr[2]+j*spr[1]+pl.srcPlaneIdx);
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i=0;i<bk.getNx();++i)
                fg.get(k*pr[2]+pl.ghostPlaneIdx*pr[1]+i) = fs.get(k*spr[2]+pl.srcPlaneIdx*spr[1]+i);
          }
        }
#ifdef MPI_ENABLED
        else {
          sendBufs.emplace_back(pl.plane);
          auto& sb = sendBufs.back();
          int kk=0;
          if (pl.xSeam) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j)
                sb[kk++] = fg.get(k*pr[2]+j*pr[1]+pl.seamPlaneIdx);
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i=0;i<bk.getNx();++i)
                sb[kk++] = fg.get(k*pr[2]+pl.seamPlaneIdx*pr[1]+i);
          }
          MPI_Request srq, rrq;
          mpi().iSend(sb.data(), pl.plane, pl.partnerRank, &srq, MF_SYNC_TAG_BASE + pl.sendTag);
          sendReqs.push_back(srq);
          recvBufs.emplace_back(pl.plane);
          recvBlock.push_back(pl.ghostBlockIdx);
          recvPlaneIdx.push_back(pl.ghostPlaneIdx);
          recvXSeam.push_back(pl.xSeam);
          mpi().iRecv(recvBufs.back().data(), pl.plane, pl.partnerRank, &rrq, MF_SYNC_TAG_BASE + pl.recvTag);
          recvReqs.push_back(rrq);
        }
#endif
      }
#ifdef MPI_ENABLED
      if (!sendReqs.empty()) {
        MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
        for (std::size_t i=0;i<recvReqs.size();++i) {
          MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
          auto& fg = field.getBlockField(recvBlock[i]);
          const auto& bk = Geo.getBlock(recvBlock[i]);
          const auto& pr = bk.getProjection();
          const auto& rb = recvBufs[i];
          int kk=0;
          if (recvXSeam[i]) {
            for (int k=0;k<bk.getNz();++k)
              for (int j=0;j<bk.getNy();++j)
                fg.get(k*pr[2]+j*pr[1]+recvPlaneIdx[i]) = rb[kk++];
          } else {
            for (int k=0;k<bk.getNz();++k)
              for (int i2=0;i2<bk.getNx();++i2)
                fg.get(k*pr[2]+recvPlaneIdx[i]*pr[1]+i2) = rb[kk++];
          }
        }
      }
#endif
    };

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
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& bk = Geo.getBlock(b);
        const bool atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
        if (!atL && !atR) continue;
        const auto& pr = bk.getProjection();
        const int nx = bk.getNx(), ny = bk.getNy(), nz = bk.getNz();
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
        const int srcI = atL ? (pnx - 2) : 1;
        const int dstI = atL ? 0 : (nx - 1);
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
              std::size_t id = kk * pr[2] + jj * pr[1] + dstI;
              MFCELL c(id, bl);
              for (unsigned k2 = 0; k2 < MFLatSet::q; ++k2)
                c[k2] = latset::w<MFLatSet>(k2) * psi_g;
              bPsi.get(id) = psi_g;
            }
        } else {
#ifdef MPI_ENABLED
          const int tag = 3000 + std::min(eL, eP) * 2 + 0;
          const int sendI = atL ? 1 : (nx - 2);
          std::vector<T> sbuf(ny * nz), rbuf(ny * nz);
          int kk = 0;
          for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
              sbuf[kk++] = bPsi.get(k * pr[2] + j * pr[1] + sendI);
          mpi().sendRecv(sbuf.data(), rbuf.data(), ny * nz, pRank, pRank, tag);
          kk = 0;
          for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j) {
              T psi_g = rbuf[kk++] + shift;
              std::size_t id = k * pr[2] + j * pr[1] + dstI;
              MFCELL c(id, bl);
              for (unsigned k2 = 0; k2 < MFLatSet::q; ++k2)
                c[k2] = latset::w<MFLatSet>(k2) * psi_g;
              bPsi.get(id) = psi_g;
            }
#endif
        }
      }
    };

    // 0c-0e+: psi-solver sub-iterations (PsiSolver_Iter repeats of
    // collision+stream+macro). Each sub-iteration is one D3Q7 MRT diffusion
    // step with the boosted relaxation omega_psi = 1/(0.5+K*mu); ~100
    // iterations converge the solve to its (smooth) fixed point, removing the
    // grid-scale |H| oscillation that would facet the interface.
    for(int sub=0;sub<PsiSolver_Iter;++sub){
      // 0c: MF collision (direct iteration)
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int k=ov;k<bk.getNz()-ov;++k){
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              MFCELL c(k*pr[2]+j*pr[1]+i,bl);
              collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
            }
          }
        }
      }
      MF_Per.Apply();
      SyncMFPops();   // re-exchange the seam ghost POPs (SetXSeamPsiPops overrides x-seam)

      // 0d: Stream (先施加 x 缝斜坡平移周期 BC: 提供正确的入流分布)
      MFLattice.NormalFullCommunicate();
      SetXSeamPsiPops();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

      // 0e: PSI = Σg_i (interior only, walls keep psi_bc)
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
          for(int k=ov;k<bk.getNz()-ov;++k){
            for(int j=ov;j<bk.getNy()-ov;++j){
              for(int i=ov;i<bk.getNx()-ov;++i){
                std::size_t id=k*pr[2]+j*pr[1]+i; MFCELL c(id,bl);
                T psi=0; for(unsigned k2=0;k2<MFLatSet::q;++k2)psi+=c[k2];
                bPsi.get(id)=psi;
              }
            }
          }
        }
      }
      CommunicatePSI<T>(MFLattice);

      // 0e+: re-pin walls (0e overwrote them with the drifted Σg) so the
      // stored field is exactly -H0*z at the walls/halo before H is computed.
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
    }
    // One sync per step is enough: the sub-iterations never read neighbor PSI
    // (collision uses own pops, 0e uses own cells), and the H ghosts are only
    // read by A4 below, after the sync that follows 0f.
    // PSI 缝同步: y 向普通周期 (SyncMFPeriodicGhosts), x 向用斜坡平移周期 BC 覆写
    // (SyncMFPeriodicGhosts 会把斜坡按普通周期卷绕而错位, SetXSeamPsiPops 修正 x)。
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>());
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
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HZ<T>>());
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>());

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
      const T band=T{16.0};
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
              MFMagneticForce3D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
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
        T tstar = T(t()) * std::sqrt(gravity/L_ref);
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
