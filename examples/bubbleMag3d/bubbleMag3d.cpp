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

// bubble
T Bubble_Radius;
Vector<T, 3> Bubble_Center;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Eo, Re, U_g, Tau_ns;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bom, DeltaRho;

int PsiSolver_Iter; T PsiSolver_K;

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("bubbleMag3d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Nk = r.getValue<int>("Mesh","Nk");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Bubble_Radius = r.getValue<T>("Bubble","Radius");
  Bubble_Center[0]=r.getValue<T>("Bubble","CenterX");
  Bubble_Center[1]=r.getValue<T>("Bubble","CenterY");
  Bubble_Center[2]=r.getValue<T>("Bubble","CenterZ");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  Eo=r.getValue<T>("Two_Phase","Eo");
  Re=r.getValue<T>("Two_Phase","Re");
  U_g=r.getValue<T>("Two_Phase","U_g");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  Bom=r.getValue<T>("Magnetic_Field","Bom");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");

  eta_l=T(0.0568)/T(100.0); eta_h=T(0.0568);
  eta_l=r.getValue<T>("Two_Phase","eta_l",eta_l);
  eta_h=r.getValue<T>("Two_Phase","eta_h",eta_h);
  // Compute sigma, gravity from Eo and Re (Guo 2025 definitions)
  // Re = sqrt(|Gy|)*rho_h*D^(3/2)/eta_h
  // Eo = |Gy|*rho_h*D^2/sigma
  T D_bubble = T{2} * Bubble_Radius;
  T gy_sqrt = Re * eta_h / (rho_h * std::pow(D_bubble, T{1.5}));
  gravity = gy_sqrt * gy_sqrt;
  sigma = gravity * rho_h * D_bubble * D_bubble / Eo;
  // H0 from Bo_m: Bo_m = μ₀ * H₀² * D / (2*σ),  μ₀=1 in LBM units
  H0 = std::sqrt(T(2.0) * Bom * sigma / D_bubble);
  Beta=T(12.0)*sigma/Interface_Width;
  Kappa=T(3.0)*Interface_Width*sigma*T(0.5);
  Tau_phi=T(3.0)*Mobility+T(0.5); Omega_phi=T(1.0)/Tau_phi;
  Tau_ns=T(0.5)+eta_h/rho_h/LatSet::cs2;
  DeltaRho=rho_h-rho_l;

  T cs=std::sqrt(LatSet::cs2), Ma=U_g/cs;
  MPI_RANK(0){
    printf("---- Bubble Rising in Ferrofluid (3D) ----\n");
    printf("Mesh: %dx%dx%d  BlockCellLen=%d\n",Ni,Nj,Nk,BlockCellLen);
    printf("Bubble: R=%.1f center=(%.0f,%.0f,%.0f)\n",Bubble_Radius,Bubble_Center[0],Bubble_Center[1],Bubble_Center[2]);
    printf("rho: l=%.4f h=%.1f  eta: l=%.6f h=%.6f  sigma=%.2e g=%.2e\n",rho_l,rho_h,eta_l,eta_h,sigma,gravity);
    printf("Eo=%.0f Re=%.0f W=%.1f M=%.3f tau_phi=%.3f\n",Eo,Re,Interface_Width,Mobility,Tau_phi);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) H0=%.3f Bom=%.3f\n",chi_l,chi_h,mu_l,mu_h,H0,Bom);
    printf("PsiSolver: K=%.3f iter=%d (omega_mu9=%.3f omega_mu1=%.3f)\n",PsiSolver_K,PsiSolver_Iter, T{1}/(T{0.5}+PsiSolver_K*mu_h),T{1}/(T{0.5}+PsiSolver_K*mu_l));
    printf("U_g=%.3f Ma=%.3f\n",U_g,Ma);
    printf("------------------------------------------\n");
  }
  if(Ma>T(0.2)) MPI_RANK(0){ fprintf(stderr,"[Warn] Ma=%.3f > 0.2\n",Ma); }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Bubble Rising in Ferrofluid (3D)..."));
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
  // Z is vertical (H_global = Nk*Cell_Len), walls in Z, X & Y are periodic.
  AABB<T,3> domain({0,0,0},{T(Ni*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> left  ({T(-Cell_Len),0,0},{0,T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> right ({T(Ni*Cell_Len),0,0},{T((Ni+1)*Cell_Len),T(Nj*Cell_Len),T(Nk*Cell_Len)});
  AABB<T,3> front ({0,T(-Cell_Len),0},{T(Ni*Cell_Len),0,T(Nk*Cell_Len)});
  AABB<T,3> back  ({0,T(Nj*Cell_Len),0},{T(Ni*Cell_Len),T((Nj+1)*Cell_Len),T(Nk*Cell_Len)});
  BlockGeometryHelper3D<T> GeoHelper(Ni,Nj,Nk,domain,Cell_Len,BlockCellLen);
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
  // FIX(边界伪影, 同 rtMag3d/rosenMag3d): SetupBoundary 会把 x/y 周期侧面上
  // 的域边界格点误设为 Bounceback 无滑移壁。把四个侧面的首/末物理列/行重置为
  // BulkFlag, 仅保留 z 顶底壁为 Bounceback。
  {
    AABB<T,3> left_col ({T{0},T{0},Cell_Len},{Cell_Len,T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
    AABB<T,3> right_col({T((Ni-1)*Cell_Len),T{0},Cell_Len},{T(Ni*Cell_Len),T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
    AABB<T,3> front_col({T{0},T{0},Cell_Len},{T(Ni*Cell_Len),Cell_Len,T((Nk-1)*Cell_Len)});
    AABB<T,3> back_col ({T{0},T((Nj-1)*Cell_Len),Cell_Len},{T(Ni*Cell_Len),T(Nj*Cell_Len),T((Nk-1)*Cell_Len)});
    FlagFM.forEach(left_col, [&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
    FlagFM.forEach(right_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
    FlagFM.forEach(front_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
    FlagFM.forEach(back_col, [&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  }

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

  // -- init phi (3D tanh spherical bubble) --
  T R_phys=Bubble_Radius*Cell_Len, xc=Bubble_Center[0]*Cell_Len;
  T yc=Bubble_Center[1]*Cell_Len, zc=Bubble_Center[2]*Cell_Len;
  T W_phys=Interface_Width*Cell_Len;
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1],mz=bk.getMin()[2];
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k){
      T z=mz+T(k)*vs, dz=z-zc;
      for(int j=ov;j<bk.getNy()-ov;++j){
        T y=my+T(j)*vs, dy=y-yc;
        for(int i=ov;i<bk.getNx()-ov;++i){
          T x=mx+T(i)*vs, dx=x-xc;
          T dist=std::sqrt(dx*dx+dy*dy+dz*dz);
          T phi=T{0.5}+T{0.5}*std::tanh(T{2.0}*(dist-R_phys)/W_phys);
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

  // init NS pops (p=0, u=0)
  Vector<T,3> uz{0,0,0}; T pz=0;
  for(int b=0;b<Geo.getBlockNum();++b){
    auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    int ov=0;
    for(int k=ov;k<bk.getNz()-ov;++k)
      for(int j=ov;j<bk.getNy()-ov;++j)
        for(int i=ov;i<bk.getNx()-ov;++i){
          std::size_t id=k*pr[2]+j*pr[1]+i; NSCELL c(id,bl);
          for(unsigned kk=0;kk<LatSet::q;++kk){
            T uc=uz*latset::c<LatSet>(kk);
            c[kk]=latset::w<LatSet>(kk)*(pz+LatSet::InvCs2*uc+uc*uc*T{0.5}*LatSet::InvCs4-LatSet::InvCs2*T{0});
          }
        }
  }

  // init MF: psi = -H0*z (uniform field along +z), g_k = w_k*psi
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),mz=bk.getMin()[2];
      for(int k=0;k<bk.getNz();++k){
        T z=mz+T(k)*vs, psi=-H0*z;
        for(int j=0;j<bk.getNy();++j){
          for(int i=0;i<bk.getNx();++i){
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
  MF_Per.Setup(left,NbrDirection::XN,right,NbrDirection::XP);
  MF_Per.Setup(right,NbrDirection::XP,left,NbrDirection::XN);
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
  vtmo::ScalarWriter Hmw("HMAG",MFLattice.getField<HMAG<T>>());
  vtmo::vtmWriter<T,3> MW("bubbleMag3d",Geo);
  MW.addWriterSet(PW,Hmw);

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

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*z  and wall pops = feq(psi_bc)
    // NOTE: lattice row kk holds physical z = minZ + (kk-ov)*vs (ov halo
    // rows below the block). Pinning by absolute z (halo + wall rows) gives
    // the exact linear profile, so H == H0 everywhere in the far field.
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      const T nwall=Cell_Len*T{3.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
        for(int kk=0;kk<nz;++kk){
          T z=minZ+T(kk-ov)*vs;
          if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
            T psi_w=-H0*z;
            for(int jj=0;jj<ny;++jj){
              for(int ii=0;ii<nx;++ii){
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
    auto SyncMFPeriodicGhosts = [&](auto& field, int fidx) {
      const T xLeft = T{0};
      const T xRight = T(Ni) * Cell_Len;
      const T yFront = T{0};
      const T yBack = T(Nj) * Cell_Len;
      const int nEntry = int(GlobalBlockTable.size() / 11);
      const int myRank = mpi().getRank();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& bk = Geo.getBlock(b);
        const bool atL = bk.getMin()[0] < xLeft + Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > xRight - Cell_Len * T{0.5};
        const bool atF = bk.getMin()[1] < yFront + Cell_Len * T{0.5};
        const bool atB = bk.getMax()[1] > yBack - Cell_Len * T{0.5};
        if (!atL && !atR && !atF && !atB) continue;
        const auto& pr = bk.getProjection();
        const int nx = bk.getNx(), ny = bk.getNy(), nz = bk.getNz();
        auto& f1 = field.getBlockField(b);
        // ---- x 方向: partner 同 (minY,minZ), 对侧 x 边 ----
        if (atL || atR) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[2] != bk.getMin()[1] || en[3] != bk.getMin()[2]) continue;
            if (atL && en[8] > T{0.5}) eP = e;  // partner atR
            if (atR && en[7] > T{0.5}) eP = e;  // partner atL
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int pnx = int(enP[4]);
            const int plane = ny * nz;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const int pnx2 = Geo.getBlock(pb).getNx();
                if (atL) for (int kk = 0; kk < nz; ++kk)
                  for (int jj = 0; jj < ny; ++jj)
                    f1.get(kk * pr[2] + jj * pr[1] + 0) = f2.get(kk * pr[2] + jj * pr[1] + pnx2 - 2);
                if (atR) for (int kk = 0; kk < nz; ++kk)
                  for (int jj = 0; jj < ny; ++jj)
                    f1.get(kk * pr[2] + jj * pr[1] + nx - 1) = f2.get(kk * pr[2] + jj * pr[1] + 1);
              }
            } else {
              const int tag = fidx * 1000 + std::min(eL, eP) * 2 + 0;
              const int srcI = atL ? 1 : (nx - 2);
              const int dstI = atL ? 0 : (nx - 1);
              std::vector<T> sbuf(plane), rbuf(plane);
              int kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                  sbuf[kk++] = f1.get(k * pr[2] + j * pr[1] + srcI);
              mpi().sendRecv(sbuf.data(), rbuf.data(), plane, pRank, pRank, tag);
              kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int j = 0; j < ny; ++j)
                  f1.get(k * pr[2] + j * pr[1] + dstI) = rbuf[kk++];
            }
          }
        }
        // ---- y 方向: partner 同 (minX,minZ), 对侧 y 边 ----
        if (atF || atB) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[1] != bk.getMin()[0] || en[3] != bk.getMin()[2]) continue;
            if (atF && en[10] > T{0.5}) eP = e;  // partner atB
            if (atB && en[9]  > T{0.5}) eP = e;  // partner atF
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
            const int pny = int(enP[5]);
            const int plane = nx * nz;
            if (pRank == myRank) {
              int pb = -1;
              for (int bb = 0; bb < Geo.getBlockNum(); ++bb) {
                const auto& bkk = Geo.getBlock(bb);
                if (bkk.getMin()[0] == enP[1] && bkk.getMin()[1] == enP[2] &&
                    bkk.getMin()[2] == enP[3]) { pb = bb; break; }
              }
              if (pb >= 0) {
                auto& f2 = field.getBlockField(pb);
                const int pny2 = Geo.getBlock(pb).getNy();
                if (atF) for (int kk = 0; kk < nz; ++kk)
                  for (int ii = 0; ii < nx; ++ii)
                    f1.get(kk * pr[2] + 0 * pr[1] + ii) = f2.get(kk * pr[2] + (pny2 - 2) * pr[1] + ii);
                if (atB) for (int kk = 0; kk < nz; ++kk)
                  for (int ii = 0; ii < nx; ++ii)
                    f1.get(kk * pr[2] + (ny - 1) * pr[1] + ii) = f2.get(kk * pr[2] + 1 * pr[1] + ii);
              }
            } else {
              const int tag = fidx * 1000 + std::min(eL, eP) * 2 + 1;
              const int srcJ = atF ? 1 : (ny - 2);
              const int dstJ = atF ? 0 : (ny - 1);
              std::vector<T> sbuf(plane), rbuf(plane);
              int kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                  sbuf[kk++] = f1.get(k * pr[2] + srcJ * pr[1] + i);
              mpi().sendRecv(sbuf.data(), rbuf.data(), plane, pRank, pRank, tag);
              kk = 0;
              for (int k = 0; k < nz; ++k)
                for (int i = 0; i < nx; ++i)
                  f1.get(k * pr[2] + dstJ * pr[1] + i) = rbuf[kk++];
            }
          }
        }
      }
    };

    // 0c-0e+: psi-solver sub-iterations (PsiSolver_Iter repeats of
    // collision+stream+macro). Each sub-iteration is one D3Q7 MRT diffusion
    // step with the boosted relaxation omega_psi = 1/(0.5+K*mu); ~100
    // iterations converge the solve to its (smooth) fixed point, removing the
    // grid-scale |H| oscillation that would facet the bubble surface.
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
          T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
          for(int kk=0;kk<nz;++kk){
            T z=minZ+T(kk-ov)*vs;
            if((z<=nwall&&z>=-nwall)||(z<=H_global+nwall&&z>=H_global-nwall)){
              T psi_w=-H0*z;
              for(int jj=0;jj<ny;++jj){
                for(int ii=0;ii<nx;++ii){
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
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>(), 0);

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
    // |z|<=16 and |z-H_global|<=16. (Physically the far-field Kelvin force
    // is ~0 anyway: |H| is uniform, ∇|H|≈0. 3D 三重角伪影比 2D 更宽, 故
    // band 取 16 而非 2D 的 12.)
    if(Bom>T{0}){
      const T band=T{16.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
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

    // wall phi=1 BC
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T mz=bk.getMin()[2],Mz=bk.getMax()[2];
        if(mz<Cell_Len*T{1.5}) for(int jj=0;jj<ny;++jj) for(int i=0;i<nx;++i)bP.get(ov*pr[2]+jj*pr[1]+i)=T{1};
        if(Mz>H_global-Cell_Len*T{1.5}){int kk=nz-1-ov;for(int jj=0;jj<ny;++jj)for(int i=0;i<nx;++i)bP.get(kk*pr[2]+jj*pr[1]+i)=T{1};}
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
        if(Mz>H_global-Cell_Len*T{1.5}) for(int k=nz-ov;k<nz;++k) for(int jj=0;jj<ny;++jj) for(int i=0;i<nx;++i) bP.get(k*pr[2]+jj*pr[1]+i)=T{1};
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
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
