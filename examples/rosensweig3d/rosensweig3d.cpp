// rosensweig3d.cpp — 3D Rosensweig instability in a ferrofluid layer
// (Guo et al., Phys. Fluids 37, 022148 (2025), Section III.D, extended to 3D)
//
// Setup: a rectangular cavity, bottom 1/3 ferrofluid (heavy, phi=1) and top 2/3
// organic solvent (light, phi=0), a vertical uniform magnetic field H0, gravity
// downward. Above the critical field Hc the interface breaks into a hexagonal
// array of peaks (normal-field / Rosensweig instability).
//
// Built on the bubbleMag3d framework: phase-field (Allen-Cahn, D3Q19) + NS
// (velocity-based, D3Q19 MRT) + magnetic scalar potential (D3Q7 MRT) solved by
// Lax-Wendroff-style finite-difference LBM. The only physics changes are the
// stratified initial condition (with a small interface perturbation that seeds
// the instability) and the top-wall phase-field boundary (phi=0, organic).
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

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (Rosensweig): organic l / ferrofluid h
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, DeltaRho, Tau_ns;

// Rosensweig specific
T Interface_Z, LambdaC, GravityOverride, Perturb_Amp, Hc;
std::string Perturb_Mode;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0;

int PsiSolver_Iter; T PsiSolver_K;

int MaxStep, OutputStep;
std::string work_dir;

static constexpr T PI = T{3.14159265358979323846};

// ---- initial interface height (physical units) at (x,y) ----
// seeds the instability: 'hex' = hexagonal array at the critical wavelength,
// 'xwave' = parallel ridges (2D-like), 'noise' = deterministic low-mode mix.
static inline T interfaceHeight(T x, T y) {
  T h = Interface_Z;
  const T kc = T{2.0} * PI / LambdaC;
  if (Perturb_Mode == "xwave") {
    h += Perturb_Amp * std::cos(kc * x);
  } else if (Perturb_Mode == "noise") {
    // deterministic sum of near-critical box modes (reproducible across ranks)
    h += Perturb_Amp * 0.5 * (
        std::sin(T{2.0} * PI * (3.0 * x / T(Ni) + 2.0 * y / T(Nj)) + T{0.7})
      + std::sin(T{2.0} * PI * (4.0 * x / T(Ni) - 3.0 * y / T(Nj)) + T{1.3})
      + std::sin(T{2.0} * PI * (5.0 * x / T(Ni) + 1.0 * y / T(Nj)) + T{2.1})
      + std::sin(T{2.0} * PI * (3.0 * x / T(Ni) + 4.0 * y / T(Nj)) + T{2.9})
      + std::sin(T{2.0} * PI * (6.0 * x / T(Ni) - 2.0 * y / T(Nj)) + T{3.7})
      + std::sin(T{2.0} * PI * (2.0 * x / T(Ni) + 5.0 * y / T(Nj)) + T{4.3})
      + std::sin(T{2.0} * PI * (4.0 * x / T(Ni) + 4.0 * y / T(Nj)) + T{5.1})
      + std::sin(T{2.0} * PI * (5.0 * x / T(Ni) - 4.0 * y / T(Nj)) + T{5.9}));
  } else {  // hex: three plane waves at 60° -> hexagonal peak array
    const T c30 = T{0.8660254037844386};
    h += Perturb_Amp * (std::cos(kc * x)
                        + std::cos(kc * (T{0.5} * x + c30 * y))
                        + std::cos(kc * (T{0.5} * x - c30 * y)));
  }
  return h;
}

void readParam() {
  iniReader r("rosensweig3d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Nk = r.getValue<int>("Mesh","Nk");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  eta_l=r.getValue<T>("Two_Phase","eta_l");
  eta_h=r.getValue<T>("Two_Phase","eta_h");
  sigma=r.getValue<T>("Two_Phase","sigma");
  Interface_Z=r.getValue<T>("Rosensweig","Interface_Z");
  LambdaC=r.getValue<T>("Rosensweig","LambdaC");
  GravityOverride=r.getValue<T>("Rosensweig","Gravity");
  Perturb_Amp=r.getValue<T>("Rosensweig","Perturb_Amp");
  Perturb_Mode=r.getValue<std::string>("Rosensweig","Perturb_Mode");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  H0=r.getValue<T>("Magnetic_Field","H0");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");

  DeltaRho=rho_h-rho_l;
  // Gravity from the critical wavelength (Cowley-Rosensweig, paper Eq. 71):
  //   lambda_c = 2*pi*sqrt(sigma/(g*deltaRho))  ->  g = sigma/deltaRho*(2*pi/lambda_c)^2
  if(GravityOverride>T{0}) gravity=GravityOverride;
  else gravity = sigma/DeltaRho * std::pow(T{2.0}*PI/LambdaC, T{2});
  // Critical field (paper Eq. 71, mu0=1):
  //   H_c = sqrt( (2/mu0) * (mu0/mu_h+1)/(mu0/mu_h-1)^2 ) * (sigma*g*deltaRho)^(1/4)
  {
    T r = T{1}/mu_h;  // mu0/mu
    T Hc_pref = std::sqrt(T{2.0}*(r+T{1}) / ((r-T{1})*(r-T{1})));
    Hc = Hc_pref * std::pow(sigma*gravity*DeltaRho, T{0.25});
  }
  Beta=T(12.0)*sigma/Interface_Width;
  Kappa=T(3.0)*Interface_Width*sigma*T(0.5);
  Tau_phi=T(3.0)*Mobility+T(0.5); Omega_phi=T(1.0)/Tau_phi;
  Tau_ns=T(0.5)+eta_h/rho_h/LatSet::cs2;

  MPI_RANK(0){
    printf("---- Rosensweig Instability in Ferrofluid (3D) ----\n");
    printf("Mesh: %dx%dx%d  BlockCellLen=%d\n",Ni,Nj,Nk,BlockCellLen);
    printf("Layer: ferrofluid z<%.1f, organic above  (W=%.1f M=%.3f)\n",Interface_Z,Interface_Width,Mobility);
    printf("rho: l=%.4f h=%.4f  eta: l=%.6f h=%.6f  sigma=%.2e\n",rho_l,rho_h,eta_l,eta_h,sigma);
    printf("gravity=%.3e  lambda_c=%.1f cells  perturb=%s A=%.2f\n",gravity,LambdaC,Perturb_Mode.c_str(),Perturb_Amp);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f)  Hc=%.4f  H0=%.4f  H0/Hc=%.3f\n",
           chi_l,chi_h,mu_l,mu_h,Hc,H0,H0/Hc);
    printf("PsiSolver: K=%.3f iter=%d\n",PsiSolver_K,PsiSolver_Iter);
    printf("TotalStep=%d OutputStep=%d\n",MaxStep,OutputStep);
    printf("-----------------------------------------------\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Rosensweig Instability (3D)..."));
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
  // divide so each block is ~BlockCellLen^3 (scales with the mesh size)
  {
    int bx = std::max(1, Ni / BlockCellLen);
    int by = std::max(1, Nj / BlockCellLen);
    int bz = std::max(1, Nk / BlockCellLen);
    GeoHelper.CreateBlocks(bx, by, bz);
  }
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

  // -- init phi (3D stratified layer: ferrofluid phi=1 below interface, organic phi=0 above) --
  T W_phys=Interface_Width*Cell_Len;
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1],mz=bk.getMin()[2];
    int ov=bk.getOverlap();
    for(int k=ov;k<bk.getNz()-ov;++k){
      T z=mz+T(k)*vs;
      for(int j=ov;j<bk.getNy()-ov;++j){
        T y=my+T(j)*vs;
        for(int i=ov;i<bk.getNx()-ov;++i){
          T x=mx+T(i)*vs;
          T zif=interfaceHeight(x,y);
          T phi=T{0.5}-T{0.5}*std::tanh(T{2.0}*(z-zif)/W_phys);
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

  // init NS pops with the hydrostatic pressure profile (p = g·∫_z^H ρ dz'),
  // u=0. This removes the settling transient a p=0 start would create, so the
  // interface begins near equilibrium and the magnetic instability can develop
  // from a quiet base.
  Vector<T,3> uz{0,0,0};
  {
    auto& dF=NSLattice.getField<DENSITY<T>>();
    auto& pF=NSLattice.getField<PRESSURE<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPhi=phiField.getBlockField(b);
      auto& bD=dF.getBlockField(b); auto& bP=pF.getBlockField(b);
      int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
      T vs=bk.getVoxelSize(), mz=bk.getMin()[2];
      for(int j=ov;j<ny-ov;++j) for(int i=ov;i<nx-ov;++i){
        // walk k from top to bottom, accumulating the hydrostatic pressure
        T pacc=0;
        for(int k=nz-1-ov;k>=ov;--k){
          std::size_t id=k*pr[2]+j*pr[1]+i;
          T phi=bPhi.get(id);
          T rho=rho_l+phi*(rho_h-rho_l);
          // half-cell integration: p = g * rho * (H - z) accumulated top-down
          T p = pacc + T{0.5}*gravity*rho*vs;
          pacc += gravity*rho*vs;
          NSCELL c(id,bl);
          for(unsigned kk=0;kk<LatSet::q;++kk){
            T uc=uz*latset::c<LatSet>(kk);
            c[kk]=latset::w<LatSet>(kk)*(p+LatSet::InvCs2*uc+uc*uc*T{0.5}*LatSet::InvCs4-LatSet::InvCs2*T{0});
          }
          bD.get(id)=rho; bP.get(id)=p;
        }
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
  vtmo::ScalarWriter PSw("PSI",MFLattice.getField<PSI<T>>());
  vtmo::VectorWriter VW("Velocity",NSLattice.getField<VELOCITY<T,3>>());
  vtmo::vtmWriter<T,3> MW("rosensweig3d",Geo);
  MW.addWriterSet(PW,Hmw,PSw,VW);

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
    // field into the opposite x/y-wrapped ghost planes (needed because
    // MF_Per.Apply() only copies pops + GenericRho=PHI, leaving the wrapped
    // PSI/HX/HY/HZ/HMAG ghosts stale). 3D/MPI: partners may be on other ranks,
    // located via the global block table (see bubbleMag3d for details).
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
        if (atL || atR) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[2] != bk.getMin()[1] || en[3] != bk.getMin()[2]) continue;
            if (atL && en[8] > T{0.5}) eP = e;
            if (atR && en[7] > T{0.5}) eP = e;
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
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
        if (atF || atB) {
          int eL = -1, eP = -1;
          for (int e = 0; e < nEntry; ++e) {
            const double* en = &GlobalBlockTable[e * 11];
            if (en[1] != bk.getMin()[0] || en[3] != bk.getMin()[2]) continue;
            if (atF && en[10] > T{0.5}) eP = e;
            if (atB && en[9]  > T{0.5}) eP = e;
            if (int(en[0]) == myRank && en[1] == bk.getMin()[0] &&
                en[2] == bk.getMin()[1] && en[3] == bk.getMin()[2])
              eL = e;
          }
          if (eP >= 0 && eL >= 0) {
            const double* enP = &GlobalBlockTable[eP * 11];
            const int pRank = int(enP[0]);
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

    // 0c-0e+: psi-solver sub-iterations (collision+stream+macro per iteration)
    for(int sub=0;sub<PsiSolver_Iter;++sub){
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

      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

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

      // re-pin walls so the stored field is exactly -H0*z there before H is computed
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

    // A4: Magnetic force (paper Eq. (8): F_m = (μ0 χ/2)∇(|H|²) = χ|H|∇|H|, μ0=1).
    // Force-free wall bands: the z-wall ψ-pins leave a spurious |H| spike within
    // ~6 rows of each wall (solver artifact). The interface sits at z=Interface_Z,
    // far from both walls, so zeroing the wall bands is safe.
    if(H0>T{0}){
      const T band=T{8.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        T minZ=bk.getMin()[2],vs=bk.getVoxelSize();
        for(int k=ov;k<bk.getNz()-ov;++k){
          T z=minZ+T(k)*vs;
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

    // wall phi BC: bottom wall ferrofluid (phi=1), top wall organic (phi=0)
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

    // ===== Diagnostics: interface peak growth =====
    ++t; ++ot;
    if(t()%OutputStep==0){
      // Global max/min interface height (phi=0.5 crossing along z) via MPI allreduce.
      // Local scan: each interior (i,j) column crosses phi 1->0 once (ferrofluid
      // below, organic above); interpolate the crossing z and track the global
      // max (peak crest) and min (valley).
      T zmax=-1e30, zmin=1e30;
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        T vs=bk.getVoxelSize(), mz=bk.getMin()[2];
        for(int j=ov;j<ny-ov;++j) for(int i=ov;i<nx-ov;++i){
          T prev=0, prevz=0;
          // scan k from the ghost below to the ghost above so crossings that
          // fall on a block z-boundary (interface near a block edge) are caught
          for(int k=ov-1;k<=nz-ov;++k){
            T ph=bP.get(k*pr[2]+j*pr[1]+i);
            T zz=mz+T(k)*vs;
            if(k>ov-1 && prev>T{0.5} && ph<T{0.5}){
              T w=(T{0.5}-prev)/(ph-prev);
              T cr=prevz+w*(zz-prevz);
              if(cr>zmax) zmax=cr;
              if(cr<zmin) zmin=cr;
            }
            prev=ph;prevz=zz;
          }
        }
      }
      T umax=0, hmax=0; double phisum=0, phicnt=0;
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int nx=bk.getNx(),ny=bk.getNy(),nz=bk.getNz(),ov=bk.getOverlap();
        auto& bP=pF.getBlockField(b);
        auto& bV=NSLattice.getField<VELOCITY<T,3>>().getBlockField(b);
        auto& bH=MFLattice.getField<HMAG<T>>().getBlockField(b);
        for(int k=ov;k<nz-ov;++k) for(int j=ov;j<ny-ov;++j) for(int i=ov;i<nx-ov;++i){
          std::size_t id=k*pr[2]+j*pr[1]+i;
          auto v=bV.get(id); T um=std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); if(um>umax)umax=um;
          T hm=bH.get(id); if(hm>hmax)hmax=hm;
          phisum+=bP.get(id); phicnt+=1.0;
        }
      }
#ifdef MPI_ENABLED
      T gzmax=zmax,gzmin=zmin,gumax=umax,ghmax=hmax;
      MPI_Allreduce(&zmax,&gzmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      MPI_Allreduce(&zmin,&gzmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
      MPI_Allreduce(&umax,&gumax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      MPI_Allreduce(&hmax,&ghmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
      double gphisum=phisum,gphicnt=phicnt;
      MPI_Allreduce(&phisum,&gphisum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&phicnt,&gphicnt,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#else
      T gzmax=zmax,gzmin=zmin,gumax=umax,ghmax=hmax;
      double gphisum=phisum,gphicnt=phicnt;
#endif
      IF_MPI_RANK(0){
        ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
        T hnorm = H0>T{0} ? ghmax/H0 : ghmax;
        printf("  step=%6d  z_if: max=%7.2f min=%7.2f amp=%7.2f  |u|max=%.4f  |H|max/H0=%.3f  <phi>=%.3f\n",
               (int)t(),gzmax,gzmin,gzmax-gzmin,gumax,hnorm, gphicnt>0?gphisum/gphicnt:T{0});
        Printer::Endl();
      }
      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
