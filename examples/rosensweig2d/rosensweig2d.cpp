// rosensweig2d.cpp — 2D Rosensweig instability of ferrofluid (Phase field + NS + Magnetic)
// Reference: Guo et al. (2025) Phys. Fluids 37, 022148, Sec. III.D (Case D)
//
// Reuses the bubbleMag2d framework:
//   - Phase field: second-order conservative Allen-Cahn (D2Q9 MRT)
//   - Flow field:   velocity-based incompressible NS      (D2Q9 MRT, He-Luo)
//   - Magnetic field: scalar potential ψ with pseudo-time  (D2Q5 MRT, μ₀≡1)
//
// Differences from bubbleMag2d (Case C → Case D):
//   * Geometry: 3:1 rectangular cavity (Lx=21mm, Ly=7mm). No bubble.
//   * Initial φ: heavy ferrofluid (φ=1) fills bottom 1/3, light organic (φ=0)
//                fills top 2/3. Interface y₀ = 2Ly/3 is perturbed by a small
//                cosine of N modes across Lx to seed the instability.
//   * Wall φ BC: bottom = 1 (heavy), top = 0 (light).  (bubbleMag2d: both = 1)
//   * Gravity:   derived from Cowley-Rosensweig critical wavelength
//                λ_c = 2π√(σ/(g·Δρ)), NOT from Re/Eo.
//   * H0:        derived from physical H₀ (kA/m) via critical-field consistency
//                Bo_m_c = μ₀ H_c² λ_c / (2σ) preserved between physical & lattice.
//
// Boundary conditions (paper Sec. III.D):
//   - Top/bottom: no-slip wall (NS), φ fixed (PF), ψ = -H₀·y Dirichlet (MF)
//   - Left/right: periodic for all three fields
//
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

// Rosensweig geometry & initial perturbation
T Interface_Y_Ratio;   // y_interface = ratio * Ly   (default 2/3: heavy at bottom 1/3)
T Perturb_Amplitude;   // amplitude in lattice cells (small, e.g. 0.5)
int Perturb_Modes;     // number of cosine peaks across Lx

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (lattice units, Guo2025 Case D Table)
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Tau_ns, DeltaRho;

// magnetic field (Guo2025 Case D)
T chi_l, chi_h, mu_l, mu_h, H0, H0_phys_kAm, Bom;

// physical reference for unit conversion (Cowley-Rosensweig)
T Lx_phys_mm, Ly_phys_mm;
T sigma_phys_mNm, Hc_phys_kAm, lambda_c_phys_mm;

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("rosensweig2d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");

  Interface_Y_Ratio = r.getValue<T>("Rosensweig","Interface_Y_Ratio");
  Perturb_Amplitude = r.getValue<T>("Rosensweig","Perturb_Amplitude");
  Perturb_Modes = r.getValue<int>("Rosensweig","Perturb_Modes");

  Interface_Width = r.getValue<T>("Phase_Field","Interface_Width");
  Mobility = r.getValue<T>("Phase_Field","Mobility");
  rho_l = r.getValue<T>("Two_Phase","rho_l");
  rho_h = r.getValue<T>("Two_Phase","rho_h");
  eta_l = r.getValue<T>("Two_Phase","eta_l");
  eta_h = r.getValue<T>("Two_Phase","eta_h");
  sigma = r.getValue<T>("Two_Phase","sigma");

  chi_l = r.getValue<T>("Magnetic_Field","chi_l");
  chi_h = r.getValue<T>("Magnetic_Field","chi_h");
  mu_l = r.getValue<T>("Magnetic_Field","mu_l");
  mu_h = r.getValue<T>("Magnetic_Field","mu_h");
  H0_phys_kAm = r.getValue<T>("Magnetic_Field","H0_phys_kAm");

  Lx_phys_mm = r.getValue<T>("Physical","Lx_mm");
  Ly_phys_mm = r.getValue<T>("Physical","Ly_mm");
  sigma_phys_mNm = r.getValue<T>("Physical","sigma_mNm");
  Hc_phys_kAm = r.getValue<T>("Physical","Hc_kAm");
  lambda_c_phys_mm = r.getValue<T>("Physical","lambda_c_mm");

  MaxStep = r.getValue<int>("Simulation_Settings","TotalStep");
  OutputStep = r.getValue<int>("Simulation_Settings","OutputStep");

  // ---- Derived: lattice gravity from Cowley-Rosensweig λ_c ----
  // λ_c = 2π √(σ / (g · Δρ))
  // Physical λ_c (mm) → lattice cells: λ_c_lattice = λ_c_phys × Ni / Lx_phys
  DeltaRho = rho_h - rho_l;
  T lambda_c_lattice = (lambda_c_phys_mm * T(Ni)) / Lx_phys_mm;
  T half_lc_over_2pi = lambda_c_lattice / (T{2} * T{M_PI});
  gravity = sigma / (DeltaRho * half_lc_over_2pi * half_lc_over_2pi);

  // ---- Derived: lattice H0 from physical H0 (kA/m) ----
  // Use critical-field consistency: Bo_m_c is dimensionless, preserved across units.
  //   Bo_m_c_phys = μ₀ H_c² λ_c / (2 σ)       (SI)
  //   Bo_m_c_lat  = H_c_lat² λ_c_lat / (2 σ_lat)   (μ₀ ≡ 1 in lattice)
  // Equating → H_c_lat = H_c_phys × √( μ₀ σ_lat λ_c_phys / (σ_phys λ_c_lat) )
  T mu0_phys = T{4} * T{M_PI} * T{1e-7};
  T lambda_c_phys_m = lambda_c_phys_mm * T{1e-3};
  T sigma_phys = sigma_phys_mNm * T{1e-3};        // N/m
  T Hc_phys_Am = Hc_phys_kAm * T{1e3};            // A/m
  T Hc_lattice = std::sqrt(mu0_phys * Hc_phys_Am * Hc_phys_Am *
                           lambda_c_phys_m * sigma /
                           (sigma_phys * lambda_c_lattice));
  T H_conv = Hc_lattice / Hc_phys_kAm;            // lattice units per (kA/m)
  H0 = H0_phys_kAm * H_conv;
  // Reporting Bo_m using λ_c as the characteristic length
  Bom = H0 * H0 * lambda_c_lattice / (T{2} * sigma);

  Beta = T{12.0} * sigma / Interface_Width;
  Kappa = T{3.0} * Interface_Width * sigma * T{0.5};
  Tau_phi = T{3.0} * Mobility + T{0.5}; Omega_phi = T{1.0} / Tau_phi;
  Tau_ns = T{0.5} + eta_h / rho_h / LatSet::cs2;

  MPI_RANK(0){
    Printer::Print_BigBanner("Rosensweig Instability (Guo2025 Case D)");
    Printer::PrintTitle("Mesh");
    Printer::Print("Ni", Ni);
    Printer::Print("Nj", Nj);
    Printer::Print("BlockCellLen", BlockCellLen);
    Printer::Print("Aspect", T(Ni)/T(Nj));
    Printer::Endl();
    Printer::PrintTitle("Interface");
    Printer::Print("y0_ratio", Interface_Y_Ratio);
    Printer::Print("Perturb_A", Perturb_Amplitude);
    Printer::Print("Modes", Perturb_Modes);
    Printer::Endl();
    Printer::PrintTitle("Phase Field");
    Printer::Print("W", Interface_Width);
    Printer::Print("M", Mobility);
    Printer::Print("tau_phi", Tau_phi);
    Printer::Print("Beta", Beta);
    Printer::Print("Kappa", Kappa);
    Printer::Endl();
    Printer::PrintTitle("Two Phase");
    Printer::Print("rho_l", rho_l);
    Printer::Print("rho_h", rho_h);
    Printer::Print("eta_l", eta_l);
    Printer::Print("eta_h", eta_h);
    Printer::Print("sigma", sigma);
    Printer::Print("d_rho", DeltaRho);
    Printer::Endl();
    Printer::PrintTitle("Magnetic");
    Printer::Print("chi_l", chi_l);
    Printer::Print("chi_h", chi_h);
    Printer::Print("mu_l", mu_l);
    Printer::Print("mu_h", mu_h);
    Printer::Endl();
    Printer::PrintTitle("Physical Reference");
    Printer::Print("Lx_mm", Lx_phys_mm);
    Printer::Print("Ly_mm", Ly_phys_mm);
    Printer::Print("sigma_mNm", sigma_phys_mNm);
    Printer::Print("Hc_kAm", Hc_phys_kAm);
    Printer::Print("lambda_c_mm", lambda_c_phys_mm);
    Printer::Endl();
    Printer::PrintTitle("Derived");
    Printer::Print("gravity", gravity);
    Printer::Print("H0_lattice", H0);
    Printer::Print("H0_phys_kAm", H0_phys_kAm);
    Printer::Print("Bo_m", Bom);
    Printer::Print("lambda_c_lat", lambda_c_lattice);
    Printer::Endl();
    Printer::Print("tau_ns", Tau_ns);
    Printer::Print("Hc_lattice", Hc_lattice);
    Printer::Print("H_conv", H_conv);
    Printer::Endl();
    Printer::PrintTitle("Steps");
    Printer::Print("Total", MaxStep);
    Printer::Print("Output", OutputStep);
    Printer::Endl();
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(("Initializing Rosensweig Instability..."));
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

  // -- flag: bulk interior + periodic left/right + bounceback top/bottom --
  BlockFieldManager<FLAG,T,2> FlagFM(Geo,VoidFlag);
  FlagFM.forEach(domain,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(left,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.forEach(right,[&](FLAG&f,std::size_t id){f.SetField(id,PeriodicFlag);});
  FlagFM.template SetupBoundary<LatSet>(domain,BouncebackFlag);

  // FIX: SetupBoundary marks ALL domain-boundary cells (including left/right
  // columns x=0 and x=Ni-1) as BouncebackFlag because their periodic ghost
  // neighbors lie outside the domain AABB. This would incorrectly treat
  // left/right boundaries as no-slip walls instead of periodic.
  // Reset left/right boundary columns back to BulkFlag (keeping top/bottom
  // rows as BouncebackFlag for the no-slip walls).
  // AABBs exclude y=0 and y=Nj-1 to preserve top/bottom wall BC at corners.
  AABB<T,2> left_col({T(0),            Cell_Len},
                     {Cell_Len,         T((Nj-1)*Cell_Len)});
  AABB<T,2> right_col({T((Ni-1)*Cell_Len), Cell_Len},
                      {T(Ni*Cell_Len),     T((Nj-1)*Cell_Len)});
  FlagFM.forEach(left_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});
  FlagFM.forEach(right_col,[&](FLAG&f,std::size_t id){f.SetField(id,BulkFlag);});

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
    MU_L<T>,MU_H<T>,CHI_L<T>,CHI_H<T>,H_0<T>>;
  using MFREF=TypePack<PHI<T>>;
  using MFPACK=TypePack<MFFIELDS,MFREF>;
  ValuePack MFI(T{},T{1.0},T{mu_l},T{chi_l},T{},T{},T{},T{},
    mu_l,mu_h,chi_l,chi_h,H0);
  using MFCELL=Cell<T,MFLatSet,ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T,MFLatSet,MFPACK> MFLattice(Geo,MFI,MFBaseConv,
    &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams<T>(MFLattice,mu_l,mu_h,chi_l,chi_h,H0);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- init phi: layered heavy/light with cosine perturbation at y_interface --
  //   y_interface(x) = Interface_Y_Ratio * Ly + A * cos(2π * n * x / Lx)
  //   φ(x,y) = 0.5 + 0.5 * tanh(2 * (y_interface - y) / W)
  //   → φ=1 below interface (heavy ferrofluid), φ=0 above (light organic)
  T Ly_phys = T(Nj) * Cell_Len;
  T Lx_phys = T(Ni) * Cell_Len;
  T y0 = Interface_Y_Ratio * Ly_phys;
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
        T y_int = y0 + Perturb_Amplitude * Cell_Len *
                  std::cos(T{2}*T{M_PI}*T(Perturb_Modes)*x/Lx_phys);
        // phi=1 below interface (heavy), phi=0 above (light)
        T phi=T{0.5}+T{0.5}*std::tanh(T{2.0}*(y_int-y)/W_phys);
        bPhi.get(j*pr[1]+i)=phi;
      }
    }
  }
  // init PF pops (feq with u=0)
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

  // init MF: flat-interface solution + analytical evanescent perturbation correction
  //
  // The perturbation-induced field correction is an evanescent mode
  // (exp(-k|y-y0|)) that is an EXACT FIXED POINT (rho=1.0) of the D2Q5 LBM
  // diffusion solver. The solver CANNOT develop it from a zero initial
  // condition, but it WILL preserve it if seeded. By initialising psi with the
  // analytical correction, the correct magnetic force variation is present from
  // step 0 and maintained throughout, driving the Rosensweig instability.
  //
  // Flat-interface solution (consistent with wall BCs):
  //   Below (y<y0, ferrofluid): psi = -Hbelow*y,  Hbelow = H0*mu_l/mu_h
  //   Above (y>y0, light fluid): psi = psi(y0) - H0*(y-y0)
  //
  // Evanescent correction (eta = A*cos(kx), k = 2*pi*n/Lx):
  //   Above: dpsi = Ca*cos(kx)*exp(-k(y-y0))
  //   Below: dpsi = Cb*cos(kx)*exp( k(y-y0))
  //   Ca = A*H0*(1-1/mu_r)*mu_h/(mu_h+mu_l)
  //   Cb = -(mu_l/mu_h)*Ca
  {
    T Hbelow = H0 * mu_l / mu_h;
    T psi_y0 = -Hbelow * y0;
    T wk = T{2} * T{M_PI} * T(Perturb_Modes) / Lx_phys;
    T mu_ratio = mu_l / mu_h;
    T Ca = Perturb_Amplitude * Cell_Len * H0 * (T{1} - mu_ratio)
           * mu_h / (mu_h + mu_l);
    T Cb = -mu_ratio * Ca;
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),my=bk.getMin()[1],mx=bk.getMin()[0];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs;
        for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs;
          T psi_flat, dpsi;
          if(y >= y0){
            psi_flat = psi_y0 - H0 * (y - y0);
            dpsi = Ca * std::cos(wk * x) * std::exp(-wk * (y - y0));
          } else {
            psi_flat = -Hbelow * y;
            dpsi = Cb * std::cos(wk * x) * std::exp(wk * (y - y0));
          }
          T psi = psi_flat + dpsi;
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

  // -- Coupling tasks (PF↔NS) --
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

  // MF coupling: PF→MF (per-cell μ, χ, ω_ψ update)
  using MCT=tmp::Key_TypePair<BulkFlag,MFUpdateCoeffs2D<PFCELL,MFCELL>>;
  using MCSel=CoupledTaskSelector<std::uint8_t,PFCELL,MFCELL,MCT>;
  BlockLatManagerCoupling MCC(PFLattice,MFLattice);

  // Writers
  vtmo::ScalarWriter PW("PHI",PFLattice.getField<PHI<T>>());
  vtmo::vtmWriter<T,2> MW("rosensweig2d",Geo);
  MW.addWriterSet(PW);

  // ===== initial setup =====
  PFLattice.NormalFullCommunicate(); NSLattice.NormalFullCommunicate(); MFLattice.NormalFullCommunicate();
  NS_Per.Apply(); PF_Per.Apply(); MF_Per.Apply();

  PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
  PFLattice.getField<NORMAL<T,2>>().Communicate();
  PFLattice.getField<GRAD<T,2>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);
  RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(0,FlagFM);
  MW.WriteBinary(0);

  Printer::Print_BigBanner(("Start Calculation..."));
  Timer t; Timer ot;
  T H_global=T(Nj)*Cell_Len;
  // Flat-interface MF constants for wall BCs
  T y0_lat=Interface_Y_Ratio*H_global;
  T Hbelow=H0*mu_l/mu_h;
  T psi_at_y0=-Hbelow*y0_lat;
  // Evanescent correction constants for wall BCs (must match init block)
  // Prevents the wall BC from draining the seeded evanescent perturbation
  T wk_bc=T{2}*T{M_PI}*T(Perturb_Modes)/(T(Ni)*Cell_Len);
  T mu_ratio_bc=mu_l/mu_h;
  T Ca_bc=Perturb_Amplitude*Cell_Len*H0*(T{1}-mu_ratio_bc)*mu_h/(mu_h+mu_l);
  T Cb_bc=-mu_ratio_bc*Ca_bc;

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b-0e: MF pseudo-time iteration (multiple sub-iterations per step)
    // Each sub-iteration: set wall BC → collision → stream → update PSI
    // Multiple sub-iterations ensure the magnetic field converges to the
    // steady-state corresponding to the current phi distribution, which is
    // critical for accurate magnetic force computation.
    constexpr int MF_SUBSTEPS = 30;
    for(int sub=0; sub<MF_SUBSTEPS; ++sub){
      // 0b: Set wall psi using flat-interface solution + evanescent correction
      // Bottom wall (y≈0, ferrofluid): ψ = -Hbelow·y + Cb·cos(kx)·exp(k(y-y0))
      // Top wall (y≈Ly, light fluid): ψ = ψ(y0)-H0·(y-y0) + Ca·cos(kx)·exp(-k(y-y0))
      // The evanescent correction prevents the wall BC from draining the seeded
      // perturbation, improving correction retention from ~51% toward ~100%.
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
          T minY=bk.getMin()[1],maxY=bk.getMax()[1],vs=bk.getVoxelSize(),mX=bk.getMin()[0];
          if(minY<Cell_Len*T{1.5}){
            for(int jj=0;jj<=ov;++jj){
              T y=minY+T(jj)*vs;
              for(int ii=0;ii<nx;++ii){
                T x=mX+T(ii)*vs;
                T psi_w=-Hbelow*y + Cb_bc*std::cos(wk_bc*x)*std::exp(wk_bc*(y-y0_lat));
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
          if(maxY>H_global-Cell_Len*T{1.5}){
            for(int jj=ny-ov-1;jj<ny;++jj){
              T y=minY+T(jj)*vs;
              for(int ii=0;ii<nx;++ii){
                T x=mX+T(ii)*vs;
                T psi_w=psi_at_y0-H0*(y-y0_lat) + Ca_bc*std::cos(wk_bc*x)*std::exp(-wk_bc*(y-y0_lat));
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }

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
      MF_Per.Apply();

      // 0d: Stream (pre-stream comm syncs ghost pops; post-stream comm omitted —
      // ghost pops are not read until next substep's pre-stream comm overwrites them)
      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();

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
      // CommunicatePSI omitted within substep loop: collision & wall BC don't
      // read ghost PSI; ghost PSI is properly synced after the loop (lines below).
    }

    // Periodic ghost PSI sync: Apply() only copies POP+PHI; MFComputeH2D reads
    // neighbor PSI, so left/right periodic ghost cells must mirror the opposite
    // interior PSI instead of holding the stale initial value (-H0*y).
    MF_Per.ApplyField(MFLattice.getField<PSI<T>>());

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

    // Periodic ghost H-field sync: MFMagneticForce2D reads neighbor HMAG, so the
    // left/right periodic ghost HX/HY/HMAG must mirror the opposite interior
    // values. Without this, x=0 cells read a stale HMAG=0 at the x=-1 ghost.
    MF_Per.ApplyField(MFLattice.getField<HX<T>>());
    MF_Per.ApplyField(MFLattice.getField<HY<T>>());
    MF_Per.ApplyField(MFLattice.getField<HMAG<T>>());

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic Kelvin force F_m = χ|H|∇|H|  (only when H0 > 0)
    if(H0>T{0}){
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
    // PrC re-enabled: pressure correction for incompressible LBM with density contrast.
    // Without PrC, peaks grow unbounded (4.45mm, domain-limited) at both W=4 and W=8.
    // PrC is essential for correct nonlinear saturation. With W=12, PrC damping (~1/W)
    // is reduced enough to allow slow growth while still providing saturation damping.
    PrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,PrT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().Communicate();

    // ===== Phase B-C: PF + NS collision =====
    PFLattice.template ApplyInnerCellDynamics<PFSel>(FlagFM);
    PF_Per.Apply(); PFLattice.NormalFullCommunicate();

    // ViC disabled: paper (Guo2025) does not use visco-force correction.
    // Testing if this reduces the additional damping that prevents growth at H0=8.2.
    // ViC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,ViT>>(t(),FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSSel>(FlagFM);
    NS_Per.Apply(); NSLattice.NormalFullCommunicate();

    // ===== Phase D: Streaming =====
    PF_BB.Apply(t());
    // PF Y ghost pop fix: copy wall pops to ghost (no-flux BC, same as bubbleMag2d)
    // The phi FIELD is already set to Dirichlet values (BOTTOM=1, TOP=0) in the
    // wall phi BC section below. For the POPULATIONS, we use the same ghost-pop
    // copy as bubbleMag2d, which creates a no-flux (∂φ/∂n=0) BC for the pops.
    // This is sufficient because the ff2d.hh mask fix prevents bulk phi leakage
    // (the original motivation for the Dirichlet pop override).
    {
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

    // wall phi BC: BOTTOM=1 (heavy ferrofluid), TOP=0 (light organic) — differs from bubbleMag2d
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=T{1};
        if(My>H_global-Cell_Len*T{1.5}){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=T{0};}
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi ghost fix (BOTTOM=1, TOP=0)
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{1};
        if(My>H_global-Cell_Len*T{1.5}) for(int j=ny-ov;j<ny;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{0};
      }
    }

    // gradients (normal, laplacian, chem potential)
    PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
    PFLattice.getField<NORMAL<T,2>>().Communicate();
    PFLattice.getField<GRAD<T,2>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // wall grad/chempot fix (one-sided 2nd-order extrapolation; same form as bubbleMag2d)
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

      // Diagnostic: check HMAG range to verify MF solver convergence
      // H0_lattice is the expected |H| for a uniform field. If min/max HMAG
      // are both ≈ H0, the field is uniform (no interface distortion captured).
      // Deviation from H0 indicates the MF solver is capturing interface effects.
      {
        T hmag_max=0, hmag_min=1e10;
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& mf_bl=MFLattice.getBlockLat(b);
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          int ov=bk.getOverlap();
          for(int j=ov;j<bk.getNy()-ov;++j)
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=j*pr[1]+i;
              MFCELL mf(id,mf_bl);
              T hm=mf.template get<HMAG<T>>();
              if(hm>hmag_max)hmag_max=hm;
              if(hm<hmag_min)hmag_min=hm;
            }
        }
        Printer::Print("HMAG_min", hmag_min);
        Printer::Print("HMAG_max", hmag_max);
        Printer::Print("H0_ref", H0);

        // Diagnostic: HMAG variation at interface height (measures field
        // concentration that drives the instability) and interface amplitude
        {
          T hmag_if_max=-1e10, hmag_if_min=1e10;
          T phi_peak=0, phi_trough=1;
          int j_if=int(y0_lat+T{0.5});
          for(int b=0;b<Geo.getBlockNum();++b){
            auto& mf_bl=MFLattice.getBlockLat(b);
            auto& pf_bl=PFLattice.getBlockLat(b);
            const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
            int ov=bk.getOverlap();
            T my=bk.getMin()[1];
            int j_local=int(y0_lat-my+T{0.5});
            if(j_local>=ov && j_local<bk.getNy()-ov){
              for(int i=ov;i<bk.getNx()-ov;++i){
                std::size_t id=j_local*pr[1]+i;
                MFCELL mf(id,mf_bl);
                PFCELL pf(id,pf_bl);
                T hm=mf.template get<HMAG<T>>();
                T ph=pf.template get<PHI<T>>();
                if(hm>hmag_if_max)hmag_if_max=hm;
                if(hm<hmag_if_min)hmag_if_min=hm;
                if(ph>phi_peak)phi_peak=ph;
                if(ph<phi_trough)phi_trough=ph;
              }
            }
          }
          T hmag_var=hmag_if_max-hmag_if_min;
          Printer::Print("HMAG_if_var", hmag_var);
          Printer::Print("HMAG_if_max", hmag_if_max);
          Printer::Print("HMAG_if_min", hmag_if_min);
          Printer::Print("Phi_peak", phi_peak);
          Printer::Print("Phi_trough", phi_trough);
        }
      }

      // Sync periodic ghost cells right before VTK output.
      // Phase A–E above only refresh INTERIOR cells; the left/right periodic
      // ghost cells of every lattice (PF, NS, MF) still hold one-step-lagged
      // values from the previous MF_Per.ApplyField / PF_Per.Apply calls in the
      // loop. NormalFullCommunicate()/Communicate() in Phase B–E sync overlap
      // between blocks but can overwrite periodic ghosts on cross-process
      // boundaries, so we re-apply periodic ghost sync for all visualized
      // fields here to ensure VTK output is one-step-consistent.
      PF_Per.ApplyField(PFLattice.getField<PHI<T>>());

      MF_Per.ApplyField(MFLattice.getField<PSI<T>>());
      MF_Per.ApplyField(MFLattice.getField<HX<T>>());
      MF_Per.ApplyField(MFLattice.getField<HY<T>>());
      MF_Per.ApplyField(MFLattice.getField<HMAG<T>>());

      MW.WriteBinary(t());
    }
  }

  Printer::Print_BigBanner(("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
