// dropletMag2d.cpp — ferrofluid droplet deformation in a uniform magnetic
// field (paper Sec. III.B, Flament experiment benchmark): validates the
// magnetic Kelvin force. A circular ferrofluid droplet (mu2, chi2) deforms
// into a prolate shape along the field; the aspect ratio b/a vs H0 is
// compared with the literature.
// Phase field + NS + Magnetic, mirroring the rosensweig2d framework.
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"
#include <cstring>

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

// ---- Simulation Parameters ----
int Ni, Nj;
T Cell_Len;
int BlockCellLen, Thread_Num;

// droplet (ferrofluid circle at the center, radius R)
T Droplet_Radius, Droplet_Cx, Droplet_Cy;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase (lattice units, paper Sec. III.D)
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Tau_ns, DeltaRho;

// magnetic field (lattice units; mu, chi are relative / per-cell)
T chi_l, chi_h, mu_l, mu_h, H0, Hc_lat;
T H0_kAm, Hc_kAm;  // physical field strengths (kA/m)

// psi-solver: boosted relaxation tau = 0.5 + PsiSolver_K*mu
int PsiSolver_Iter; T PsiSolver_K;

int MaxStep, OutputStep;
std::string work_dir;

void readParam(int argc, char* argv[]) {
  std::string ininame = "dropletMag2d.ini";
  T H0_kAm_override = T{-1};
  if (argc > 1) ininame = argv[1];
  if (argc > 2) H0_kAm_override = std::atof(argv[2]);
  iniReader r(ininame);
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Droplet_Radius = r.getValue<T>("Droplet","Radius");
  Droplet_Cx = r.getValue<T>("Droplet","CenterX");
  Droplet_Cy = r.getValue<T>("Droplet","CenterY");
  Interface_Width=r.getValue<T>("Phase_Field","Interface_Width");
  Mobility=r.getValue<T>("Phase_Field","Mobility");
  rho_l=r.getValue<T>("Two_Phase","rho_l");
  rho_h=r.getValue<T>("Two_Phase","rho_h");
  eta_l=r.getValue<T>("Two_Phase","eta_l");
  eta_h=r.getValue<T>("Two_Phase","eta_h");
  sigma=r.getValue<T>("Two_Phase","sigma");
  MaxStep=r.getValue<int>("Simulation_Settings","TotalStep");
  if (argc > 3) MaxStep = std::atoi(argv[3]);  // TotalStep override
  OutputStep=r.getValue<int>("Simulation_Settings","OutputStep");
  // Magnetic
  chi_l=r.getValue<T>("Magnetic_Field","chi_l");
  chi_h=r.getValue<T>("Magnetic_Field","chi_h");
  mu_l=r.getValue<T>("Magnetic_Field","mu_l");
  mu_h=r.getValue<T>("Magnetic_Field","mu_h");
  H0_kAm=r.getValue<T>("Magnetic_Field","H0_kAm");
  Hc_kAm=r.getValue<T>("Magnetic_Field","Hc_kAm");
  PsiSolver_Iter=r.getValue<int>("Magnetic_Field","PsiSolver_Iter");
  PsiSolver_K=r.getValue<T>("Magnetic_Field","PsiSolver_K");
  if(H0_kAm_override>=T{0}) H0_kAm=H0_kAm_override;
  // Physical reference (Cowley-Rosensweig unit conversion)
  T Lx_phys_mm = r.getValue<T>("Physical","Lx_mm");
  T sigma_phys_mNm = r.getValue<T>("Physical","sigma_mNm");
  T lambda_c_phys_mm = r.getValue<T>("Physical","lambda_c_mm");

  DeltaRho=rho_h-rho_l;
  Beta=T(12.0)*sigma/Interface_Width;
  Kappa=T(3.0)*Interface_Width*sigma*T(0.5);
  Tau_phi=T(3.0)*Mobility+T(0.5); Omega_phi=T(1.0)/Tau_phi;
  Tau_ns=T(0.5)+eta_h/rho_h/LatSet::cs2;

  // gravity from the Cowley-Rosensweig critical wavelength:
  //   lambda_c = 2*pi*sqrt(sigma/(g*DeltaRho))  ->  g = sigma/DeltaRho*(2*pi/lambda_c)^2
  // with lambda_c_lattice = lambda_c_phys * Ni / Lx_phys (same physical lambda_c).
  T lambda_c_lattice = lambda_c_phys_mm * T(Ni) / Lx_phys_mm;
  gravity = sigma / (DeltaRho * std::pow(lambda_c_lattice / (T{2} * T{M_PI}), T{2}));

  // H0 from the physical field strength via the DIRECT dimensional
  // conversion (mesh-independent): the lattice magnetic field unit follows
  // from the Kelvin-force consistency mu0_conv*H_conv^2 = rho_conv*U^2,
  //   H_conv = sqrt(rho_conv*U^2/mu0_phys)
  // with rho_conv = rho_l_phys/rho_l, U = eta_conv/(rho_conv*dx),
  // eta_conv = eta_l_phys/eta_l, dx = Lx_phys/Ni.
  T rho_l_phys = r.getValue<T>("Physical","rho_l_phys");   // kg/m^3
  T eta_l_phys = r.getValue<T>("Physical","eta_l_phys");   // Pa*s
  T mu0_phys = T{4} * T{M_PI} * T{1e-7};
  T rho_conv = rho_l_phys / rho_l;
  T eta_conv = eta_l_phys / eta_l;
  T dx_phys = Lx_phys_mm * T{1e-3} / T(Ni);
  T U_phys = eta_conv / (rho_conv * dx_phys);
  T H_conv = std::sqrt(rho_conv * U_phys * U_phys / mu0_phys);  // A/m per lattice unit
  Hc_lat = H_conv * Hc_kAm * T{1e3};   // lattice critical field (report only)
  H0 = H0_kAm * T{1e3} / H_conv;

  MPI_RANK(0){
    printf("---- Ferrofluid Droplet Deformation (Paper Sec. III.B) ----\n");
    printf("Mesh: %dx%d  BlockCellLen=%d\n",Ni,Nj,BlockCellLen);
    printf("Droplet: R=%.1f center=(%.0f,%.0f)\n",Droplet_Radius,Droplet_Cx,Droplet_Cy);
    printf("rho: l=%.3f h=%.3f  eta: l=%.5f h=%.5f  sigma=%.5f\n",rho_l,rho_h,eta_l,eta_h,sigma);
    printf("W=%.1f M=%.3f tau_phi=%.3f tau_ns=%.3f\n",Interface_Width,Mobility,Tau_phi,Tau_ns);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) Hc_lat=%.4f\n",chi_l,chi_h,mu_l,mu_h,Hc_lat);
    printf("H0: phys=%.2f kA/m (Hc=%.2f) -> lat=%.4f\n",H0_kAm,Hc_kAm,H0);
    printf("PsiSolver: K=%.3f iter=%d (omega_mu=%.3f omega_1=%.3f)\n",PsiSolver_K,PsiSolver_Iter,
           T{1}/(T{0.5}+PsiSolver_K*mu_h),T{1}/(T{0.5}+PsiSolver_K*mu_l));
    printf("SeamSync: MPI-aware field sync (v2)\n");
    printf("---------------------------------------\n");
  }
}

// ---- PreForceScaled2D ----
// F_p = -(p/3) * DeltaRho * grad_phi * PrC_SCALE
// Example-local copy of ff::FFPreForce2D with the reference-validated
// PrC_SCALE=0.6 (uniform-mesh interface is much thicker than the paper's AMR
// finest level, so the full PrC=1.0 over-damps and suppresses peak growth;
// PrC=0.6 gives the correct nonlinear saturation toward the theoretical peak
// height eta_peak = (lambda_c/4*pi)*sqrt(2*(H0/Hc)^2-2) = 0.454mm).
// The framework functor stays untouched so bubbleMag2d keeps PrC=1.0.
template <typename PFCELL, typename NSCELL>
struct PreForceScaled2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    constexpr T PrC_SCALE = T{0.3};
    T p = ns_cell.template get<PRESSURE<T>>();
    T delta_rho = pf_cell.template get<ff::DELTARHO<T>>();
    const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
    T coeff = -p / T{3} * delta_rho * PrC_SCALE;
    auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
    ns_force[0] += coeff * grad_phi[0];
    ns_force[1] += coeff * grad_phi[1];
  }
};

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Rosensweig Instability..."));
  readParam(argc,argv);

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
  // FIX: SetupBoundary marks ALL domain-boundary cells (including the left/
  // right columns x=0 and x=Ni-1) as BouncebackFlag because their periodic
  // ghost neighbors lie outside the domain AABB. This would treat the left/
  // right boundaries as no-slip walls instead of periodic, causing boundary
  // artifacts and premature peak growth at the seam. Reset those two columns
  // to BulkFlag (the corner rows keep the top/bottom wall BC).
  AABB<T,2> left_col({T{0},Cell_Len},{Cell_Len,T((Nj-1)*Cell_Len)});
  AABB<T,2> right_col({T((Ni-1)*Cell_Len),Cell_Len},{T(Ni*Cell_Len),T((Nj-1)*Cell_Len)});
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

  // -- init phi: circular ferrofluid droplet (phi=1 inside, phi=0 outside) --
  T R_phys=Droplet_Radius*Cell_Len, W_phys=Interface_Width*Cell_Len;
  T xc_d=Droplet_Cx*Cell_Len, yc_d=Droplet_Cy*Cell_Len;
  const T Pi=T{3.14159265358979323846};
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=my+T(j)*vs;
      for(int i=ov;i<bk.getNx()-ov;++i){
        T x=mx+T(i)*vs;
        T dist=std::sqrt((x-xc_d)*(x-xc_d)+(y-yc_d)*(y-yc_d));
        T phi=T{0.5}-T{0.5}*std::tanh(T{2.0}*(dist-R_phys)/W_phys);
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

  // init MF: uniform far field psi = -H0*y (walls are in the outer solvent;
  // the droplet is local and distorts the field only in its vicinity).
  {
    auto& psiF=MFLattice.getField<PSI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bPsi=psiF.getBlockField(b);
      T vs=bk.getVoxelSize(),my=bk.getMin()[1];
      for(int j=0;j<bk.getNy();++j){
        T y=my+T(j)*vs, psi=-H0*y;
        for(int i=0;i<bk.getNx();++i){
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
  using PrT=tmp::Key_TypePair<BulkFlag,PreForceScaled2D<PFCELL,NSCELL>>;
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
  // Keep the default output (1 ghost column per side): ParaView needs the
  // overlapping columns to tile the blocks seamlessly (no gaps between blocks).
  vtmo::vtmWriter<T,2> MW("dropletMag2d",Geo);
  MW.addWriterSet(PW,PS,VW,Dw,Fw,Hxw,Hyw);

  // -- droplet aspect-ratio probe: b/a of the phi=0.5 contour (b: vertical,
  // a: horizontal). Writes droplet_probe.dat (step, a_cells, b_cells, b/a).
  FILE* probeF=nullptr;
  IF_MPI_RANK(0){ probeF=std::fopen("droplet_probe.dat","w");
    std::fprintf(probeF,"# step  a  b  b_over_a\n"); }
  auto WriteProbe = [&](T step){
    T amin=T{1e30},amax=-T{1e30},bmin=T{1e30},bmax=-T{1e30};
    auto& pF=PFLattice.getField<PHI<T>>();
    for(int b=0;b<Geo.getBlockNum();++b){
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      auto& bP=pF.getBlockField(b);
      T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
      int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
      for(int j=ov;j<ny-ov;++j){
        T y=my+T(j)*vs;
        for(int i=ov;i<nx-ov;++i){
          T x=mx+T(i)*vs;
          if(bP.get(j*pr[1]+i)>=T{0.5}){
            if(x<amin)amin=x; if(x>amax)amax=x;
            if(y<bmin)bmin=y; if(y>bmax)bmax=y;
          }
        }
      }
    }
    T a=amax-amin, b=bmax-bmin;
#ifdef MPI_ENABLED
    mpi().reduceAndBcast<T>(amin,MPI_MIN,0);
    mpi().reduceAndBcast<T>(amax,MPI_MAX,0);
    mpi().reduceAndBcast<T>(bmin,MPI_MIN,0);
    mpi().reduceAndBcast<T>(bmax,MPI_MAX,0);
    a=amax-amin; b=bmax-bmin;
#endif
    IF_MPI_RANK(0){
      std::fprintf(probeF,"%g %g %g %g\n",step,a,b,b/a);
      std::fflush(probeF);
    }
  };

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
  const T phiBottom=T{0};  // solvent at the bottom wall
  const T phiTop=T{0};     // solvent at the top wall

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*y (uniform far field, walls in the outer solvent)
    // and wall pops = feq(psi_bc) (Dirichlet pinning)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      const T nwall=Cell_Len*T{3.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minY=bk.getMin()[1],vs=bk.getVoxelSize();
        for(int jj=0;jj<ny;++jj){
          T y=minY+T(jj-ov)*vs;
          if((y<=nwall&&y>=-nwall)||(y<=H_global+nwall&&y>=H_global-nwall)){
            T psi_w=-H0*y;
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
              bPsi.get(id)=psi_w;
            }
          }
        }
      }
    }

    // SyncSeamField: mirror the physical edge columns of a per-cell field
    // across the periodic seam into the opposite x-wrapped ghost columns.
    // The framework's FixedPeriodicBoundaryManager copies only pops + GenericRho,
    // so the wrapped ghost columns of PSI/HX/HY/HMAG/VELOCITY/... stay frozen;
    // and a local-only mirror (bubbleMag2d's SyncMFPeriodicGhosts) fails when
    // the seam blocks live on different ranks. This version handles both:
    // same-rank pairs are copied directly, cross-rank pairs exchange via MPI.
    // (Block-interior ghosts are already filled by the per-field Communicate.)
    auto SyncSeamField = [&](auto& field) {
      using ValT = typename std::decay_t<decltype(field.getBlockField(0))>::value_type;
      const int valSz = static_cast<int>(sizeof(ValT));
      const T Lx = T(Ni) * Cell_Len;
      const int myRank = mpi().getRank();
      const auto& gGeo = GeoHelper.getBlockGeometry();
      const int nb = Geo.getBlockNum();
#ifdef MPI_ENABLED
      std::vector<std::vector<char>> sendBufs, recvBufs;
      std::vector<MPI_Request> sendReqs, recvReqs;
      std::vector<std::pair<int, int>> recvJobs;  // (localBlock, ghostCol0)
#endif
      for (int b = 0; b < nb; ++b) {
        const auto& bk = Geo.getBlock(b);
        const bool atL = bk.getMin()[0] < Cell_Len * T{0.5};
        const bool atR = bk.getMax()[0] > Lx - Cell_Len * T{0.5};
        if (!atL && !atR) continue;
        const int nx = bk.getNx(), ny = bk.getNy(), ov = bk.getOverlap();
        const auto& pr = bk.getProjection();
        // counterpart: global edge block on the opposite side with same y-range
        const T ymid = (bk.getMin()[1] + bk.getMax()[1]) * T{0.5};
        const T xprobe = atL ? (Lx - Cell_Len * T{0.5}) : (Cell_Len * T{0.5});
        int srcGid = -1;
        for (std::size_t sbi = 0; sbi < gGeo.getBlockNum(); ++sbi) {
          const auto& gb = gGeo.getBlock(sbi);
          if (gb.getSelfBlock().isInside(Vector<T, 2>{xprobe, ymid})) {
            srcGid = gb.getBlockId(); break;
          }
        }
        if (srcGid < 0) continue;
        const int ghostCol0 = atL ? 0 : nx - ov;
        const int theirPhysCol = atL ? (nx - 1 - ov) : ov;
        if (GeoHelper.whichRank(srcGid) == myRank) {
          // same-rank: direct copy (counterpart is a local block)
          const int srcBlock = Geo.findBlockIndex(srcGid);
          auto& bF = field.getBlockField(b);
          auto& sF = field.getBlockField(srcBlock);
          for (int j = 0; j < ny; ++j) {
            const ValT v = sF.get(j * pr[1] + theirPhysCol);
            for (int c = 0; c < ov; ++c) bF.get(j * pr[1] + ghostCol0 + c) = v;
          }
        } else {
          // cross-rank: exchange my physical edge column with the counterpart's
#ifdef MPI_ENABLED
          const int myPhysCol = atL ? ov : nx - 1 - ov;
          std::vector<char> snd(ny * valSz);
          auto& bF = field.getBlockField(b);
          for (int j = 0; j < ny; ++j) {
            const ValT v = bF.get(j * pr[1] + myPhysCol);
            std::memcpy(&snd[j * valSz], &v, valSz);
          }
          sendBufs.emplace_back(std::move(snd));
          MPI_Request sreq;
          mpi().iSend(sendBufs.back().data(), static_cast<int>(sendBufs.back().size()),
                      GeoHelper.whichRank(srcGid), &sreq, 9500 + bk.getBlockId());
          sendReqs.push_back(sreq);
          recvBufs.emplace_back(ny * valSz);
          recvJobs.emplace_back(b, ghostCol0);
          MPI_Request rreq;
          mpi().iRecv(recvBufs.back().data(), static_cast<int>(recvBufs.back().size()),
                      GeoHelper.whichRank(srcGid), &rreq, 9500 + srcGid);
          recvReqs.push_back(rreq);
#endif
        }
      }
#ifdef MPI_ENABLED
      MPI_Waitall(static_cast<int>(sendReqs.size()), sendReqs.data(), MPI_STATUSES_IGNORE);
      for (std::size_t i = 0; i < recvReqs.size(); ++i) {
        MPI_Wait(&recvReqs[i], MPI_STATUS_IGNORE);
        auto& bF = field.getBlockField(recvJobs[i].first);
        const auto& bk = Geo.getBlock(recvJobs[i].first);
        const auto& pr = bk.getProjection();
        const int ov = bk.getOverlap();
        const int gCol0 = recvJobs[i].second;
        const auto& buf = recvBufs[i];
        for (int j = 0; j < bk.getNy(); ++j) {
          ValT v;
          std::memcpy(&v, &buf[j * valSz], valSz);
          for (int c = 0; c < ov; ++c) bF.get(j * pr[1] + gCol0 + c) = v;
        }
      }
#endif
    };

    // 0c-0e+: psi-solver sub-iterations (see bubbleMag2d notes)
    for(int sub=0;sub<PsiSolver_Iter;++sub){
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

      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

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

      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        const T nwall=Cell_Len*T{3.0};
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
          T minY=bk.getMin()[1],vs=bk.getVoxelSize();
          for(int jj=0;jj<ny;++jj){
            T y=minY+T(jj-ov)*vs;
            if((y<=nwall&&y>=-nwall)||(y<=H_global+nwall&&y>=H_global-nwall)){
              T psi_w=-H0*y;
              for(int ii=0;ii<nx;++ii){
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }
    }
    SyncSeamField(MFLattice.getField<PSI<T>>());

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
    SyncSeamField(MFLattice.getField<HX<T>>());
    SyncSeamField(MFLattice.getField<HY<T>>());
    SyncSeamField(MFLattice.getField<HMAG<T>>());

    // ===== Phase A: Force setup =====
    RoC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,RoT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    STC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,STT>>(t(),FlagFM);

    // A4: Magnetic force (Kelvin form, mu0=1) — zeroed in wall bands where the
    // pinned psi solver leaves spurious |H| spikes (see bubbleMag2d notes).
    if(H0>T{0}){
      const T band=T{12.0};
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        auto& ns_bl=NSLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        T minY=bk.getMin()[1],vs=bk.getVoxelSize();
        for(int j=ov;j<bk.getNy()-ov;++j){
          T y=minY+T(j-ov)*vs;
          if(y<=band||y>=H_global-band) continue;
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

    // ViC (visco-force correction) disabled: the reference paper model does
    // not use it, and it adds damping that suppresses peak growth (validated
    // in the reference Rosen implementation).
    // ViC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,ViT>>(t(),FlagFM);
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
            // phi floor/cap at 0.01/0.999: prevents trough collapse (numerical
            // divergence when a peak trough becomes extremely thin). Far from
            // the interface transition (phi 0.3-0.7), so no physics is affected.
            if(pn<T{0.01})pn=T{0.01}; if(pn>T{0.999})pn=T{0.999}; bP.get(id)=pn;
          }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi BC: ferrofluid (1) at bottom, solvent (0) at top
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=phiBottom;
        if(My>H_global-Cell_Len*T{1.5}){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=phiTop;}
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // wall phi ghost fix
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=phiBottom;
        if(My>H_global-Cell_Len*T{1.5}) for(int j=ny-ov;j<ny;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=phiTop;
      }
    }
    // seam PHI: the periodic manager's GenericRho copy runs before the macro
    // update, so re-sync the x-wrapped ghosts to the current field
    SyncSeamField(PFLattice.getField<PHI<T>>());

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

    // keep the seam ghost columns of the NS fields consistent for the VTK output
    SyncSeamField(NSLattice.getField<VELOCITY<T,2>>());
    SyncSeamField(NSLattice.getField<PRESSURE<T>>());
    SyncSeamField(NSLattice.getField<FORCE<T,2>>());

    ++t; ++ot;
    if(t()%OutputStep==0){
      ot.Print_InnerLoopPerformance(Geo.getTotalCellNum(),OutputStep);
      Printer::Endl();
      MW.WriteBinary(t());
      WriteProbe(t());
    }
  }

  IF_MPI_RANK(0){ if(probeF) std::fclose(probeF); }
  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  t.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
