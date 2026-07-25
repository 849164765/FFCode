// bubbleMag2d.cpp — 2D bubble rising in ferrofluid (Phase field + NS + Magnetic)
// Force formula selector:
//   0 = Full Maxwell: (H·∇χ)H - (1/2)|H|²∇χ   *** RECOMMENDED ***
//         This is the ONLY formula that produces real magnetic deformation.
//         The (H·∇χ)H term is NOT a gradient and cannot be absorbed by pressure.
//   1 = Interfacial only: -(1/2)|H|²∇χ
//         Equivalent to Paper (differs by a pressure-absorbed gradient).
//         Does NOT produce magnetic deformation in incompressible flow.
//   2 = Paper Eq.(8): (χ/2)∇|H|²
//         Equivalent to Interfacial (Paper = Interfacial + (1/2)∇(χ|H|²)).
//         Does NOT produce magnetic deformation in incompressible flow.
//
// IMPORTANT: Parallel testing at Bom=1.94 (t=6000) showed:
//   Maxwell     -> elong=1.0521 (real magnetic deformation)
//   Paper       -> elong=1.0204 (== Interfacial, == buoyancy-only baseline)
//   Interfacial -> elong=1.0204 (== Paper, no magnetic deformation)
// See CHANGELOG.md "修改 11" for the mathematical proof and full analysis.
// Can be overridden at compile time: -DFORCE_FORMULA=0/1/2
#ifndef FORCE_FORMULA
#define FORCE_FORMULA 0
#endif

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

// bubble
T Bubble_Radius;
Vector<T, 2> Bubble_Center;

// phase field
T Interface_Width, Mobility, Tau_phi, Omega_phi, Kappa, Beta;

// two-phase
T rho_l, rho_h, eta_l, eta_h, sigma, gravity, Eo, Re, U_g, Tau_ns;

// magnetic field
T chi_l, chi_h, mu_l, mu_h, H0, Bom, DeltaRho;

int MaxStep, OutputStep;
std::string work_dir;

void readParam() {
  iniReader r("bubbleMag2d.ini");
  work_dir = r.getValue<std::string>("workdir","workdir_");
  Thread_Num = r.getValue<int>("parallel","thread_num");
  Ni = r.getValue<int>("Mesh","Ni"); Nj = r.getValue<int>("Mesh","Nj");
  Cell_Len = r.getValue<T>("Mesh","Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh","BlockCellLen");
  Bubble_Radius = r.getValue<T>("Bubble","Radius");
  Bubble_Center[0]=r.getValue<T>("Bubble","CenterX");
  Bubble_Center[1]=r.getValue<T>("Bubble","CenterY");
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

  eta_l=T(0.6)/T(3500.0); eta_h=T(0.6)/T(35.0);
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
    printf("---- Bubble Rising in Ferrofluid ----\n");
    printf("Mesh: %dx%d  BlockCellLen=%d\n",Ni,Nj,BlockCellLen);
    printf("Bubble: R=%.1f center=(%.0f,%.0f)\n",Bubble_Radius,Bubble_Center[0],Bubble_Center[1]);
    printf("rho: l=%.4f h=%.1f  eta: l=%.6f h=%.6f  sigma=%.2e g=%.2e\n",rho_l,rho_h,eta_l,eta_h,sigma,gravity);
    printf("Eo=%.0f Re=%.0f W=%.1f M=%.3f tau_phi=%.3f\n",Eo,Re,Interface_Width,Mobility,Tau_phi);
    printf("Magnetic: chi=(%.1f,%.1f) mu=(%.1f,%.1f) H0=%.3f Bom=%.3f\n",chi_l,chi_h,mu_l,mu_h,H0,Bom);
    printf("U_g=%.3f Ma=%.3f\n",U_g,Ma);
    printf("---------------------------------------\n");
  }
  if(Ma>T(0.2)) MPI_RANK(0){ fprintf(stderr,"[Warn] Ma=%.3f > 0.2\n",Ma); }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag=1,BulkFlag=2,BouncebackFlag=4,PeriodicFlag=8;
  mpi().init(&argc,&argv); MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing Bubble Rising in Ferrofluid..."));
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

  // -- init phi (tanh bubble) --
  T R_phys=Bubble_Radius*Cell_Len, xc=Bubble_Center[0]*Cell_Len, yc=Bubble_Center[1]*Cell_Len;
  T W_phys=Interface_Width*Cell_Len;
  auto& phiField=PFLattice.getField<PHI<T>>();
  for(int b=0;b<Geo.getBlockNum();++b){
    const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
    auto& bPhi=phiField.getBlockField(b); auto& bl=PFLattice.getBlockLat(b);
    T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
    int ov=0;
    for(int j=ov;j<bk.getNy()-ov;++j){
      T y=my+T(j)*vs, dy=y-yc;
      for(int i=ov;i<bk.getNx()-ov;++i){
        T x=mx+T(i)*vs, dx=x-xc;
        T dist=std::sqrt(dx*dx+dy*dy);
        T phi=T{0.5}+T{0.5}*std::tanh(T{2.0}*(dist-R_phys)/W_phys);
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

  // init MF: psi = -H0*y (uniform field), g_k = w_k*psi
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
  vtmo::vtmWriter<T,2> MW("bubbleMag2d",Geo);
  MW.addWriterSet(PW,PS,VW,Dw,Fw,Hxw,Hyw);

  // ===== Pre-converge magnetic field to steady state for diagnostics =====
  {
    MPI_RANK(0) printf("Pre-converging magnetic field (200 iters)...\n"); fflush(stdout);
    for(int iter=0;iter<200;++iter){
      // 0a: Update coeffs
      MCC.ApplyInnerCellDynamics<MCSel>(0,FlagFM);
      CommunicateOMEGAPSI<T>(MFLattice);
      // 0b: Set wall psi = -H0*y
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b);
          int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
          T minY=bk.getMin()[1],maxY=bk.getMax()[1],vs=bk.getVoxelSize();
          T H_global=T(Nj)*Cell_Len;
          if(minY<Cell_Len*T{1.5}){
            for(int jj=0;jj<=ov;++jj){
              T y=minY+T(jj)*vs, psi_w=-H0*y;
              for(int ii=0;ii<nx;++ii){
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
          if(maxY>H_global-Cell_Len*T{1.5}){
            for(int jj=ny-ov-1;jj<ny;++jj){
              T y=minY+T(jj)*vs, psi_w=-H0*y;
              for(int ii=0;ii<nx;++ii){
                std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
                for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
                bPsi.get(id)=psi_w;
              }
            }
          }
        }
      }
      // 0c: Collision
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          collision::MRTDiffusion<MFCELL,OMEGA_PSI<T>>::apply(c);
        }
      }
      MF_Per.Apply();
      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();
      // 0e: PSI
      {
        auto& psiF=MFLattice.getField<PSI<T>>();
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          auto& bPsi=psiF.getBlockField(b); int ov=bk.getOverlap();
          for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
            std::size_t id=j*pr[1]+i; MFCELL c(id,bl);
            T psi=0; for(unsigned k=0;k<MFLatSet::q;++k)psi+=c[k];
            bPsi.get(id)=psi;
          }
        }
      }
      CommunicatePSI<T>(MFLattice);
      // 0f: Compute H
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
          MFCELL c(j*pr[1]+i,bl);
          MFComputeH2D<MFCELL>::apply(c);
        }
      }
      CommunicateAllMFFields<T>(MFLattice);
    }
    MPI_RANK(0) printf("Pre-convergence done.\n");
  }

  // Force diagnostic after pre-convergence
  {
    T maxF_mag=0, maxF_inter=0, maxF_body=0;
    T maxHx=0, maxHy=0, maxDchi=0, maxDHsq=0, maxChi=0;
    T maxProd=0;
    int nz_body=0, nz_inter=0, total_interface=0;
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl);
        T chi_l=mf.template get<CHI_L<T>>(), chi_h=mf.template get<CHI_H<T>>();
        T chi=mf.template get<CHI_PERCELL<T>>();
        T Hx=mf.template get<HX<T>>(), Hy=mf.template get<HY<T>>(), Hmag=mf.template get<HMAG<T>>();
        T dchi_x=0,dchi_y=0;
        for(unsigned k=1;k<MFLatSet::q;++k){
          T pk=mf.getNeighbor(k).template get<PHI<T>>();
          T wk=latset::w<MFLatSet>(k); const auto& ck=latset::c<MFLatSet>(k);
          dchi_x+=wk*ck[0]*pk; dchi_y+=wk*ck[1]*pk;
        }
        dchi_x=(chi_h-chi_l)*dchi_x/MFLatSet::cs2;
        dchi_y=(chi_h-chi_l)*dchi_y/MFLatSet::cs2;
        T dHsq_x=0,dHsq_y=0;
        for(unsigned k=1;k<MFLatSet::q;++k){
          T Hk=mf.getNeighbor(k).template get<HMAG<T>>();
          T Hsqk=Hk*Hk; T wk=latset::w<MFLatSet>(k); const auto& ck=latset::c<MFLatSet>(k);
          dHsq_x+=wk*ck[0]*Hsqk; dHsq_y+=wk*ck[1]*Hsqk;
        }
        dHsq_x/=MFLatSet::cs2; dHsq_y/=MFLatSet::cs2;
        T HdotDchi=Hx*dchi_x+Hy*dchi_y;
        T half_Hsq=T{0.5}*Hmag*Hmag;
        T Fmx=Hx*HdotDchi-half_Hsq*dchi_x, Fmy=Hy*HdotDchi-half_Hsq*dchi_y;
        T Fmag=std::sqrt(Fmx*Fmx+Fmy*Fmy);
        T Fix=-half_Hsq*dchi_x, Fiy=-half_Hsq*dchi_y;
        T Fi=std::sqrt(Fix*Fix+Fiy*Fiy);
        T half_chi=T{0.5}*chi;
        T Fbx=half_chi*dHsq_x, Fby=half_chi*dHsq_y;
        T Fb=std::sqrt(Fbx*Fbx+Fby*Fby);
        if(Fmag>maxF_mag)maxF_mag=Fmag;
        if(Fi>maxF_inter)maxF_inter=Fi;
        if(Fb>maxF_body)maxF_body=Fb;
        T hx=std::abs(Hx); if(hx>maxHx)maxHx=hx;
        T hy=std::abs(Hy); if(hy>maxHy)maxHy=hy;
        T dc=std::sqrt(dchi_x*dchi_x+dchi_y*dchi_y); if(dc>maxDchi)maxDchi=dc;
        T dh=std::sqrt(dHsq_x*dHsq_x+dHsq_y*dHsq_y); if(dh>maxDHsq)maxDHsq=dh;
        if(chi>maxChi)maxChi=chi;
        T prod=chi*dh; if(prod>maxProd)maxProd=prod;
        if(dc>1e-4) total_interface++;
        if(Fb>1e-15) nz_body++;
        if(Fi>1e-15) nz_inter++;
      }
    }
    MPI_RANK(0){
      printf("[Post-converge] |H|_max=(%.2e,%.2e) |∇χ|_max=%.2e |∇|H|²|_max=%.2e χ_max=%.2e\n",
             maxHy,maxHx,maxDchi,maxDHsq,maxChi);
      printf("  F_mag(Maxwell)=%.2e  F_inter=%.2e  F_body(paper)=%.2e\n",
             maxF_mag,maxF_inter,maxF_body);
      printf("  max(χ*|∇|H|²|)=%.2e  interface_cells=%d  nz_body=%d  nz_inter=%d\n",
             maxProd,total_interface,nz_body,nz_inter);
      fflush(stdout);
    }
  }

  // ===== Direct NS force comparison: apply each formula and measure total force =====
  {
    MPI_RANK(0) printf("Direct NS force comparison (Bom=%.3f)...\n", Bom); fflush(stdout);
    // Test 1: Full Maxwell
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
      auto& ns_bl=NSLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
        MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
      }
    }
    T sumF_maxwell=0, maxF_maxwell=0;
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& ns_bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL ns(id,ns_bl);
        auto F=ns.template get<FORCE<T,2>>();
        T f=std::sqrt(F[0]*F[0]+F[1]*F[1]);
        sumF_maxwell+=f; if(f>maxF_maxwell)maxF_maxwell=f;
      }
    }

    // Test 2: Paper formula
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
      auto& ns_bl=NSLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
        MFMagneticForcePaper2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
      }
    }
    T sumF_paper=0, maxF_paper=0;
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& ns_bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL ns(id,ns_bl);
        auto F=ns.template get<FORCE<T,2>>();
        T f=std::sqrt(F[0]*F[0]+F[1]*F[1]);
        sumF_paper+=f; if(f>maxF_paper)maxF_paper=f;
      }
    }

    // Test 3: Interfacial only
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
      auto& ns_bl=NSLattice.getBlockLat(b);
      const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i;
        PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl); NSCELL ns(id,ns_bl);
        MFMagneticForceInterfacial2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
      }
    }
    T sumF_inter=0, maxF_inter=0;
    for(int b=0;b<Geo.getBlockNum();++b){
      auto& ns_bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
      int ov=bk.getOverlap();
      for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
        std::size_t id=j*pr[1]+i; NSCELL ns(id,ns_bl);
        auto F=ns.template get<FORCE<T,2>>();
        T f=std::sqrt(F[0]*F[0]+F[1]*F[1]);
        sumF_inter+=f; if(f>maxF_inter)maxF_inter=f;
      }
    }
    MPI_RANK(0){
      printf("NS Force comparison:\n");
      printf("  Maxwell:   sum|F|=%.2e  max|F|=%.2e\n", sumF_maxwell, maxF_maxwell);
      printf("  Interfacial: sum|F|=%.2e  max|F|=%.2e\n", sumF_inter, maxF_inter);
      printf("  Paper:     sum|F|=%.2e  max|F|=%.2e\n", sumF_paper, maxF_paper);
      printf("  Ratio paper/Maxwell sum=%.2f  max=%.2f\n",
             sumF_paper/maxF_paper, maxF_paper/maxF_maxwell);
      fflush(stdout);
    }
    // Reset force to zero for clean start
    NSLattice.getField<FORCE<T,2>>().InitValue(Vector<T,2>{0,0});
  }

  // Write initial state
  MW.WriteBinary(0);

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  Timer t; Timer ot;
  T H_global=T(Nj)*Cell_Len;

  while(t()<MaxStep){
    // ===== Phase 0: Magnetic field solve =====
    // 0a: Update per-cell mu, chi, omega_psi from phi
    MCC.ApplyInnerCellDynamics<MCSel>(t(),FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);

    // 0b: Set wall psi = -H0*y  and wall pops = feq(psi_bc)
    {
      auto& psiF=MFLattice.getField<PSI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& bl=MFLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bPsi=psiF.getBlockField(b);
        int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T minY=bk.getMin()[1],maxY=bk.getMax()[1],vs=bk.getVoxelSize();
        if(minY<Cell_Len*T{1.5}){
          for(int jj=0;jj<=ov;++jj){
            T y=minY+T(jj)*vs, psi_w=-H0*y;
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
              bPsi.get(id)=psi_w;
            }
          }
        }
        if(maxY>H_global-Cell_Len*T{1.5}){
          for(int jj=ny-ov-1;jj<ny;++jj){
            T y=minY+T(jj)*vs, psi_w=-H0*y;
            for(int ii=0;ii<nx;++ii){
              std::size_t id=jj*pr[1]+ii; MFCELL c(id,bl);
              for(unsigned k=0;k<MFLatSet::q;++k) c[k]=latset::w<MFLatSet>(k)*psi_w;
              bPsi.get(id)=psi_w;
            }
          }
        }
      }
    }

    // 0c: MF collision
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
#if FORCE_FORMULA == 0
            // Full Maxwell: (H·∇χ)H - (1/2)|H|²∇χ
            MFMagneticForce2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#elif FORCE_FORMULA == 1
            // Interfacial only: -(1/2)|H|²∇χ
            MFMagneticForceInterfacial2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#else
            // Paper Eq.(8): (χ/2)∇|H|²
            MFMagneticForcePaper2D<PFCELL,MFCELL,NSCELL>::apply(pf,mf,ns);
#endif
          }
        }
      }
    }

    // Force diagnostic at step 0 and every 10000 steps
    if(t()==0||t()%10000==0){
      T maxF_mag=0, maxF_inter=0, maxF_body=0;
      T maxHx=0, maxHy=0, maxDchi=0, maxDHsq=0, maxChi=0;
      // Also track where the max values occur and the product at interface
      T maxProd=0; // max of chi * |∇|H|²|
      int nz_body=0, nz_inter=0, total_interface=0;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& pf_bl=PFLattice.getBlockLat(b); auto& mf_bl=MFLattice.getBlockLat(b);
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
          std::size_t id=j*pr[1]+i;
          PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl);
          T chi_l=mf.template get<CHI_L<T>>(), chi_h=mf.template get<CHI_H<T>>();
          T chi=mf.template get<CHI_PERCELL<T>>();
          T Hx=mf.template get<HX<T>>(), Hy=mf.template get<HY<T>>(), Hmag=mf.template get<HMAG<T>>();
          // ∇χ from φ
          T dchi_x=0,dchi_y=0;
          for(unsigned k=1;k<MFLatSet::q;++k){
            T pk=mf.getNeighbor(k).template get<PHI<T>>();
            T wk=latset::w<MFLatSet>(k); const auto& ck=latset::c<MFLatSet>(k);
            dchi_x+=wk*ck[0]*pk; dchi_y+=wk*ck[1]*pk;
          }
          dchi_x=(chi_h-chi_l)*dchi_x/MFLatSet::cs2;
          dchi_y=(chi_h-chi_l)*dchi_y/MFLatSet::cs2;
          // ∇|H|²
          T dHsq_x=0,dHsq_y=0;
          for(unsigned k=1;k<MFLatSet::q;++k){
            T Hk=mf.getNeighbor(k).template get<HMAG<T>>();
            T Hsqk=Hk*Hk; T wk=latset::w<MFLatSet>(k); const auto& ck=latset::c<MFLatSet>(k);
            dHsq_x+=wk*ck[0]*Hsqk; dHsq_y+=wk*ck[1]*Hsqk;
          }
          dHsq_x/=MFLatSet::cs2; dHsq_y/=MFLatSet::cs2;
          T HdotDchi=Hx*dchi_x+Hy*dchi_y;
          T half_Hsq=T{0.5}*Hmag*Hmag;
          // Maxwell: F_m = (H·∇χ)H - (1/2)|H|²∇χ
          T Fmx=Hx*HdotDchi-half_Hsq*dchi_x, Fmy=Hy*HdotDchi-half_Hsq*dchi_y;
          T Fmag=std::sqrt(Fmx*Fmx+Fmy*Fmy);
          // Interfacial only: F_m = -(1/2)|H|²∇χ
          T Fix=-half_Hsq*dchi_x, Fiy=-half_Hsq*dchi_y;
          T Fi=std::sqrt(Fix*Fix+Fiy*Fiy);
          // Body: F_m = (χ/2)∇|H|²
          T half_chi=T{0.5}*chi;
          T Fbx=half_chi*dHsq_x, Fby=half_chi*dHsq_y;
          T Fb=std::sqrt(Fbx*Fbx+Fby*Fby);
          if(Fmag>maxF_mag)maxF_mag=Fmag;
          if(Fi>maxF_inter)maxF_inter=Fi;
          if(Fb>maxF_body)maxF_body=Fb;
          T hx=std::abs(Hx); if(hx>maxHx)maxHx=hx;
          T hy=std::abs(Hy); if(hy>maxHy)maxHy=hy;
          T dc=std::sqrt(dchi_x*dchi_x+dchi_y*dchi_y); if(dc>maxDchi)maxDchi=dc;
          T dh=std::sqrt(dHsq_x*dHsq_x+dHsq_y*dHsq_y); if(dh>maxDHsq)maxDHsq=dh;
          if(chi>maxChi)maxChi=chi;
          T prod=chi*dh; if(prod>maxProd)maxProd=prod;
          if(dc>1e-4) total_interface++;
          if(Fb>1e-15) nz_body++;
          if(Fi>1e-15) nz_inter++;
        }
      }
      MPI_RANK(0){
        printf("[Diag t=%d] |H|_max=(%.2e,%.2e) |∇χ|_max=%.2e |∇|H|²|_max=%.2e χ_max=%.2e\n",
               t(),maxHy,maxHx,maxDchi,maxDHsq,maxChi);
        printf("  F_mag(Maxwell)=%.2e  F_inter=%.2e  F_body(paper)=%.2e\n",
               maxF_mag,maxF_inter,maxF_body);
        printf("  max(χ*|∇|H|²|)=%.2e  interface_cells=%d  nz_body=%d  nz_inter=%d\n",
               maxProd,total_interface,nz_body,nz_inter);
        fflush(stdout);
      }
    }

    // --- Direct NS force check at interface points ---
    if(t()==0){
      T nsF_at_top=0, nsF_at_side=0, nsF_at_bottom=0, nsF_at_center=0;
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& ns_bl=NSLattice.getBlockLat(b);
        T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
        T cx=T(128)*Cell_Len, cy=T(128)*Cell_Len, R=T(50)*Cell_Len;
        for(int j=0;j<bk.getNy();++j) for(int i=0;i<bk.getNx();++i){
          T x=mx+T(i)*vs, y=my+T(j)*vs;
          T dx=x-cx, dy=y-cy, dist=std::sqrt(dx*dx+dy*dy);
          std::size_t id=j*pr[1]+i;
          // Top pole: (128, 78) = center + (0, -R)
          if(std::abs(x-cx)<vs && std::abs(y-(cy-R))<vs){
            NSCELL ns(id,ns_bl); auto F=ns.template get<FORCE<T,2>>();
            nsF_at_top=std::sqrt(F[0]*F[0]+F[1]*F[1]); }
          // Side: (178, 128) = center + (R, 0)
          if(std::abs(x-(cx+R))<vs && std::abs(y-cy)<vs){
            NSCELL ns(id,ns_bl); auto F=ns.template get<FORCE<T,2>>();
            nsF_at_side=std::sqrt(F[0]*F[0]+F[1]*F[1]); }
          // Bottom: (128, 178) = center + (0, +R)
          if(std::abs(x-cx)<vs && std::abs(y-(cy+R))<vs){
            NSCELL ns(id,ns_bl); auto F=ns.template get<FORCE<T,2>>();
            nsF_at_bottom=std::sqrt(F[0]*F[0]+F[1]*F[1]); }
          // Center: (128, 128)
          if(std::abs(x-cx)<vs && std::abs(y-cy)<vs){
            NSCELL ns(id,ns_bl); auto F=ns.template get<FORCE<T,2>>();
            nsF_at_center=std::sqrt(F[0]*F[0]+F[1]*F[1]); }
        }
      }
      MPI_RANK(0) printf("[NS Force] top=%.2e  side=%.2e  bottom=%.2e  center=%.2e\n",
             nsF_at_top, nsF_at_side, nsF_at_bottom, nsF_at_center);
      fflush(stdout);
    }

    GrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,GrT>>(t(),FlagFM);
    PrC.ApplyInnerCellDynamics<CoupledTaskSelector<std::uint8_t,PFCELL,NSCELL,PrT>>(t(),FlagFM);
    NSLattice.getField<FORCE<T,2>>().Communicate();

    // === NS force sum diagnostic after ALL forces applied (before NS collision) ===
    if(t()%200==0){
      T totalNSForce=0, maxNSForce=0;
      T sumFx=0, sumFy=0;
      int nonZero=0;
      for(int b=0;b<Geo.getBlockNum();++b){
        auto& ns_bl=NSLattice.getBlockLat(b); const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        int ov=bk.getOverlap();
        for(int j=ov;j<bk.getNy()-ov;++j) for(int i=ov;i<bk.getNx()-ov;++i){
          std::size_t id=j*pr[1]+i; NSCELL ns(id,ns_bl);
          auto F=ns.template get<FORCE<T,2>>();
          T f=std::sqrt(F[0]*F[0]+F[1]*F[1]);
          totalNSForce+=f; if(f>maxNSForce)maxNSForce=f;
          sumFx+=F[0]; sumFy+=F[1];
          if(f>1e-15) nonZero++;
        }
      }
      MPI_RANK(0){
        printf("[NS Force t=%d] sum|F|=%.2e max|F|=%.2e sumF=(%.2e,%.2e) nz=%d\n",
               t(), totalNSForce, maxNSForce, sumFx, sumFy, nonZero);
        fflush(stdout);
      }
    }

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

    // wall phi=1 BC
    {
      auto& pF=PFLattice.getField<PHI<T>>();
      for(int b=0;b<Geo.getBlockNum();++b){
        const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
        auto& bP=pF.getBlockField(b); int nx=bk.getNx(),ny=bk.getNy(),ov=bk.getOverlap();
        T my=bk.getMin()[1],My=bk.getMax()[1];
        if(my<Cell_Len*T{1.5}) for(int i=0;i<nx;++i)bP.get(ov*pr[1]+i)=T{1};
        if(My>H_global-Cell_Len*T{1.5}){int jj=ny-1-ov;for(int i=0;i<nx;++i)bP.get(jj*pr[1]+i)=T{1};}
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
        if(my<Cell_Len*T{1.5}) for(int j=0;j<ov;++j) for(int i=0;i<nx;++i) bP.get(j*pr[1]+i)=T{1};
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
      // ===== Geometry & Field Enhancement Diagnostic (Phase 1) =====
      // Measures: elongation ratio (Dy/Dx), centroid, area, field enhancement Henh
      // Henh = avg|H| inside bubble / H0, should be ~1.33 if field correctly enhanced
      {
        T sumW=T{0}, sumX=T{0}, sumY=T{0};
        T sumH_in=T{0}; int n_in=0;
        T xmin=1e30, xmax=-1e30, ymin=1e30, ymax=-1e30;
        for(int b=0;b<Geo.getBlockNum();++b){
          auto& pf_bl=PFLattice.getBlockLat(b);
          auto& mf_bl=MFLattice.getBlockLat(b);
          const auto& bk=Geo.getBlock(b); const auto& pr=bk.getProjection();
          int ov=bk.getOverlap();
          T vs=bk.getVoxelSize(),mx=bk.getMin()[0],my=bk.getMin()[1];
          for(int j=ov;j<bk.getNy()-ov;++j){
            for(int i=ov;i<bk.getNx()-ov;++i){
              std::size_t id=j*pr[1]+i;
              PFCELL pf(id,pf_bl); MFCELL mf(id,mf_bl);
              T phi=pf.template get<PHI<T>>();
              T x=mx+T(i)*vs, y=my+T(j)*vs;
              if(phi < T{0.5}){ // inside bubble
                T w=T{1}-phi;
                sumW+=w; sumX+=x*w; sumY+=y*w;
                if(x<xmin)xmin=x; if(x>xmax)xmax=x;
                if(y<ymin)ymin=y; if(y>ymax)ymax=y;
              }
              if(phi < T{0.1}){ // deep inside (for field enhancement)
                T Hmag=mf.template get<HMAG<T>>();
                sumH_in+=Hmag; n_in++;
              }
            }
          }
        }
        T xc = (sumW>T{0}) ? sumX/sumW : T{0};
        T yc = (sumW>T{0}) ? sumY/sumW : T{0};
        T Dx = xmax-xmin, Dy = ymax-ymin;
        T elong = (Dx>T{0}) ? Dy/Dx : T{0};
        T Henh = (n_in>0 && H0>T{0}) ? (sumH_in/T(n_in))/H0 : T{0};
        MPI_RANK(0){
          printf("[Geom t=%d] elong=%.4f center=(%.2f,%.2f) area=%.2f Dx=%.2f Dy=%.2f Henh=%.4f (expect~1.33)\n",
                 t(), elong, xc, yc, sumW, Dx, Dy, Henh);
          fflush(stdout);
        }
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
