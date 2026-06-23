// ferro_bubble.cpp
// Ferrofluid bubble rise — PF+NS+Mag three-solver FDLBM with MPI
// Pattern: mirrors addMagField branch's ferrodrop2d.cpp (single FlagFM)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)

#include "freelb.h"
#include "freelb.hh"

using T = FLOAT;
using PFLatSet  = D2Q9<T>;
using NSLatSet  = D2Q9<T>;
using MagLatSet = D2Q9<T>;

/*----------------------------------------------
            Simulation Parameters
-----------------------------------------------*/
int Ni, Nj, BlockCellLen, Thread_Num;
T Cell_Len;
T Bubble_Radius;
Vector<T, 2> Bubble_Center;

T Interface_Width, Mobility;
T Tau_phi, Kappa, Beta;

T rho_l, rho_h, eta_l, eta_h, mu_l, mu_h, mu0, H0;
T sigma, Eo, Re, U_g, gravity;

int MaxStep, OutputStep;

void readParam() {
  iniReader reader("ferro_bubble.ini");
  Thread_Num = reader.getValue<int>("parallel", "thread_num");
  Ni = reader.getValue<int>("Mesh", "Ni");
  Nj = reader.getValue<int>("Mesh", "Nj");
  Cell_Len = reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = reader.getValue<int>("Mesh", "BlockCellLen");
  Bubble_Radius = reader.getValue<T>("Bubble", "Radius");
  Bubble_Center[0] = reader.getValue<T>("Bubble", "CenterX");
  Bubble_Center[1] = reader.getValue<T>("Bubble", "CenterY");
  rho_l = reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = reader.getValue<T>("Two_Phase", "rho_h");
  eta_l = reader.getValue<T>("Two_Phase", "eta_l");
  eta_h = reader.getValue<T>("Two_Phase", "eta_h");
  mu_l = reader.getValue<T>("Two_Phase", "mu_l");
  mu_h = reader.getValue<T>("Two_Phase", "mu_h");
  mu0  = reader.getValue<T>("Two_Phase", "mu0");
  Eo = reader.getValue<T>("Two_Phase", "Eo");
  Re = reader.getValue<T>("Two_Phase", "Re");
  Interface_Width = reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = reader.getValue<T>("Phase_Field", "Mobility");
  sigma = reader.getValue<T>("Phase_Field", "sigma");
  MaxStep = reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = reader.getValue<int>("Simulation_Settings", "OutputStep");

  // 5-step parameter design
  T D = T(2.0) * Bubble_Radius;
  U_g = T(0.04);
  T g_abs = U_g * U_g / D;
  gravity = -g_abs;
  T nu = U_g * D / Re;
  eta_h = nu * rho_h;
  eta_l = eta_h / T(10);
  T DeltaRho = rho_h - rho_l;
  sigma = DeltaRho * g_abs * D * D / Eo;

  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  Beta  = T(12.0) * sigma / Interface_Width;
  Tau_phi = T(0.5) + Mobility / PFLatSet::cs2;

  MPI_RANK(0) {
    std::cout << "------ Ferrofluid Bubble Rise FDLBM ------\n";
    std::cout << "[Mesh]:" << Ni << "x" << Nj
              << " BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]:R=" << Bubble_Radius
              << " Ctr=(" << Bubble_Center[0] << "," << Bubble_Center[1] << ")\n";
    std::cout << "[TwoPhase]:rho_l=" << rho_l << " rho_h=" << rho_h
              << " eta_l=" << eta_l << " eta_h=" << eta_h << "\n";
    std::cout << "[TwoPhase]:mu_l=" << mu_l << " mu_h=" << mu_h
              << " Eo=" << Eo << " Re=" << Re << "\n";
    std::cout << "[Phase]:W=" << Interface_Width << " M=" << Mobility
              << " beta=" << Beta << " kappa=" << Kappa << "\n";
    std::cout << "[Sim]:MaxStep=" << MaxStep
              << " OutputStep=" << OutputStep << "\n";
#ifdef MPI_ENABLED
    std::cout << "[Parallel]:" << mpi().getSize() << " MPI processes\n";
#endif
    std::cout << "------------------------------------------\n";
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner("Initializing Ferrofluid Bubble Rise FDLBM...");
  readParam();

  // ---------- converters ----------
  BaseConverter<T> BaseConv(PFLatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), T(0.55));

  BaseConverter<T> PFBaseConv(PFLatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  T lat_mu0 = T(0.01);
  T mu_avg = (mu_l + mu_h) * T(0.5);
  T mag_tau = T(0.5) + mu_avg / PFLatSet::cs2;
  if (mag_tau < T(0.55)) mag_tau = T(0.55);
  BaseConverter<T> MagBaseConv(PFLatSet::cs2);
  MagBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), mag_tau);

  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ---------- geometry ----------
  AABB<T, 2> domain(Vector<T, 2>(T{0}, T{0}),
                     Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ---------- single shared flag field ----------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<PFLatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Ferro", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ==============================================================
  // NS lattice
  // ==============================================================
  using NSFIELDS = TypePack<
      RHO<T>, VELOCITY<T, 2>, POP<T, 9>, FORCE<T, 2>,
      PRESSURE<T>, VISCOSITY<T>,
      RHOL<T>, RHOH<T>, ETAL<T>, ETAH<T>,
      CONSTFORCE<T, 2>>;
  ValuePack NSInit(
      BaseConv.getLatRhoInit(),
      Vector<T, 2>{T{0}, T{0}},
      T{},
      Vector<T, 2>{T{0}, T{0}},
      T{},
      eta_h,
      rho_l, rho_h, eta_l, eta_h,
      Vector<T, 2>{T{0}, gravity});
  using NSCELL = Cell<T, NSLatSet, NSFIELDS>;
  BlockLatticeManager<T, NSLatSet, NSFIELDS> NSLattice(Geo, NSInit, BaseConv);

  // ==============================================================
  // PF lattice
  // ==============================================================
  using PFFIELDS = TypePack<
      PHI<T>, POP<T, 9>, GRAD<T, 2>, PHASENORMAL<T, 2>,
      INTERFACEWIDTH<T>, PHASELAPLACIAN<T>, PHASECHEMPOTENTIAL<T>,
      GBETA<T>, TAUPHI<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, 2>>;
  using PFFPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInit(
      T{},
      T{},
      Vector<T, 2>{},
      Vector<T, 2>{},
      Interface_Width,
      T{}, T{},
      Beta, Tau_phi);
  using PFCELL = Cell<T, PFLatSet, ExtractFieldPack<PFFPACK>::mergedpack>;
  BlockLatticeManager<T, PFLatSet, PFFPACK> PFLattice(
      Geo, PFInit, PFBaseConv,
      &NSLattice.getField<VELOCITY<T, 2>>());

  // ==============================================================
  // Mag lattice
  // ==============================================================
  using MAGFIELDS = TypePack<
      PSI<T>, POP<T, 9>, MAGH<T, 2>, MAGHSQ<T>,
      MAGPERMEABILITY<T>, FORCE<T, 2>,
      MUL<T>, MUH<T>, FLAG>;
  using MAGFIELDREFS = TypePack<PHI<T>>;
  using MAGFPACK = TypePack<MAGFIELDS, MAGFIELDREFS>;
  ValuePack MagInit(
      T{}, T{}, Vector<T, 2>{},
      T{}, T{}, Vector<T, 2>{},
      mu_l, mu_h, BulkFlag);
  using MAGCELL = Cell<T, MagLatSet, ExtractFieldPack<MAGFPACK>::mergedpack>;
  BlockLatticeManager<T, MagLatSet, MAGFPACK> MagLattice(
      Geo, MagInit, MagBaseConv,
      &PFLattice.getField<PHI<T>>());

  // ==============================================================
  // Initialize PHI (tanh bubble profile)
  // ==============================================================
  T R_phys = Bubble_Radius * Cell_Len;
  T xc = Bubble_Center[0] * Cell_Len;
  T yc = Bubble_Center[1] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    const auto& blk = Geo.getBlock(b);
    auto& bphi = phiField.getBlockField(b);
    T dx = blk.getVoxelSize();
    for (int j = 0; j < blk.getNy(); ++j) {
      for (int i = 0; i < blk.getNx(); ++i) {
        std::size_t id = j * blk.getProjection()[1] + i;
        T x = blk.getMin()[0] + T(i) * dx;
        T y = blk.getMin()[1] + T(j) * dx;
        T dist = std::sqrt((x - xc) * (x - xc) + (y - yc) * (y - yc));
        bphi.get(id) = T(0.5) + T(0.5)
          * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
      }
    }
  }

  // ==============================================================
  // Initialize PF pops (zero-velocity equilibrium)
  // ==============================================================
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& lat = PFLattice.getBlockLat(b);
    auto& bphi = phiField.getBlockField(b);
    const auto& blk = Geo.getBlock(b);
    for (int j = 0; j < blk.getNy(); ++j) {
      for (int i = 0; i < blk.getNx(); ++i) {
        std::size_t id = j * blk.getProjection()[1] + i;
        PFCELL cell(id, lat);
        for (int k = 0; k < 9; ++k) {
          cell[k] = latset::w<PFLatSet>(k) * bphi.get(id);
        }
      }
    }
  }

  // ==============================================================
  // Initialize NS pops (quiescent equilibrium for velocity-based NS)
  // For velocity-based NS (Eq.31), zero-moment = p/(rho*cs^2).
  // At rest with p = rho*cs^2, g_eq[k] = w_k. So init g[k] = w_k.
  // Also init PRESSURE = rho*cs^2 so NSEquilibrium has valid input.
  // ==============================================================
  T rho_init = BaseConv.getLatRhoInit();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& lat = NSLattice.getBlockLat(b);
    const auto& blk = Geo.getBlock(b);
    for (int j = 0; j < blk.getNy(); ++j) {
      for (int i = 0; i < blk.getNx(); ++i) {
        std::size_t id = j * blk.getProjection()[1] + i;
        NSCELL cell(id, lat);
        for (int k = 0; k < 9; ++k) {
          cell[k] = latset::w<NSLatSet>(k);
        }
        cell.template get<PRESSURE<T>>() = rho_init * NSLatSet::cs2;
      }
    }
  }

  // ==============================================================
  // Initialize Mag psi = -H0*y, pops = w*psi
  // ==============================================================
  T H0_lat = T{1};
  auto& psiField = MagLattice.getField<PSI<T>>();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& lat = MagLattice.getBlockLat(b);
    auto& bpsi = psiField.getBlockField(b);
    const auto& blk = Geo.getBlock(b);
    T dx = blk.getVoxelSize();
    for (int j = 0; j < blk.getNy(); ++j) {
      for (int i = 0; i < blk.getNx(); ++i) {
        std::size_t id = j * blk.getProjection()[1] + i;
        T y = blk.getMin()[1] + T(j) * dx;
        T psi = -H0_lat * y;
        bpsi.get(id) = psi;
        MAGCELL cell(id, lat);
        for (int k = 0; k < 9; ++k) {
          cell[k] = latset::w<MagLatSet>(k) * psi;
        }
      }
    }
  }

  // ==============================================================
  // Set field constants
  // ==============================================================
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);
  PFLattice.getField<GBETA<T>>().InitValue(Beta);

  // ==============================================================
  // Boundary conditions (single shared FlagFM)
  // ==============================================================
  using PFBLK  = BlockLatticeManager<T, PFLatSet,  PFFPACK>;
  using NSBLK  = BlockLatticeManager<T, NSLatSet,  NSFIELDS>;
  using MAGBLK = BlockLatticeManager<T, MagLatSet, MAGFPACK>;

  BBLikeFixedBlockBdManager<
      bounceback::normal<PFCELL>, PFBLK, BlockFieldManager<FLAG, T, 2>>
      PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeFixedBlockBdManager<
      bounceback::normal<NSCELL>, NSBLK, BlockFieldManager<FLAG, T, 2>>
      NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeFixedBlockBdManager<
      bounceback::normal<MAGCELL>, MAGBLK, BlockFieldManager<FLAG, T, 2>>
      Mag_BB("Mag_BB", MagLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ==============================================================
  // Task definitions
  // ==============================================================

  // --- PF pre-coupling tasks ---
  using PFGradT = tmp::Key_TypePair<BulkFlag, IsotropicGradient<PFCELL>>;
  using PFNrmT  = tmp::Key_TypePair<BulkFlag, ComputeNormal<PFCELL>>;
  using PFLapT  = tmp::Key_TypePair<BulkFlag, IsotropicLaplacian<PFCELL>>;
  using PFChemT = tmp::Key_TypePair<BulkFlag, ChemicalPotential<PFCELL>>;

  using PFGradS = TaskSelector<std::uint8_t, PFCELL, PFGradT>;
  using PFNrmS  = TaskSelector<std::uint8_t, PFCELL, PFNrmT>;
  using PFLapS  = TaskSelector<std::uint8_t, PFCELL, PFLapT>;
  using PFChemS = TaskSelector<std::uint8_t, PFCELL, PFChemT>;

  // --- PF collision task ---
  using PFCollT = tmp::Key_TypePair<BulkFlag, PhaseFieldMRT<PFCELL>>;
  using PFCollS = TaskSelector<std::uint8_t, PFCELL, PFCollT>;

  // --- NS moment task (VELOCITY/PRESSURE, runs BEFORE collision) ---
  using NSMomT  = tmp::Key_TypePair<BulkFlag, moment::NSMomentum<NSCELL>>;
  using NSMomS  = TaskSelector<std::uint8_t, NSCELL, NSMomT>;
  // --- NS collision task ---
  using NSCollT = tmp::Key_TypePair<BulkFlag, NSMRT<NSCELL>>;
  using NSCollS = TaskSelector<std::uint8_t, NSCELL, NSCollT>;

  // --- Mag tasks ---
  using MagBCT  = tmp::Key_TypePair<BulkFlag, MagNeumannBC<MAGCELL>>;
  using MagMRTT = tmp::Key_TypePair<BulkFlag, MagMRT<MAGCELL>>;
  using MagGrdT = tmp::Key_TypePair<BulkFlag, MagGradient<MAGCELL>>;
  using MagHSqT = tmp::Key_TypePair<BulkFlag, MagHSq<MAGCELL>>;
  using MagFrcT = tmp::Key_TypePair<BulkFlag, MagForce<MAGCELL>>;

  using MagBCS  = TaskSelector<std::uint8_t, MAGCELL, MagBCT>;
  using MagMRTS = TaskSelector<std::uint8_t, MAGCELL, MagMRTT>;
  using MagGrdS = TaskSelector<std::uint8_t, MAGCELL, MagGrdT>;
  using MagHSqS = TaskSelector<std::uint8_t, MAGCELL, MagHSqT>;
  using MagFrcS = TaskSelector<std::uint8_t, MAGCELL, MagFrcT>;

  // ==============================================================
  // Coupling tasks (BlockLatManagerCoupling)
  // ==============================================================
  using PropT = tmp::Key_TypePair<BulkFlag,
      PFtoNS_Properties<PFCELL, NSCELL>>;
  using PropTS = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PropT>;
  BlockLatManagerCoupling PropC(PFLattice, NSLattice);

  using STT = tmp::Key_TypePair<BulkFlag,
      PFtoNS_Forces<PFCELL, NSCELL>>;
  using STTS = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STT>;
  BlockLatManagerCoupling STC(PFLattice, NSLattice);

  using MagPT = tmp::Key_TypePair<BulkFlag,
      PFtoMag_Permeability<PFCELL, MAGCELL>>;
  using MagPTS = CoupledTaskSelector<std::uint8_t, PFCELL, MAGCELL, MagPT>;
  BlockLatManagerCoupling MagPC(PFLattice, MagLattice);

  using MagNT = tmp::Key_TypePair<BulkFlag,
      MagtoNS_Force<MAGCELL, NSCELL>>;
  using MagNTS = CoupledTaskSelector<std::uint8_t, MAGCELL, NSCELL, MagNT>;
  BlockLatManagerCoupling MagNC(MagLattice, NSLattice);

  using VelT = tmp::Key_TypePair<BulkFlag,
      NStoPF_Velocity<NSCELL, PFCELL>>;
  using VelTS = CoupledTaskSelector<std::uint8_t, NSCELL, PFCELL, VelT>;
  BlockLatManagerCoupling VelC(NSLattice, PFLattice);

  // ==============================================================
  // VTK Writers
  // ==============================================================
  vtmo::ScalarWriter PHIW("phi", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GradW("grad", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter VelW("u", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoW("rho", NSLattice.getField<RHO<T>>());
  vtmo::ScalarWriter PsiW("psi", MagLattice.getField<PSI<T>>());
  vtmo::VectorWriter HW("H", MagLattice.getField<MAGH<T, 2>>());
  vtmo::vtmWriter<T, 2> MainWriter("ferro_bubble", Geo);
  MainWriter.addWriterSet(PHIW, GradW, VelW, RhoW, PsiW, HW);

  // ==============================================================
  // Initial communication
  // ==============================================================
  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  MagLattice.NormalCommunicate();
  MagLattice.getField<PSI<T>>().Communicate();
  MagLattice.getField<MAGH<T, 2>>().Communicate();
  MainWriter.WriteBinary(0);

  // ==============================================================
  // Main loop
  // ==============================================================
  Timer MainLoopTimer;
  Timer OutputTimer;
  Printer::Print_BigBanner("Starting Ferrofluid Bubble Rise Simulation...");

  while (MainLoopTimer() < MaxStep) {

    // ----------------------------------------------------------
    // Step 1: Update PHI (sum + clamp)
    // ----------------------------------------------------------
    {
      auto& pf = PFLattice.getField<PHI<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& blk = Geo.getBlock(b);
        auto& bphi = pf.getBlockField(b);
        auto& lat = PFLattice.getBlockLat(b);
        int ol = blk.getOverlap();
        for (int j = ol; j < blk.getNy() - ol; ++j) {
          for (int i = ol; i < blk.getNx() - ol; ++i) {
            std::size_t id = j * blk.getProjection()[1] + i;
            PFCELL cell(id, lat);
            T p = T{0};
            for (int k = 0; k < 9; ++k) p += cell[k];
            if (p < T{0}) p = T{0};
            if (p > T{1}) p = T{1};
            bphi.get(id) = p;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Step 2: PF pre-coupling
    PFLattice.template ApplyCellDynamics<PFGradS>(FlagFM);
    PFLattice.template ApplyCellDynamics<PFNrmS>(FlagFM);
    PFLattice.template ApplyCellDynamics<PFLapS>(FlagFM);
    PFLattice.template ApplyCellDynamics<PFChemS>(FlagFM);

    PFLattice.getField<GRAD<T, 2>>().Communicate();
    PFLattice.getField<PHASENORMAL<T, 2>>().Communicate();
    PFLattice.getField<PHASELAPLACIAN<T>>().Communicate();
    PFLattice.getField<PHASECHEMPOTENTIAL<T>>().Communicate();

    // Step 3: PF -> Mag coupling
    MagPC.ApplyCellDynamics<MagPTS>(MainLoopTimer(), FlagFM);

    // Step 4: Mag solve
    MagLattice.template ApplyCellDynamics<MagBCS>(FlagFM);
    MagLattice.template ApplyCellDynamics<MagMRTS>(FlagFM);
    MagLattice.getField<PSI<T>>().Communicate();
    Mag_BB.Apply(MainLoopTimer());
    MagLattice.Stream();
    MagLattice.NormalCommunicate();
    // Sync Mag fields after stream (NormalCommunicate only syncs POP)
    MagLattice.getField<PSI<T>>().Communicate();
    MagLattice.getField<MAGPERMEABILITY<T>>().Communicate();
    MagLattice.getField<MAGH<T, 2>>().Communicate();

    

    MagLattice.template ApplyCellDynamics<MagGrdS>(FlagFM);
    MagLattice.template ApplyCellDynamics<MagHSqS>(FlagFM);
    MagLattice.template ApplyCellDynamics<MagFrcS>(FlagFM);

    // Step 5: NS force accumulation
    // F_total = F_s (surface tension) + F_b (buoyancy) + F_m (magnetic)
    // F_b = rho*g is handled inside PFtoNS_Forces coupling operator
    NSLattice.getField<FORCE<T, 2>>().InitValue(
        Vector<T, 2>{T{0}, T{0}});
    NSLattice.getField<FORCE<T, 2>>().Communicate();

    PropC.ApplyCellDynamics<PropTS>(MainLoopTimer(), FlagFM);
    STC.ApplyCellDynamics<STTS>(MainLoopTimer(), FlagFM);
    MagNC.ApplyCellDynamics<MagNTS>(MainLoopTimer(), FlagFM);

    NSLattice.getField<FORCE<T, 2>>().Communicate();

    // Step 6: NS solve
    NSLattice.template ApplyCellDynamics<NSMomS>(FlagFM);
    NSLattice.template ApplyCellDynamics<NSCollS>(FlagFM);
    NSLattice.getField<FORCE<T, 2>>().Communicate();
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();
    NSLattice.getField<RHO<T>>().Communicate();
    NSLattice.getField<VELOCITY<T, 2>>().Communicate();

    // Step 7: NS -> PF velocity
    VelC.ApplyCellDynamics<VelTS>(MainLoopTimer(), FlagFM);

    // Step 8: PF solve
    PFLattice.template ApplyCellDynamics<PFCollS>(FlagFM);
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();
    PFLattice.getField<PHI<T>>().Communicate();
    PFLattice.getField<PHASENORMAL<T, 2>>().Communicate();
    PFLattice.getField<GRAD<T, 2>>().Communicate();

    // Timer + Output
    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      // Sync all VTK output fields before writing
      PFLattice.getField<PHI<T>>().Communicate();
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();
      NSLattice.getField<RHO<T>>().Communicate();
      MagLattice.getField<PSI<T>>().Communicate();
      MagLattice.getField<MAGH<T, 2>>().Communicate();
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner("Ferrofluid Bubble Rise Simulation Complete!");
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}
