// bubblecoalescence2d.cpp
// 2D coaxial bubble coalescence simulation based on:
// 郝文杰 et al., "基于VOF方法的同轴气泡聚并数值模拟研究",
// Atomic Energy Science and Technology, 2025, 59(6):1262-1271
//
// Physical setup:
//   - Two coaxial (vertically aligned) equal-diameter bubbles
//   - Initial gap d0 = 2mm (~0.5*D) between bubbles
//   - Gravity drives buoyancy: lower bubble catches up to upper one
//   - Isothermal, no phase change, 2D domain
//   - Wall boundaries on all sides, bubble interior is light fluid

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
            Simulation Parameters
-----------------------------------------------*/
int Ni, Nj;
T Cell_Len;
int BlockCellLen;
int Thread_Num;

// bubble parameters
T Bubble1_Radius;
Vector<T, 2> Bubble1_Center;
T Bubble2_Radius;
Vector<T, 2> Bubble2_Center;
T Initial_Gap;

// phase field parameters
T Interface_Width;
T Mobility;

// two-phase parameters
T rho_l, rho_h;
T eta_l, eta_h;
T sigma;
T gravity;
T Eo, Re, U_g;

// simulation
int MaxStep;
int OutputStep;

void readParam() {
  iniReader param_reader("bubblecoalescence2d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Bubble1_Radius = param_reader.getValue<T>("Bubble1", "Radius");
  Bubble1_Center[0] = param_reader.getValue<T>("Bubble1", "CenterX");
  Bubble1_Center[1] = param_reader.getValue<T>("Bubble1", "CenterY");

  Bubble2_Radius = param_reader.getValue<T>("Bubble2", "Radius");
  Bubble2_Center[0] = param_reader.getValue<T>("Bubble2", "CenterX");
  Bubble2_Center[1] = param_reader.getValue<T>("Bubble2", "CenterY");
  Initial_Gap = std::abs(Bubble1_Center[1] - Bubble2_Center[1])
              - Bubble1_Radius - Bubble2_Radius;

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  Re = param_reader.getValue<T>("Two_Phase", "Re");
  U_g = param_reader.getValue<T>("Two_Phase", "U_g");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ===== 5-step parameter design =====
  T D = T(2.0) * Bubble1_Radius;

  // Step 3: gravity = U_g^2 / D
  T g_abs = U_g * U_g / D;
  gravity = -g_abs;  // negative = downward

  // Step 4: kinematic viscosity from Re
  T nu = U_g * D / Re;
  eta_h = nu * rho_h;
  eta_l = eta_h / T(10);  // gas is much less viscous

  // Step 5: surface tension from Eo
  T DeltaRho = rho_h - rho_l;
  sigma = DeltaRho * g_abs * D * D / Eo;

  // Note: Beta, Kappa, Tau_phi, Tau_ns will be computed by Converter in main()

  T cs = std::sqrt(LatSet::cs2);
  T CFL = U_g + cs;
  T Ma = U_g / cs;
  T Ca = eta_h * U_g / sigma;

  MPI_RANK(0) {
    T Mo = gravity * nu * nu * nu * rho_h * rho_h / (sigma * sigma * sigma);
    if (Mo < T(0)) Mo = -Mo;
    std::cout << "===== Bubble Coalescence Simulation =====\n"
              << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n"
              << "[Bubble]: D=" << D << "  R=" << Bubble1_Radius
              << "  Gap=" << Initial_Gap << " lu\n"
              << "  Upper: (" << Bubble1_Center[0] << "," << Bubble1_Center[1] << ")\n"
              << "  Lower: (" << Bubble2_Center[0] << "," << Bubble2_Center[1] << ")\n"
              << "[Design]: U_g=" << U_g << " (Ma=" << Ma << ")\n"
              << "  g=" << gravity << "  nu=" << nu << "  sigma=" << sigma << "\n"
              << "[Phase]: W=" << Interface_Width << "  M=" << Mobility << "\n"
              << " (Converter-derived tau/beta/kappa printed in main())\n"
              << "[Scaling]: Re=" << Re << "  Eo=" << Eo
              << "  Mo=" << Mo << "  Ca=" << Ca << "\n"
              << "[CFL]: " << CFL << (CFL < T(1.2) ? "  OK" : "  WARN")
              << "\n  rho_l=" << rho_l << "  rho_h=" << rho_h
              << "  eta_l=" << eta_l << "  eta_h=" << eta_h << "\n"
              << "[Sim]: MaxStep=" << MaxStep << "  OutputStep=" << OutputStep << "\n"
              << "==========================================\n"
              << std::endl;
  }

  if (CFL >= T(1.2)) {
    MPI_RANK(0) {
      std::cerr << "[ERROR] CFL=" << CFL << " > 1.2! Reduce U_g.\n";
    }
    exit(1);
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner("Initializing Bubble Coalescence Simulation...");
  readParam();

  // ------------------ converters ------------------
  T D = T(2.0) * Bubble1_Radius;

  // NS BaseConverter: ConvertFromTimeStep auto-computes Tau_ns = 0.5 + nu/(cs2*dx²)*dt
  // dx=1, dt=1 (standard LBM), charU=U_g for correct Lattice_charU stability check
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromTimeStep(T(1), T(1), rho_h, D, U_g, nu);

  // PhaseFieldConverter: auto-computes beta, kappa, tau_phi, omega_phi
  PhaseFieldConverter<T> PhaseConv(BaseConv, LatSet::cs2);
  PhaseConv.fromLattice(Interface_Width, Mobility, sigma);  // tau_phi computed internally

  // Register both converters for validation
  UnitConvManager<T> ConvManager(&BaseConv, &PhaseConv);
  ConvManager.Check_and_Print();

  // ----- Converter-derived parameters -----
  T Tau_ns   = BaseConv.getLattice_RT();
  T Omega_ns = BaseConv.getOMEGA();
  T Tau_phi  = PhaseConv.getLattice_RT();
  T Omega_phi = PhaseConv.getOMEGA();
  T Beta     = PhaseConv.getLatticeBeta();
  T Kappa    = PhaseConv.getLatticeKappa();

  // Validate Tau_ns from converter
  if (Tau_ns > T(1.5)) {
    MPI_RANK(0) {
      std::cerr << "[ERROR] tau_ns=" << Tau_ns << " > 1.5! Increase Re.\n";
    }
    exit(1);
  }

  MPI_RANK(0) {
    std::cout << "\n=== Converter-Derived Lattice Parameters ===\n";
    std::cout << "[NS]:  tau=" << Tau_ns   << "  omega=" << Omega_ns << "\n";
    std::cout << "[PF]:  tau=" << Tau_phi  << "  omega=" << Omega_phi << "\n";
    std::cout << "[PF]:  W="    << PhaseConv.getLatticeW()
              << "  M="         << PhaseConv.getLatticeMphi()
              << "  sigma="     << PhaseConv.getLatticeSigma() << "\n";
    std::cout << "[PF]:  beta=" << Beta   << "  kappa=" << Kappa << "\n";
    std::cout << "[Design]: Re=" << Re << "  Eo=" << Eo << "\n";
    std::cout << "=============================================\n\n";
  }

  // ------------------ geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ flag field ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Coalescence", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ NS lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(),
                         Vector<T, 2>{T{0}, T{0}},
                         T{},
                         Vector<T, 2>{T{0}, T{0}});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ------------------ PF lattice ------------------
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                            ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{},
                         Vector<T, 2>{T{0}, T{0}}, Vector<T, 2>{T{0}, T{0}},
                         Interface_Width, T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
      Geo, PFInitValues, PhaseConv,
      &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  ff::BroadcastAllParams<T>(PFLattice,
                            rho_l, rho_h, eta_l, eta_h,
                            gravity,
                            PhaseConv.getLatticeBeta(),
                            PhaseConv.getLatticeKappa());

  // ---- Initialize PHI field: two bubbles ----
  // phi = 0 inside bubble (light fluid), phi = 1 outside (heavy fluid)
  // Use min(phi1, phi2) to represent two disjoint bubbles
  T R1_phys = Bubble1_Radius * Cell_Len;
  T xc1_phys = Bubble1_Center[0] * Cell_Len;
  T yc1_phys = Bubble1_Center[1] * Cell_Len;
  T R2_phys = Bubble2_Radius * Cell_Len;
  T xc2_phys = Bubble2_Center[0] * Cell_Len;
  T yc2_phys = Bubble2_Center[1] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(bid);
    T vx = block.getVoxelSize();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + static_cast<T>(i) * vx;
        T y = block.getMin()[1] + static_cast<T>(j) * vx;

        T dist1 = std::sqrt((x - xc1_phys) * (x - xc1_phys)
                          + (y - yc1_phys) * (y - yc1_phys));
        T dist2 = std::sqrt((x - xc2_phys) * (x - xc2_phys)
                          + (y - yc2_phys) * (y - yc2_phys));

        // Each bubble: phi = 0.5 + 0.5*tanh(2*(dist-R)/W)
        T phi1 = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist1 - R1_phys) / W_phys);
        T phi2 = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist2 - R2_phys) / W_phys);

        // Union of two bubbles: min gives the most "bubble-like" value
        blockPhi.get(id) = (phi1 < phi2) ? phi1 : phi2;
      }
    }
  }

  // Initialize PF populations: f_i = w_i * phi
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    auto& blockLat = PFLattice.getBlockLat(bid);
    auto& blockPhi = phiField.getBlockField(bid);
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhi.get(id);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          cell[k] = latset::w<LatSet>(k) * phi;
        }
      }
    }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // Initialize NS populations to equilibrium (zero velocity)
  T ns_rho_init = BaseConv.getLatRhoInit();
  Vector<T, 2> u_zero{T{0}, T{0}};
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(bid);
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        NSCELL cell(id, blockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          cell[k] = latset::w<LatSet>(k) * ns_rho_init;
        }
      }
    }
  }

  // ------------------ BCs ------------------
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using FM = BlockFieldManager<FLAG, T, 2>;

  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>, LM_NS, FM>
      NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, LM_PF, FM>
      PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ tasks ------------------
  // NS task
  using NSBulkTask = tmp::Key_TypePair<
      BulkFlag, collision::MRTForce<NSCELL, FORCE<T, LatSet::d>>>;
  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask>;

  // PF tasks
  using FFNormTask = tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLapTask = tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemTask = tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using PFFFTasks = tmp::TupleWrapper<FFNormTask, FFLapTask, FFChemTask>;
  using PFFFTaskSel = tmp::TaskSelector<PFFFTasks, std::uint8_t, PFCELL>;

  // PF collision: MRTSource Allen-Cahn mode
  using PFCollTask = tmp::Key_TypePair<
      BulkFlag,
      collision::MRTSource<equilibrium::SecondOrder<PFCELL>,
                           NORMAL<T, LatSet::d>, true>>;
  using PFTaskSelector = TaskSelector<std::uint8_t, PFCELL, PFCollTask>;

  // Coupling tasks (PF -> NS)
  using STForceTask =
      tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using STForceSel =
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  using GravForceTask =
      tmp::Key_TypePair<BulkFlag, ff::FFGravityForce2D<PFCELL, NSCELL>>;
  using GravForceSel =
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravForceTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VelWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter ChemWriter("ChemPot",
      PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubblecoalescence2d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VelWriter, ChemWriter, ForceWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner("Start Bubble Coalescence Calculation...");

  while (MainLoopTimer() < MaxStep) {
    // ---- Phase field pre-processing ----
    // Step 1: Sum PF populations -> PHI
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
        const auto& block = Geo.getBlock(bid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(bid);
        auto& blockLat = PFLattice.getBlockLat(bid);
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            PFCELL cell(id, blockLat);
            T phi_new = T(0);
            for (unsigned int k = 0; k < LatSet::q; ++k)
              phi_new += cell[k];
            if (phi_new < T{0}) phi_new = T{0};
            if (phi_new > T{1}) phi_new = T{1};
            blockPhi.get(id) = phi_new;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Step 2: Compute GRAD, NORMAL, LAPLACIAN, CHEMICALPOTENTIAL
    PFLattice.template ApplyCellDynamics<PFFFTaskSel>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // ---- NS force accumulation ----
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(
        Vector<T, 2>{T{0}, T{0}});

    // Step 4: Surface tension force
    STCoupling.ApplyCellDynamics<STForceSel>(MainLoopTimer(), FlagFM);

    // Step 5: Gravity/buoyancy force
    GravCoupling.ApplyCellDynamics<GravForceSel>(MainLoopTimer(), FlagFM);

    // ---- NS collision and streaming ----
    // Step 6: NS MRTForce collision
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // Step 7: NS BCs + Stream + Communicate
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---- PF collision and streaming ----
    // Step 8: PF MRTSource collision (Allen-Cahn sharpening)
    PFLattice.template ApplyCellDynamics<PFTaskSelector>(FlagFM);

    // Step 9: PF BCs + Stream + Communicate
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner("Bubble Coalescence Simulation Complete!");
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}
