// rti2d.cpp
// 2D Rayleigh-Taylor Instability simulation
// Based on: Celani et al., "Phase-field model for the Rayleigh-Taylor
//           instability of immiscible fluids" (2018)
//
// Dual-lattice coupling pattern follows bubble2d.cpp:
// - PF lattice: phase field via Cahn-Hilliard (BGKSource)
// - NS lattice: Navier-Stokes with forces (BGKForce)
// - Velocity: PF references NS VELOCITY field (auto-sync)
// - Coupling: BlockLatManagerCoupling + CoupledTaskSelector

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

// fluid parameters
T rho_l;            // ρ₁ (lighter fluid, below)
T rho_h;            // ρ₂ (heavier fluid, above)
T gravity;          // negative = downward
T nu;               // kinematic viscosity
T Tau_ns;

// phase field parameters
T Interface_Width;  // W
T sigma;            // surface tension
T Mobility;
T Tau_phi;
T Kappa;
T Beta;

// perturbation
T Pert_Amp;         // h₀
int WaveNumber;     // mode number k

// simulation
int MaxStep;
int OutputStep;

void readParam() {
  iniReader param_reader("rti2d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  rho_l = param_reader.getValue<T>("Fluid", "Rho_Light");
  rho_h = param_reader.getValue<T>("Fluid", "Rho_Heavy");
  gravity = param_reader.getValue<T>("Fluid", "Gravity");
  nu = param_reader.getValue<T>("Fluid", "Kinematic_Viscosity");
  Tau_ns = T(0.5) + nu / LatSet::cs2;

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  sigma = param_reader.getValue<T>("Phase_Field", "Surface_Tension");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");
  Kappa = T(1.5) * Interface_Width * sigma;
  Beta = T(12.0) * sigma / Interface_Width;
  Tau_phi = T(0.5) + Mobility / LatSet::cs2;

  Pert_Amp = param_reader.getValue<T>("Perturbation", "Amplitude");
  WaveNumber = param_reader.getValue<int>("Perturbation", "WaveNumber");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  MPI_RANK(0) {
    T Atwood = (rho_h - rho_l) / (rho_h + rho_l);
    T k_phys = T(2.0 * M_PI * WaveNumber) / (Ni * Cell_Len);
    T k_c = std::sqrt(T(2.0) * Atwood * std::abs(gravity) * (rho_l + rho_h) * T(0.5) / sigma);
    T alpha_sq = Atwood * std::abs(gravity) * k_phys
               - sigma / (rho_l + rho_h) * k_phys * k_phys * k_phys;
    std::cout << "======= RTI 2D Simulation (Celani et al., 2018) =======\n"
              << "[Grid]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n"
              << "[Fluid]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " A=" << Atwood << " nu=" << nu << "\n"
              << "[Gravity]: " << gravity << "\n"
              << "[Phase]: W=" << Interface_Width << " sigma=" << sigma
              << " kappa=" << Kappa << " beta=" << Beta
              << " Mobility=" << Mobility << " tau_phi=" << Tau_phi << "\n"
              << "[NS]: tau_ns=" << Tau_ns << " omega=" << (T(1)/Tau_ns) << "\n"
              << "[Perturbation]: h0=" << Pert_Amp << " k=" << WaveNumber
              << " h0/lambda=" << (Pert_Amp * WaveNumber / (Ni * Cell_Len))
              << " h0/W=" << (Pert_Amp / Interface_Width) << "\n"
              << "[Dispersion]: k_phys=" << k_phys << " k_c=" << k_c
              << " alpha_sq=" << alpha_sq
              << (alpha_sq > 0 ? " UNSTABLE" : " stable") << "\n"
              << "[Sim]: MaxStep=" << MaxStep << " OutputStep=" << OutputStep << "\n"
              << "==========================================================\n"
              << std::endl;
  }
}

// ===================================================================
// Custom Boussinesq buoyancy functor for RTI
// Uses (ρ(φ) - ρ₀) * g  instead of (ρ - ρ_h) * g
// Both fluids experience opposite forces → proper RTI dynamics
// ===================================================================
template <typename PFCELL, typename NSCELL>
struct RTIBuoyancy2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    T phi = pf_cell.template get<typename PFCELL::GenericRho>();
    T _rho_l = pf_cell.template get<ff::RHO_L<T>>();
    T _rho_h = pf_cell.template get<ff::RHO_H<T>>();
    T _gravity = pf_cell.template get<ff::GRAVITY<T>>();

    T rho = (T(1) - phi) * _rho_h + phi * _rho_l;
    T rho_0 = (_rho_l + _rho_h) * T(0.5);
    // F_b = (ρ(φ) - ρ₀) * g  (Boussinesq, Eq.(2.8) with φ ∈ [0,1])
    ns_cell.template get<FORCE<T, LatSet::d>>()[1] += (rho - rho_0) * _gravity;
  }
};

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t PeriodicFlag = std::uint8_t(8);

  mpi().init(&argc, &argv);
  if (mpi().getRank() == 0) {
    std::cout << "Total MPI processes: " << mpi().getSize() << std::endl;
  }
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner("Initializing RTI 2D...");
  readParam();

  // ===================================================================
  // Converters (lattice units: dx=1, dt=1)
  // ===================================================================
  BaseConverter<T> NSConv(LatSet::cs2);
  NSConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);

  BaseConverter<T> PFConv(LatSet::cs2);
  PFConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  UnitConvManager<T> ConvManager(&NSConv);
  ConvManager.Check_and_Print();

  // ===================================================================
  // Geometry (periodic in X, no-slip walls in Y)
  // ===================================================================
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  // X-periodic ghost columns (covers interior y range only)
  AABB<T, 2> left(Vector<T, 2>(T(-Cell_Len), T(0)),
                  Vector<T, 2>(T(0), T(Nj * Cell_Len)));
  AABB<T, 2> right(Vector<T, 2>(T(Ni * Cell_Len), T(0)),
                   Vector<T, 2>(T((Ni + 1) * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ===================================================================
  // Flag field
  // X: periodic (left/right boxes), Y: no-slip walls (bounce-back)
  // ===================================================================
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(
      domain, [&](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });
  // X-periodic ghost columns (y range: interior only, not top/bottom ghost rows)
  FlagFM.forEach(
      left, [&](FLAG& f, std::size_t id) { f.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
      right, [&](FLAG& f, std::size_t id) { f.SetField(id, PeriodicFlag); });
  // Mark Y-boundary cells as BouncebackFlag
  // (top/bottom edges of domain, where neighbor is VoidFlag)
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_RTI", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ===================================================================
  // NS Lattice: Navier-Stokes with force
  // ===================================================================
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>>;
  ValuePack NSInitValues(NSConv.getLatRhoInit(),
                         Vector<T, 2>{T{0}, T{0}},
                         T{},
                         Vector<T, 2>{T{0}, T{0}});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, NSConv);

  // Initialize NS POPs to equilibrium (rho=ρ₀, u=0)
  T ns_rho_init = NSConv.getLatRhoInit();
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(bid);
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        NSCELL cell(id, blockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          T feq = latset::w<LatSet>(k) * ns_rho_init;
          cell[k] = feq;
        }
      }
    }
  }

  // ===================================================================
  // PF Lattice: Phase field (Cahn-Hilliard)
  // VELOCITY references NS VELOCITY → auto-sync (no manual copying)
  // ===================================================================
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                            ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>,
                            ff::ETA_L<T>, ff::ETA_H<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{},
                         Vector<T, 2>{T{0}, T{0}}, Vector<T, 2>{T{0}, T{0}},
                         Interface_Width, T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h,
                         T{}, T{});

  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
      Geo, PFInitValues, PFConv,
      &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  T eta_unused = T(0);
  ff::BroadcastAllParams<T>(PFLattice, rho_l, rho_h,
                            eta_unused, eta_unused,
                            gravity, Beta, Kappa);

  // ===================================================================
  // Initialize phase field: φ = 0.5*(1 - tanh(2*(y - y_interface)/W))
  // φ=1: light fluid (below), φ=0: heavy fluid (above)
  // Interface perturbation: h(x) = h₀ sin(2π k x / Lx)
  // ===================================================================
  T MidY = T(Nj * Cell_Len * 0.5);
  T Lx = T(Ni * Cell_Len);
  T k0 = T(2.0 * M_PI * WaveNumber / Lx);

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(bid);
    auto& blockLat = PFLattice.getBlockLat(bid);
    T vx = block.getVoxelSize();

    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + static_cast<T>(i) * vx;
        T y = block.getMin()[1] + static_cast<T>(j) * vx;
        T h = Pert_Amp * std::sin(k0 * x);
        T dist = y - MidY - h;
        T phi = T(0.5) * (T(1) - std::tanh(T(2) * dist / Interface_Width));
        blockPhi.get(id) = phi;
      }
    }
  }

  // Init PF POPs to equilibrium (u=0)
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    auto& blockLat = PFLattice.getBlockLat(bid);
    auto& blockPhi = phiField.getBlockField(bid);
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();

    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhi.get(id);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          cell[k] = latset::w<LatSet>(k) * phi;
        }
      }
    }
  }

  // ===================================================================
  // Boundary conditions
  // X: periodic (left ↔ right), Y: bounce-back (no-slip walls)
  // ===================================================================
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using FM = BlockFieldManager<FLAG, T, 2>;

  // X-periodic for NS
  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_Periodic("NS_Periodic", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

  // X-periodic for PF
  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_Periodic("PF_Periodic", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PF_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

#ifdef MPI_ENABLED
  NS_Periodic.SetupMPI(GeoHelper);
  PF_Periodic.SetupMPI(GeoHelper);
#endif

  // Y-direction bounce-back (no-slip walls) for NS and PF
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>, LM_NS, FM>
      NS_BB("NS_BounceBack", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, LM_PF, FM>
      PF_BB("PF_BounceBack", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ===================================================================
  // Task definitions
  // ===================================================================

  // ---- NS collision: BGKForce with vector force ----
  using NSBulkTask = tmp::Key_TypePair<
      BulkFlag,
      collision::BGKForce<
          moment::forcerhoU<NSCELL, force::Force<NSCELL>, true>,
          equilibrium::SecondOrder<NSCELL>,
          force::Force<NSCELL>>>;
  using NSPeriodicTask = tmp::Key_TypePair<
      PeriodicFlag, collision::PeriodicBoundary<NSCELL>>;
  using NSAllTasks = tmp::TupleWrapper<NSBulkTask, NSPeriodicTask>;
  using NSTaskSelector = tmp::TaskSelector<NSAllTasks, std::uint8_t, NSCELL>;

  // ---- PF pre-processing: FF2D, Laplacian, Chemical potential ----
  using FFNormTask = tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLapTask = tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemTask = tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using PFFFTasks = tmp::TupleWrapper<FFNormTask, FFLapTask, FFChemTask>;
  using PFFFTaskSel = tmp::TaskSelector<PFFFTasks, std::uint8_t, PFCELL>;

  // ---- PF chemical potential gradient → NORMAL (Cahn-Hilliard source) ----
  using FFChemGradTask =
      tmp::Key_TypePair<BulkFlag, ff::FFChemPotentialGradient2D<PFCELL>>;
  using FFChemGradSel = TaskSelector<std::uint8_t, PFCELL, FFChemGradTask>;

  // ---- PF collision: BGKSource with NORMAL (= ∇λ) as CH source ----
  using PFCollTask = tmp::Key_TypePair<
      BulkFlag,
      collision::BGKSource<equilibrium::SecondOrder<PFCELL>,
                           NORMAL<T, LatSet::d>, true>>;
  using PFPeriodicTask = tmp::Key_TypePair<
      PeriodicFlag, collision::PeriodicBoundary<PFCELL>>;
  using PFAllTasks = tmp::TupleWrapper<PFCollTask, PFPeriodicTask>;
  using PFTaskSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  // ---- Coupling tasks (PF → NS), using CoupledTaskSelector ----
  using STForceTask =
      tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using STForceSel = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  using RTIGravTask =
      tmp::Key_TypePair<BulkFlag, RTIBuoyancy2D<PFCELL, NSCELL>>;
  using RTIGravSel = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RTIGravTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  // ===================================================================
  // VTK output
  // ===================================================================
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter VELWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter ChemWriter("ChemPot", PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>());
  vtmo::vtmWriter<T, 2> MainWriter("rti2d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, VELWriter, ChemWriter);

  // ===================================================================
  // Initial communication & output
  // ===================================================================
  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  NS_Periodic.Apply();
  PF_Periodic.Apply();
  NS_BB.Apply(0);
  PF_BB.Apply(0);
  MainWriter.WriteBinary(0);

  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner("Starting RTI 2D Calculation...");

  // ===================================================================
  // Main simulation loop
  // (follows bubble2d dual-lattice coupling pattern)
  // ===================================================================
  while (MainLoopTimer() < MaxStep) {
    // ---------- Phase field pre-processing ----------
    // Step 1: Sum PF populations → PHI
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
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              phi_new += cell[k];
            }
            // Clamp to suppress numerical noise
            if (phi_new < T{0}) phi_new = T{0};
            if (phi_new > T{1}) phi_new = T{1};
            blockPhi.get(id) = phi_new;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Step 2: Compute ∇φ, ∇²φ, chemical potential μ
    PFLattice.template ApplyCellDynamics<PFFFTaskSel>(FlagFM);
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // Step 3: Compute ∇μ → NORMAL (Cahn-Hilliard source direction)
    PFLattice.template ApplyCellDynamics<FFChemGradSel>(FlagFM);

    // ---------- NS force accumulation ----------
    // Step 4: Clear NS FORCE field
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(
        Vector<T, 2>{T{0}, T{0}});

    // Step 5: Surface tension force  F_s = μ ∇φ → NS FORCE
    STCoupling.ApplyCellDynamics<STForceSel>(MainLoopTimer(), FlagFM);

    // Step 6: Boussinesq buoyancy  F_b = (ρ-ρ₀) g → NS FORCE
    GravCoupling.ApplyCellDynamics<RTIGravSel>(MainLoopTimer(), FlagFM);

    // ---------- NS collision + streaming ----------
    // Step 7: NS BGKForce collision (reads FORCE from same cell, no comm needed)
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);

    // Step 8: NS boundaries, stream, communicate
    NS_Periodic.Apply();
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---------- PF collision + streaming ----------
    // Step 9: PF BGKSource collision
    //   (reads NORMAL=∇λ from same cell, VELOCITY auto-synced via PFFIELDREFS)
    PFLattice.template ApplyCellDynamics<PFTaskSelector>(FlagFM);

    // Step 10: PF boundaries, stream, communicate
    PF_Periodic.Apply();
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
      PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>().Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner("RTI 2D Calculation Complete!");
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", NSConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}
