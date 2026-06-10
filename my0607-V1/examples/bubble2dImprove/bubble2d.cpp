// bubble2d.cpp
// 2D bubble rising simulation - Phase field + NS flow coupling
// No magnetic field, no AMR, MPI parallel

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

/*----------------------------------------------
            Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
T Cell_Len;
int BlockCellLen;
int Thread_Num;

// bubble parameters
T Bubble_Radius;
Vector<T, LatSet::d> Bubble_Center;

// phase field parameters
T Interface_Width;
T Mobility;

// two-phase parameters
T rho_l;
T rho_h;
T eta_l;
T eta_h;
T sigma;      // surface tension (lattice units)
T gravity;    // gravitational acceleration (lattice units)
T Re;         // Reynolds number
T Eo;         // Eotvos number
T U_g;        // CFL reference velocity (lattice units)
T nu;         // kinematic viscosity (lattice units)

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("bubble2d.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Bubble_Radius = param_reader.getValue<T>("Bubble", "Radius");
  Bubble_Center[0] = param_reader.getValue<T>("Bubble", "CenterX");
  Bubble_Center[1] = param_reader.getValue<T>("Bubble", "CenterY");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");
  Re = param_reader.getValue<T>("Two_Phase", "Re");
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  T viscosity_ratio = param_reader.getValue<T>("Two_Phase", "viscosity_ratio");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  T D = T(2.0) * Bubble_Radius;
  T DeltaRho = rho_h - rho_l;
  T cs = std::sqrt(LatSet::cs2);
  // --- Dimensionless design: Re, Eo, U_g → g, nu, sigma ---
  U_g     = param_reader.getValue<T>("Two_Phase", "U_g");
  T g_abs = U_g * U_g / D;
  nu      = U_g * D / Re;
  gravity = -g_abs;
  sigma   = DeltaRho * g_abs * D * D / Eo;
  eta_h   = nu * rho_h;
  eta_l   = eta_h / viscosity_ratio;

  // Derived dimensionless numbers (overwrite globals for converter output)
  T g_abs_print = std::abs(gravity);
  Re = U_g * D / nu;
  Eo = DeltaRho * g_abs_print * D * D / sigma;

  // CFL check
  T CFL = U_g + cs;
  bool cfl_ok = (CFL < T(1.2));

  MPI_RANK(0) {
    std::cout << "----------Bubble Rising Simulation----------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius << " D=" << D
              << "  Center=(" << Bubble_Center[0] << "," << Bubble_Center[1] << ")\n";
    std::cout << "[Params]: nu=" << nu << "  g=" << gravity
              << "  sigma=" << sigma << "\n";
    std::cout << "[Params]: rho_l=" << rho_l << "  rho_h=" << rho_h
              << "  eta_h=" << eta_h << "  eta_l=" << eta_l << "\n";
    std::cout << "[Ref]:   U_g=" << U_g << " (Ma=" << U_g/cs << ")"
              << "  Re=" << Re << "  Eo=" << Eo << "\n";
    std::cout << "[CFL]: (U_g+cs)*dt/dx = " << CFL
              << (cfl_ok ? "  OK" : "  WARNING: >1.2!") << "\n";
    std::cout << "--------------------------------------------\n";
  }

  if (!cfl_ok) {
    MPI_RANK(0) {
      std::cerr << "[ERROR] CFL=" << CFL << " > 1.2! Reduce U_g in ini file.\n";
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

  Printer::Print_BigBanner(std::string("Initializing Bubble Rising..."));

  readParam();

  // ------------------ define converters ------------------
  T D = T(2.0) * Bubble_Radius;

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

  IF_MPI_RANK(0) {
    std::cout << "\n=== Converter-Derived Lattice Parameters ===\n";
    std::cout << "[NS]:  tau=" << Tau_ns   << "  omega=" << Omega_ns << "\n";
    std::cout << "[PF]:  tau=" << Tau_phi  << "  omega=" << Omega_phi << "\n";
    std::cout << "[PF]:  W="    << PhaseConv.getLatticeW()
              << "  M="         << PhaseConv.getLatticeMphi()
              << "  sigma="     << PhaseConv.getLatticeSigma() << "\n";
    std::cout << "[PF]:  beta=" << Beta   << "  kappa=" << Kappa << "\n";
    std::cout << "[Design]: Re=" << Re << "  Eo=" << Eo
              << "  U_g=" << U_g << "\n";
    std::cout << "=============================================\n\n";
  }

  // ------------------ define geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1,mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // Write flag geometry for verification
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Bubble", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define NS lattice ------------------
  // Note: RHO<T> is cosmetic only (forcerhoU overwrites it with Σf≈1.0 each step)
  // Actual density variation enters via the FORCE field
  // NS fields: RHO, PRESSURE, VELOCITY, POP, FORCE, OMEGA (per-cell)
  //   RHO_L, RHO_H, ETA_L, ETA_H, GRAVITY (Data, fluid properties)
  using NSFIELDS = TypePack<RHO<T>, PRESSURE<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, OMEGA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>,
                            ff::ETA_L<T>, ff::ETA_H<T>, ff::GRAVITY<T>>;
  using NSFIELDREFS = TypePack<PHI<T>>;  // from PFLattice via addField
  using NSFIELDPACK = TypePack<NSFIELDS, NSFIELDREFS>;
  using NSALLFIELDS = TypePack<RHO<T>, PRESSURE<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                               FORCE<T, LatSet::d>, OMEGA<T>,
                               ff::RHO_L<T>, ff::RHO_H<T>,
                               ff::ETA_L<T>, ff::ETA_H<T>, ff::GRAVITY<T>,
                               PHI<T>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), T{}, Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}}, T{},
                         rho_l, rho_h, eta_l, eta_h, gravity);
  using NSCELL = Cell<T, LatSet, NSALLFIELDS>;

  BlockFieldManager<PHI<T>, T, 2>* nsPhiPtr = nullptr;
  BlockLatticeManager<T, LatSet, NSFIELDPACK> NSLattice(
    Geo, NSInitValues, BaseConv, nsPhiPtr);

  // ------------------ define PF lattice ------------------
  // PFLattice OWNS: PHI, POP, GRAD, CHEMICALPOTENTIAL
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, ff::GRAD<T>, ff::CHEMICALPOTENTIAL<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, 2>, ff::NORMAL<T>, ff::LAPLACIAN<T>, ff::ACSOURCE<T>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  using PFALLFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>,
                               ff::GRAD<T>, ff::CHEMICALPOTENTIAL<T>,
                               VELOCITY<T, 2>, ff::NORMAL<T>, ff::LAPLACIAN<T>, ff::ACSOURCE<T>>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}}, T{});
  using PFCELL = Cell<T, LatSet, PFALLFIELDS>;

  // nullptr placeholders for FF2DMgr-owned fields (filled via addField)
  BlockFieldManager<ff::NORMAL<T>, T, 2>* tempNormPtr = nullptr;
  BlockFieldManager<ff::LAPLACIAN<T>, T, 2>* tempLapPtr = nullptr;
  BlockFieldManager<ff::ACSOURCE<T>, T, 2>* tempSrcPtr = nullptr;

  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PhaseConv,
    &NSLattice.getField<VELOCITY<T, 2>>(),
    tempNormPtr, tempLapPtr, tempSrcPtr);

  // ---- Initialize PHI field ----
  // φ = 0 for light fluid (bubble interior), φ = 1 for heavy fluid (outside)
  // φ = 0.5 + 0.5*tanh(2*(r-R)/W)
  T R_phys = Bubble_Radius * Cell_Len;
  T xc_phys = Bubble_Center[0] * Cell_Len;
  T yc_phys = Bubble_Center[1] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
        T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
        T dx = x - xc_phys;
        T dy = y - yc_phys;
        T dist = std::sqrt(dx * dx + dy * dy);
        T phi = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
        blockPhi.get(id) = phi;
      }
    }
  }

  // ---- Initialize PF distribution: g_i = w_i * phi (zero-velocity equilibrium) ----
  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = PFLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    const auto& proj = block.getProjection();
    int overlap = block.getOverlap();
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhiField.get(id);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          cell[k] = latset::w<LatSet>(k) * phi;
        }
      }
    }
  }

  // === Initialize NS macro fields ===
  NSLattice.getField<VELOCITY<T, LatSet::d>>().InitValue(Vector<T, LatSet::d>{T{0}, T{0}});
  NSLattice.getField<RHO<T>>().InitValue(BaseConv.getLatRhoInit());
  NSLattice.getField<PRESSURE<T>>().InitValue(BaseConv.getLatRhoInit() * LatSet::cs2);

  // ------------------ define FF2D manager ------------------
  // FF2DManager OWNS: NORMAL, LAPLACIAN, ACSOURCE, FORCE_SF, FORCE_Grav, FORCE_Buoy, FORCE_Visc
  // FF2DManager REFERENCES: PHI, GRAD, CHEMICALPOTENTIAL, FORCE(from NS), VELOCITY(from NS)
  ValuePack FFInitValues(
      Vector<T, 2>{T{0}, T{0}}, T{}, Vector<T, 2>{T{0}, T{0}},
      Vector<T, 2>{T{0}, T{0}}, Vector<T, 2>{T{0}, T{0}},
      Vector<T, 2>{T{0}, T{0}}, Vector<T, 2>{T{0}, T{0}});
  ff::FF2DManager<T, LatSet> FF2DMgr(
      Geo, PhaseConv, FFInitValues,
      &PFLattice.getField<PHI<T>>(),
      &PFLattice.getField<ff::GRAD<T>>(),
      &PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>(),
      &NSLattice.getField<ff::RHO_L<T>>(),
      &NSLattice.getField<ff::RHO_H<T>>(),
      &NSLattice.getField<ff::ETA_L<T>>(),
      &NSLattice.getField<ff::ETA_H<T>>(),
      &NSLattice.getField<ff::GRAVITY<T>>(),
      &NSLattice.getField<RHO<T>>(),
      &NSLattice.getField<PRESSURE<T>>(),
      &NSLattice.getField<FORCE<T, LatSet::d>>(),
      &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  FF2DMgr.Init(GeoHelper);

  // PFLattice references FF2DMgr-owned fields (like SOLattice references CA::EXCESSC)
  PFLattice.template addField<ff::NORMAL<T>>(FF2DMgr.template getField<ff::NORMAL<T>>());
  PFLattice.template addField<ff::LAPLACIAN<T>>(FF2DMgr.template getField<ff::LAPLACIAN<T>>());

  // NSLattice references PHI from PFLattice (for density/viscosity interpolation)
  NSLattice.template addField<PHI<T>>(PFLattice.template getField<PHI<T>>());
  PFLattice.template addField<ff::ACSOURCE<T>>(FF2DMgr.template getField<ff::ACSOURCE<T>>());

  // ------------------ define BCs ------------------
  // NS: bounceback on all walls
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDPACK>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  // PF: bounceback on all walls
  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ define tasks ------------------
  // PF collision: g* = ω·geq + (1-ω)·g + (1-ω/2)·Si  (pure collision)
  using PFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::PhaseFieldBGKSource<
      equilibrium::PhaseFieldEquilibrium<PFCELL>,
      ff::PhaseFieldSource<PFCELL>>>;
  using PFCollisionTaskSelector = TaskSelector<std::uint8_t, PFCELL, PFCollisionTask>;

  // NS collision: f* = ω·feq + (1-ω)·f + (1-ω/2)·Gi  (pure collision)
  using NSCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::NSBGKForce<
      equilibrium::SecondOrder<NSCELL>,
      force::Force<NSCELL>>>;
  using NSCollisionTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSCollisionTask>;


  // NS density interpolation: ρ(φ) = ρ_l + φ·(ρ_h-ρ_l)
  using RhoInterpTask = tmp::Key_TypePair<
    BulkFlag,
    moment::rhoInterp<NSCELL, PHI<T>, ff::RHO_L<T>, ff::RHO_H<T>>>;
  using RhoInterpTaskSelector = TaskSelector<std::uint8_t, NSCELL, RhoInterpTask>;

  // NS omega interpolation: η(φ)→ν→ω = 1/(ν/cs²+0.5)
  using OmegaInterpTask = tmp::Key_TypePair<
    BulkFlag,
    moment::omegaInterp<NSCELL, PHI<T>,
                        ff::RHO_L<T>, ff::RHO_H<T>,
                        ff::ETA_L<T>, ff::ETA_H<T>>>;
  using OmegaInterpTaskSelector = TaskSelector<std::uint8_t, NSCELL, OmegaInterpTask>;

  // NS pressure + velocity: u = (Σc·f + 0.5F)/ρ,  p = ρ·cs²
  using PressUTask = tmp::Key_TypePair<
    BulkFlag,
    moment::pressU<NSCELL, force::Force<NSCELL>>>;
  using PressUTaskSelector = TaskSelector<std::uint8_t, NSCELL, PressUTask>;

  using PhiMomentTask = tmp::Key_TypePair<
    BulkFlag,
    moment::PhaseFieldMomenta<PFCELL>>;
  using PhiMomentTaskSelector = TaskSelector<std::uint8_t, PFCELL, PhiMomentTask>;

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<ff::GRAD<T>>());
  vtmo::VectorWriter NormalWriter("NORMAL", FF2DMgr.getField<ff::NORMAL<T>>());
  vtmo::ScalarWriter LaplacianWriter("LAPLACIAN", FF2DMgr.getField<ff::LAPLACIAN<T>>());
  vtmo::ScalarWriter ChemPotWriter("CHEMPOT", PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubble2d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter, LaplacianWriter,
                          ChemPotWriter, VecWriter, RhoWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  MainWriter.WriteBinary(MainLoopTimer());


  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep) {
    // ---- Phase 1: FF computations (gradient, chemical potential, forces) ----
    FF2DMgr.computeGradient();
    FF2DMgr.computeLaplacian();
    FF2DMgr.computeChemicalPotential();
    FF2DMgr.computeACSource();

    // Communicate FF fields needed for force computation
    PFLattice.getField<ff::GRAD<T>>().Communicate();
    PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>().Communicate();
    FF2DMgr.Communicate();  // NORMAL, LAPLACIAN

    // Clear NS FORCE, compute & sum all force components
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});
    FF2DMgr.computeForceSF();
    // FF2DMgr.computeForceP();
    // FF2DMgr.computeForceBuoy();
    // FF2DMgr.computeForceVisc();
    FF2DMgr.computeForceTotal();

    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // ---- Phase 2: Collision ----
    PFLattice.template ApplyCellDynamics<PFCollisionTaskSelector>(FlagFM);
    NSLattice.template ApplyCellDynamics<NSCollisionTaskSelector>(FlagFM);

    // ---- Phase 3: Streaming + BC ----
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();

    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---- Phase 4: Macro updates after streaming ----
    // PF: φ = Σg_i
    PFLattice.ApplyCellDynamics<PhiMomentTaskSelector>(FlagFM);
    PFLattice.getField<PHI<T>>().Communicate();

    // NS: ρ(φ) → ω(φ) → u = (Σc·f+0.5F)/ρ,  p = ρ·cs²
    NSLattice.template ApplyCellDynamics<RhoInterpTaskSelector>(FlagFM);
    NSLattice.template ApplyCellDynamics<OmegaInterpTaskSelector>(FlagFM);
    NSLattice.template ApplyCellDynamics<PressUTaskSelector>(FlagFM);

    // Communicate updated macro fields
    NSLattice.getField<RHO<T>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<VELOCITY<T, 2>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<PHI<T>>().Communicate();
      PFLattice.getField<ff::GRAD<T>>().Communicate();
      PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>().Communicate();
      FF2DMgr.Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();
      NSLattice.getField<RHO<T>>().Communicate();
      NSLattice.getField<PRESSURE<T>>().Communicate();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}
