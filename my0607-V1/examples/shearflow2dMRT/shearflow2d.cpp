// shearflow2d.cpp
// 2D Droplet Deformation and Breakup in Shear Flow
// Phase-Field + NS MRT coupling
// Reference: D. Jacqmin, "Calculation of Two-Phase Navier-Stokes Flows
//            Using Phase-Field Modeling", J. Comput. Phys. 155, 96-127 (1999)
//
// Physical setup:
//   - Periodic in x, moving walls at y=0 (bottom, -U_wall) and y=H (top, +U_wall)
//   - Shear rate gamma = 2*U_wall/H
//   - Droplet placed at center
//   - Key parameters: Ca = mu*gamma*R/sigma, Re = rho*gamma*R^2/mu

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

// ===================================================================
// Global simulation parameters
// ===================================================================
int Ni, Nj;
T Cell_Len;
int BlockCellLen;
int Thread_Num;

// droplet parameters
T Droplet_Radius;
Vector<T, 2> Droplet_Center;

// shear flow
T Wall_Velocity;           // U_wall = velocity magnitude at walls
T Shear_Rate;              // gamma = 2*U_wall/H

// fluid parameters
T rho0;
T nu;
T Viscosity_Ratio;

// phase field parameters
T Interface_Width;
T sigma;
T Mobility;

// simulation control
int MaxStep;
int OutputStep;

// derived dimensionless numbers
T Ca, Re;

void readParam() {
  iniReader param_reader("shearflow2d.ini");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Droplet_Radius = param_reader.getValue<T>("Droplet", "Radius");
  Droplet_Center[0] = param_reader.getValue<T>("Droplet", "CenterX");
  Droplet_Center[1] = param_reader.getValue<T>("Droplet", "CenterY");

  Wall_Velocity = param_reader.getValue<T>("Shear_Flow", "Wall_Velocity");
  Shear_Rate = T(2.0) * Wall_Velocity / (Ni > 0 ? T(Nj) * Cell_Len : T(1));

  rho0 = param_reader.getValue<T>("Fluid", "rho");
  nu = param_reader.getValue<T>("Fluid", "nu");
  Viscosity_Ratio = param_reader.getValue<T>("Fluid", "Viscosity_Ratio");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  sigma = param_reader.getValue<T>("Phase_Field", "Surface_Tension");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  // Dimensionless numbers
  Ca = rho0 * nu * Shear_Rate * Droplet_Radius / sigma;
  Re = rho0 * Shear_Rate * Droplet_Radius * Droplet_Radius / nu;

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  MPI_RANK(0) {
    T Ma = Wall_Velocity / std::sqrt(LatSet::cs2);
    std::cout << "===== Shear Flow Droplet Deformation (MRT) =====\n"
              << "[Grid]: " << Ni << "x" << Nj
              << "  BlockCellLen=" << BlockCellLen << "\n"
              << "[Droplet]: R=" << Droplet_Radius
              << "  Center=(" << Droplet_Center[0] << "," << Droplet_Center[1] << ")\n"
              << "[Shear]: U_wall=" << Wall_Velocity
              << "  gamma=" << Shear_Rate
              << "  Ma=" << Ma << "\n"
              << "[Fluid]: rho=" << rho0 << "  nu=" << nu << "\n"
              << "[Phase]: W=" << Interface_Width << "  sigma=" << sigma
              << "  M=" << Mobility << "\n"
              << " (Converter-derived tau/beta/kappa printed in main())\n"
              << "[Dimensionless]: Ca=" << Ca << "  Re=" << Re << "\n"
              << "[Simulation]: MaxStep=" << MaxStep
              << "  OutputStep=" << OutputStep << "\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\n";
#endif
    std::cout << "=================================================\n"
              << std::endl;
  }
}

// ===================================================================
// Main
// ===================================================================
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

  Printer::Print_BigBanner("Initializing Shear Flow Droplet (MRT)...");
  readParam();

  // ===================================================================
  // Converters (lattice units: dx=1, dt=1)
  // ConvertFromTimeStep auto-computes Tau_ns = 0.5 + nu/(cs2*dx²)*dt
  // PhaseFieldConverter auto-computes beta, kappa, tau_phi, omega_phi
  // ===================================================================
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.ConvertFromTimeStep(T(1), T(1), T(1), T(Ni), Wall_Velocity, nu);

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

  MPI_RANK(0) {
    std::cout << "\n=== Converter-Derived Lattice Parameters ===\n";
    std::cout << "[NS]:  tau=" << Tau_ns   << "  omega=" << Omega_ns << "\n";
    std::cout << "[PF]:  tau=" << Tau_phi  << "  omega=" << Omega_phi << "\n";
    std::cout << "[PF]:  W="    << PhaseConv.getLatticeW()
              << "  M="         << PhaseConv.getLatticeMphi()
              << "  sigma="     << PhaseConv.getLatticeSigma() << "\n";
    std::cout << "[PF]:  beta=" << Beta   << "  kappa=" << Kappa << "\n";
    std::cout << "==========================================\n\n";
  }

  // ===================================================================
  // Geometry: domain + periodic ghosts (left/right)
  // ===================================================================
  T H_global = T(Nj) * Cell_Len;
  T L_global = T(Ni) * Cell_Len;

  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(L_global, H_global));
  AABB<T, 2> left(Vector<T, 2>(T(-Cell_Len), T(0)),
                  Vector<T, 2>(T(0), H_global));
  AABB<T, 2> right(Vector<T, 2>(L_global, T(0)),
                   Vector<T, 2>(L_global + Cell_Len, H_global));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ===================================================================
  // Flag field
  // ===================================================================
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);

  // Interior: BulkFlag
  FlagFM.forEach(
      domain, [&](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });

  // Left/right: PeriodicFlag
  FlagFM.forEach(
      left, [&](FLAG& f, std::size_t id) { f.SetField(id, PeriodicFlag); });
  FlagFM.forEach(
      right, [&](FLAG& f, std::size_t id) { f.SetField(id, PeriodicFlag); });

  // Top/bottom: BouncebackFlag (moving wall)
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_ShearFlowMRT", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ===================================================================
  // NS Lattice (D2Q9 MRT with Guo force)
  // ===================================================================
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(),
                         Vector<T, 2>{T{0}, T{0}},
                         T{},
                         Vector<T, 2>{T{0}, T{0}});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ===================================================================
  // PF Lattice (Cahn-Hilliard)
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
                         T{0}, Beta, Kappa, rho0, rho0,
                         T{}, T{});

  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
      Geo, PFInitValues, PhaseConv,
      &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  // Broadcast phase-field parameters (rho_l=rho_h=rho0 for Boussinesq)
  T gravity_zero = T(0);
  T eta_unused = nu * rho0;  // uniform viscosity for now
  ff::BroadcastAllParams<T>(PFLattice, rho0, rho0,
                            eta_unused, eta_unused,
                            gravity_zero,
                            PhaseConv.getLatticeBeta(),
                            PhaseConv.getLatticeKappa());

  // ===================================================================
  // Initialize phase field: circular droplet, phi=1 inside, 0 outside
  //   phi = 0.5 * (1 - tanh(2*(dist - R)/W))
  // ===================================================================
  T R_phys = Droplet_Radius * Cell_Len;
  T xc_phys = Droplet_Center[0] * Cell_Len;
  T yc_phys = Droplet_Center[1] * Cell_Len;
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
        T dx = x - xc_phys;
        T dy = y - yc_phys;
        T dist = std::sqrt(dx * dx + dy * dy);
        T phi = T(0.5) * (T(1) - std::tanh(T(2) * (dist - R_phys) / W_phys));
        blockPhi.get(id) = phi;
      }
    }
  }

  // Initialize PF populations: f_i = w_i * phi (zero velocity equilibrium)
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
  // Initialize NS fields: Couette linear shear profile
  //   u_x(y) = U_wall * (2*y/H - 1)
  //   Populations: equilibrium with this velocity
  // ===================================================================
  T ns_rho_init = BaseConv.getLatRhoInit();
  auto& nsVelField = NSLattice.getField<VELOCITY<T, 2>>();

  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(bid);
    auto& blockVel = nsVelField.getBlockField(bid);
    T vx = block.getVoxelSize();

    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        T y = block.getMin()[1] + static_cast<T>(j) * vx;
        // Linear shear: u_x(y) = U_wall * (2*y/H - 1)
        T ux = Wall_Velocity * (T(2) * y / H_global - T(1));
        T uy = T(0);
        Vector<T, 2> vel{ux, uy};
        blockVel.get(id) = vel;

        // Equilibrium populations
        T u2 = ux * ux + uy * uy;
        NSCELL cell(id, blockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          T uc = latset::c<LatSet>(k) * vel;
          cell[k] = latset::w<LatSet>(k) * ns_rho_init *
                    (T{1} + LatSet::InvCs2 * uc +
                     uc * uc * T{0.5} * LatSet::InvCs4 -
                     LatSet::InvCs2 * u2 * T{0.5});
        }
      }
    }
  }

  // ===================================================================
  // Boundary conditions
  // ===================================================================
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using FM = BlockFieldManager<FLAG, T, 2>;

  // NS: periodic in x
  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_Periodic("NS_Periodic", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

  // PF: periodic in x
  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_Periodic("PF_Periodic", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PF_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

#ifdef MPI_ENABLED
  NS_Periodic.SetupMPI(GeoHelper);
  PF_Periodic.SetupMPI(GeoHelper);
#endif

  // NS: moving wall bounceback on top/bottom
  BBLikeFixedBlockBdManager<bounceback::movingwall<NSCELL>, LM_NS, FM>
      NS_MovingWall("NS_MovingWall", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  // PF: standard bounceback on top/bottom (no-flux for CH equation)
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, LM_PF, FM>
      PF_BounceBack("PF_BounceBack", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ===================================================================
  // Task definitions (MRT)
  // ===================================================================

  // ---- NS collision: MRT with Guo force ----
  using NSBulkTask = tmp::Key_TypePair<
      BulkFlag, collision::MRTForce<NSCELL, FORCE<T, LatSet::d>>>;
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

  // ---- PF collision: MRTSource with Allen-Cahn sharpening (matching bubble2dMRT) ----
  using PFCollTask = tmp::Key_TypePair<
      BulkFlag,
      collision::MRTSource<equilibrium::SecondOrder<PFCELL>,
                           NORMAL<T, LatSet::d>, true>>;
  using PFPeriodicTask = tmp::Key_TypePair<
      PeriodicFlag, collision::PeriodicBoundary<PFCELL>>;
  using PFAllTasks = tmp::TupleWrapper<PFCollTask, PFPeriodicTask>;
  using PFTaskSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  // ---- Coupling tasks (PF -> NS) ----
  using STForceTask =
      tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using STForceSel =
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  // ===================================================================
  // VTK output
  // ===================================================================
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NORMALWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VELWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter ChemWriter("ChemPot",
      PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("shearflow2dMRT", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NORMALWriter,
                          VELWriter, ChemWriter, ForceWriter);

  // ===================================================================
  // Cache wall velocities (used in both initialization and main loop)
  // ===================================================================
  Vector<T, 2> u_wall_bottom{-Wall_Velocity, T(0)};
  Vector<T, 2> u_wall_top{+Wall_Velocity, T(0)};

  // ===================================================================
  // Initial communication & output
  // ===================================================================
  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  NS_Periodic.Apply();
  PF_Periodic.Apply();

  // Set wall velocities BEFORE initial moving wall bounceback
  {
    auto& velField = NSLattice.getField<VELOCITY<T, 2>>();
    for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
      const auto& block = Geo.getBlock(bid);
      auto& blockVel = velField.getBlockField(bid);
      const auto& proj = block.getProjection();
      int overlap = block.getOverlap();
      T minY = block.getMin()[1];
      T maxY = block.getMax()[1];
      if (minY < Cell_Len * T(1.5)) {
        int j_bottom = overlap;
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j_bottom * proj[1] + i;
          blockVel.get(id) = u_wall_bottom;
        }
      }
      if (maxY > H_global - Cell_Len * T(1.5)) {
        int j_top = block.getNy() - overlap - 1;
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j_top * proj[1] + i;
          blockVel.get(id) = u_wall_top;
        }
      }
    }
  }
  NS_MovingWall.Apply(0);
  PF_BounceBack.Apply(0);
  MainWriter.WriteBinary(0);

  Timer MainLoopTimer;
  Timer OutputTimer;

  Printer::Print_BigBanner("Starting Shear Flow Droplet Simulation...");

  // ===================================================================
  // Main simulation loop
  // ===================================================================
  while (MainLoopTimer() < MaxStep) {
    // ---------- Phase field pre-processing ----------
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
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              phi_new += cell[k];
            }
            // Clamp to [0, 1]
            if (phi_new < T{0}) phi_new = T{0};
            if (phi_new > T{1}) phi_new = T{1};
            blockPhi.get(id) = phi_new;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Step 2: Compute GRAD, NORMAL, LAPLACIAN, CHEMICALPOTENTIAL
    // Allen-Cahn mode: NORMAL = ∇φ/|∇φ| (interface normal)
    PFLattice.template ApplyCellDynamics<PFFFTaskSel>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // ---------- NS force accumulation ----------
    // Step 3: Clear NS FORCE
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(
        Vector<T, 2>{T{0}, T{0}});

    // Step 4: Surface tension force F_s = mu * grad-phi -> NS FORCE
    STCoupling.ApplyCellDynamics<STForceSel>(MainLoopTimer(), FlagFM);

    // ---------- NS MRT collision ----------
    // Step 5: NS MRTForce collision
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // Step 6: NS boundaries (periodic + moving wall) + stream + communicate
    // Set wall velocities BEFORE moving wall bounceback
    {
      auto& velField = NSLattice.getField<VELOCITY<T, 2>>();
      for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
        const auto& block = Geo.getBlock(bid);
        auto& blockVel = velField.getBlockField(bid);
        const auto& proj = block.getProjection();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];

        // Bottom boundary: y ~ 0
        if (minY < Cell_Len * T(1.5)) {
          int j_bottom = overlap;
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j_bottom * proj[1] + i;
            blockVel.get(id) = u_wall_bottom;
          }
        }
        // Top boundary: y ~ H_global
        if (maxY > H_global - Cell_Len * T(1.5)) {
          int j_top = block.getNy() - overlap - 1;
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j_top * proj[1] + i;
            blockVel.get(id) = u_wall_top;
          }
        }
      }
    }

    // Apply periodic + moving wall + stream + communicate
    NS_Periodic.Apply();
    NS_MovingWall.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---------- PF collision and streaming ----------
    // Step 7: PF MRTSource collision (Allen-Cahn mode, UseCHRelaxation=false)
    PFLattice.template ApplyCellDynamics<PFTaskSelector>(FlagFM);

    // Step 8: PF boundaries (periodic + bounceback) + stream + communicate
    PF_Periodic.Apply();
    PF_BounceBack.Apply(MainLoopTimer());
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

  Printer::Print_BigBanner("Shear Flow Droplet Simulation Complete!");
  MainWriter.WriteBinary(MainLoopTimer());
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();

  return 0;
}
