// bubble3dMRT.cpp
// 3D bubble rising simulation - Phase field + NS flow coupling
// No magnetic field, no AMR, MPI parallel

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D3Q19<T>;

/*----------------------------------------------
            Simulation Parameters
-----------------------------------------------*/
int Ni;
int Nj;
int Nz;
T Cell_Len;
int BlockCellLen;
int Thread_Num;

// bubble parameters
T Bubble_Radius;
Vector<T, LatSet::d> Bubble_Center;

// phase field parameters
T Interface_Width;
T Mobility;
T Tau_phi;
T Omega_phi;
T Kappa;
T Beta;

// two-phase parameters
T rho_l;
T rho_h;
T eta_l;
T eta_h;
T sigma;
T gravity;
T Eo;
T Re;
T U_g;       // characteristic lattice velocity
T Tau_ns;

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("bubble3dMRT.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Nz = param_reader.getValue<int>("Mesh", "Nz");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Bubble_Radius = param_reader.getValue<T>("Bubble", "Radius");
  Bubble_Center[0] = param_reader.getValue<T>("Bubble", "CenterX");
  Bubble_Center[1] = param_reader.getValue<T>("Bubble", "CenterY");
  Bubble_Center[2] = param_reader.getValue<T>("Bubble", "CenterZ");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");
  eta_l = param_reader.getValue<T>("Two_Phase", "eta_l");
  eta_h = param_reader.getValue<T>("Two_Phase", "eta_h");
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  Re = param_reader.getValue<T>("Two_Phase", "Re");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ===== Scheme A: derive g, sigma, U_g from Re/Eo/eta (Safi et al. TC2) =====
  // Re = rho_h * sqrt(g) * (2*r0)^(3/2) / eta_h
  // Eo = 4 * rho_h * g * r0^2 / sigma
  // If g/sigma/U_g are explicitly provided, they override the derived values.

  T r0 = Bubble_Radius * Cell_Len;

  if (param_reader.hasKey("Two_Phase", "g")) {
    gravity = param_reader.getValue<T>("Two_Phase", "g");
  } else {
    gravity = std::pow(Re * eta_h / (rho_h * std::pow(T(2) * r0, T(1.5))), T(2));
  }
  if (param_reader.hasKey("Two_Phase", "sigma")) {
    sigma = param_reader.getValue<T>("Two_Phase", "sigma");
  } else {
    sigma = T(4) * rho_h * gravity * r0 * r0 / Eo;
  }
  if (param_reader.hasKey("Two_Phase", "U_g")) {
    U_g = param_reader.getValue<T>("Two_Phase", "U_g");
  } else {
    U_g = std::sqrt(gravity * T(2) * r0);
  }

  // Derived: Eq.(14) β=12σ/W, κ=3Wσ/2
  Beta = T(12.0) * sigma / Interface_Width;          // ~0.0612 (LB units)
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5); // ~0.1224 (LB units)
  // tau_phi = 3*Mobility + 0.5 = 0.53 (Fortran tau_ddg)
  Tau_phi = T(3.0) * Mobility + T(0.5);
  Omega_phi = T(1.0) / Tau_phi;

  // Base NS tau (heavy fluid): tau = 0.5 + 3*nu_h = 0.5 + 3*vis_high/rho_h
  // Fortran tau_ddf_high = 3*0.6/35 + 0.5 = 0.55143
  Tau_ns = T(0.5) + eta_h / rho_h / LatSet::cs2;

  // LBM compressibility check: Ma = U_g / cs
  T cs = std::sqrt(LatSet::cs2);
  T Ma = U_g / cs;

  MPI_RANK(0) {
    std::cout << "-------Bubble Rising Simulation (Scheme A / TC2)-------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "x" << Nz
              << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius
              << "  Center=(" << Bubble_Center[0] << ","
              << Bubble_Center[1] << "," << Bubble_Center[2] << ")\n";
    std::cout << "[Input]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " eta_l=" << eta_l << " eta_h=" << eta_h
              << " Re=" << Re << " Eo=" << Eo << "\n";
    std::cout << "[Derived]: g=" << gravity << " sigma=" << sigma
              << " U_g=" << U_g << " Ma=" << Ma << "\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " beta=" << Beta << " kappa=" << Kappa
              << " tau_phi=" << Tau_phi << "\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep << "  OutputStep=" << OutputStep << "\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\n";
#endif
    std::cout << "-------------------------------------------------------\n";
  }

  if (Ma > T(0.2)) {
    MPI_RANK(0) {
      std::cerr << "[Warning] Ma=" << Ma << " > 0.2, may affect incompressibility\n";
    }
  }

}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Bubble Rising (MRT)..."));

  readParam();

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);

  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ define geometry ------------------
  AABB<T, LatSet::d> domain(Vector<T, LatSet::d>(T(0), T(0), T(0)),
                            Vector<T, LatSet::d>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nz * Cell_Len)));

  BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nz, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize(), 1);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry3D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  // No-slip on all six boundaries (Safi et al. TC2)
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // Write flag geometry for verification
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, LatSet::d> GeoWriter("GeoFlag_Bubble", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define NS lattice ------------------
  // Note: RHO<T> is cosmetic only (forcerhoU overwrites it with Σf≈1.0 each step)
  // Actual density variation enters via the FORCE field
  using NSFIELDS = TypePack<DENSITY<T>, VELOCITY<T, LatSet::d>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, OMEGA<T>, PRESSURE<T>>;
  T omega_ns_init = T{1} / Tau_ns;
  ValuePack NSInitValues(T{1}, Vector<T, LatSet::d>{T{0}},
                         T{}, Vector<T, LatSet::d>{T{0}}, omega_ns_init, T{});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ------------------ define PF lattice ------------------
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                            ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>,
                            ff::DELTARHO<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, LatSet::d>{T{0}},
                         Vector<T, LatSet::d>{T{0}}, Interface_Width,
                         T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h,
                         rho_h - rho_l);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PFBaseConv, &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  // Type B params: broadcast from rank 0 to all ranks (MPI-safe)
  ff::BroadcastAllParams<T>(PFLattice,
                            rho_l, rho_h, eta_l, eta_h,
                            gravity, Beta, Kappa);

  // Init DELTARHO (block-level constant for variable density)
  PFLattice.template getField<ff::DELTARHO<T>>().InitValue(rho_h - rho_l);

  // ---- Initialize PHI and POP fields ----
  // φ = 0 for light fluid (bubble interior), φ = 1 for heavy fluid (outside)
  // φ = 0.5 + 0.5*tanh(2*(r-R)/W) per Eq.(11)
  T R_phys = Bubble_Radius * Cell_Len;
  T xc_phys = Bubble_Center[0] * Cell_Len;
  T yc_phys = Bubble_Center[1] * Cell_Len;
  T zc_phys = Bubble_Center[2] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    auto& blockLat = PFLattice.getBlockLat(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
          T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
          T z = block.getMin()[2] + static_cast<T>(k) * voxelSize;
          T dx = x - xc_phys;
          T dy = y - yc_phys;
          T dz = z - zc_phys;
          T dist = std::sqrt(dx * dx + dy * dy + dz * dz);
          T phi = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
          blockPhi.get(id) = phi;
        }
      }
    }
  }

  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = PFLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();

    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for(int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          PFCELL cell(id, blockLat);
          T phi = blockPhiField.get(id);

          for (unsigned int i_pop = 0; i_pop < LatSet::q; ++i_pop) {
             T cu = 0;
             cell[i_pop] = latset::w<LatSet>(i_pop) * phi * (T(1) + LatSet::InvCs2 * cu);
          }
        }
      }
    }
  }

  // Initialize PF interface width
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // Initialize NS populations to He-Luo incompressible equilibrium (p=0, u=0)
  // MRTForce now uses incompressible model: rho from DENSITY field (not Σf),
  // u = Σ(c·f) (not divided by rho). So Σf=0 (p=0) no longer causes division
  // by zero. This aligns with Fortran initHydroMacroVars2D (pressure=0, u=0).
  // Actual density enters via DENSITY field (set by FFRhoOmegaUpdate3D from phi).
  Vector<T, LatSet::d> u_zero{T{0}};
  T p_init = T{0};
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(blockid);
    int overlap = 0;
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          NSCELL cell(id, blockLat);
          std::array<T, LatSet::q> feq;
          equilibrium::IncompressibleSecondOrder<NSCELL>::apply(feq, p_init, u_zero);
          for (unsigned int k_pop = 0; k_pop < LatSet::q; ++k_pop) {
            cell[k_pop] = feq[k_pop];
          }
        }
      }
    }
  }

  // ------------------ define BCs ------------------
  // NS: bounceback on all walls
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  // PF: bounceback on all walls
  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // No periodic BCs: all boundaries are no-slip (Safi et al. TC2)

  // ------------------ define tasks ------------------
  // NS task: MRTForce with Guo force (moment-space equivalent of BGKForce)
  using NSBulkTask = tmp::Key_TypePair<BulkFlag, collision::MRTForce<NSCELL, FORCE<T, LatSet::d>>>;
  using NSAllTasks = tmp::TupleWrapper<NSBulkTask>;
  using NSTaskSelector = tmp::TaskSelector<NSAllTasks, std::uint8_t, NSCELL>;

  // PF tasks: FF3D (∇φ + n), FFLaplacian3D (∇²φ), FFChemPotential3D (λ)
  using FFNormalTask =
    tmp::Key_TypePair<BulkFlag, ff::FF3D<PFCELL>>;
  using FFLaplacianTask =
    tmp::Key_TypePair<BulkFlag, ff::FFLaplacian3D<PFCELL>>;
  using FFChemPotTask =
    tmp::Key_TypePair<BulkFlag, ff::FFChemPotential3D<PFCELL>>;
  using FFNormalSelector = TaskSelector<std::uint8_t, PFCELL, FFNormalTask>;
  using FFLaplacianSelector = TaskSelector<std::uint8_t, PFCELL, FFLaplacianTask>;
  using FFChemPotSelector = TaskSelector<std::uint8_t, PFCELL, FFChemPotTask>;

  // PF collision: MRTSource with NORMAL (= unit normal) as Allen-Cahn source (matches Fortran)
  using PFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::MRTSource<
      equilibrium::FirstOrder<PFCELL>,
      NORMAL<T, LatSet::d>,
      true,    // WriteToField
      true>>;  // UseCHRelaxation (Fortran-aligned rtvec)
  using PFAllTasks = tmp::TupleWrapper<PFCollisionTask>;
  using PFCollisionTaskSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  // ---- Coupling tasks (PF → NS) ----
  using STForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension3D<PFCELL, NSCELL>>;
  using STForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  using GravForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFGravityForce3D<PFCELL, NSCELL>>;
  using GravForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravForceTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  using RhoOmegaTask =
    tmp::Key_TypePair<BulkFlag, ff::FFRhoOmegaUpdate3D<PFCELL, NSCELL>>;
  using RhoOmegaTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RhoOmegaTask>;
  BlockLatManagerCoupling RhoOmegaCoupling(PFLattice, NSLattice);

  using PreForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFPreForce3D<PFCELL, NSCELL>>;
  using PreForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PreForceTask>;
  BlockLatManagerCoupling PreForceCoupling(PFLattice, NSLattice);

  using ViscoForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFViscoForce3D<PFCELL, NSCELL>>;
  using ViscoForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, ViscoForceTask>;
  BlockLatManagerCoupling ViscoForceCoupling(PFLattice, NSLattice);

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, LatSet::d>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, LatSet::d>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, LatSet::d>>());
  vtmo::ScalarWriter DensityWriter("Density", NSLattice.getField<DENSITY<T>>());
  vtmo::ScalarWriter PressureWriter("Pressure", NSLattice.getField<PRESSURE<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, LatSet::d>>());

  vtmo::vtmWriter<T, LatSet::d> MainWriter("bubble3d", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VecWriter, PressureWriter, DensityWriter, ForceWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  
  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();

  // Compute initial phi gradients, normal, laplacian, chempot (Fortran initHydroMacroVars)
  PFLattice.template ApplyCellDynamics<FFNormalSelector>(FlagFM);
  PFLattice.template ApplyCellDynamics<FFLaplacianSelector>(FlagFM);
  PFLattice.template ApplyCellDynamics<FFChemPotSelector>(FlagFM);
  PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
  PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);

  // Initial NS rho and omega from phi
  RhoOmegaCoupling.ApplyCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);

  MainWriter.WriteBinary(MainLoopTimer());

  std::cout << "gravity = " << gravity << " kappa = " << Kappa << std::endl;
  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep) {
    // ===== Fortran-aligned: force → PF collision → NS collision → stream → macro =====

    // ---- Phase A: Force setup (Fortran computeForce) ----
    // A1: Update per-cell rho and omega from phi
    RhoOmegaCoupling.ApplyCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);

    // A2: Clear NS FORCE to zero
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, LatSet::d>{T{0}});

    // A3: F_s = λ*∇φ (Fortran: surtenforce = chpoten*nablaphi)
    STCoupling.ApplyCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A4: F_b = -ρ*g (Fortran: bodyforcey = -rho * 1e-4/60)
    GravCoupling.ApplyCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A5: F_p = -(p/3)*Δρ*∇φ (Fortran: preforce from prev computeMacro2D)
    PreForceCoupling.ApplyCellDynamics<PreForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A6: Communicate FORCE to ghost cells (aligned with bubble2d/shearflow2dMRT)
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // ---- Phase B: PF collision (Fortran collisionOrderDF2D) ----
    PFLattice.template ApplyCellDynamics<PFCollisionTaskSelector>(FlagFM);

    // ---- Phase C: NS collision (Fortran collisionDenDF2D) ----
    // C1: Viscous force F_v from non-equilibrium moments
    ViscoForceCoupling.ApplyCellDynamics<ViscoForceTaskSelector>(MainLoopTimer(), FlagFM);
    // C2: MRTForce with total accumulated FORCE (F_s+F_b+F_p+F_v)
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);

    // ---- Phase D: Streaming (Fortran streamOrderDF + streamDenDF) ----
    // D1: PF bounceback + stream
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();
    // D2: NS bounceback + stream
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---- Phase E: Macro update (Fortran: computeOrderMacro + computeMacro + setMacroOrderBC) ----
    // E1: phi from PF pops (Fortran: phi = Σ ddg)
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        auto& blockLat = PFLattice.getBlockLat(blockid);
        int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              PFCELL cell(id, blockLat);
              T phi_new = T(0);
              for (unsigned int k_pop = 0; k_pop < LatSet::q; ++k_pop) {
                phi_new += cell[k_pop];
              }
              if (phi_new < T{0}) phi_new = T{0};
              if (phi_new > T{1}) phi_new = T{1};
              blockPhi.get(id) = phi_new;
            }
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // E2: phi=1 at all no-slip walls
    // Only apply to blocks that actually touch a global domain boundary
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      T Lx_global = T(Ni) * Cell_Len;
      T Ly_global = T(Nj) * Cell_Len;
      T Lz_global = T(Nz) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int nz = block.getNz();
        int overlap = block.getOverlap();
        T minX = block.getMin()[0];
        T maxX = block.getMax()[0];
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        T minZ = block.getMin()[2];
        T maxZ = block.getMax()[2];
        // x = 0 face
        if (minX < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id = k * proj[2] + j * proj[1] + overlap;
              blockPhi.get(id) = T{1};
            }
          }
        }
        // x = Ni-1 face
        if (maxX > Lx_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id = k * proj[2] + j * proj[1] + (nx - 1 - overlap);
              blockPhi.get(id) = T{1};
            }
          }
        }
        // y = 0 face
        if (minY < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id = k * proj[2] + overlap * proj[1] + i;
              blockPhi.get(id) = T{1};
            }
          }
        }
        // y = Nj-1 face
        if (maxY > Ly_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id = k * proj[2] + (ny - 1 - overlap) * proj[1] + i;
              blockPhi.get(id) = T{1};
            }
          }
        }
        // z = 0 face
        if (minZ < Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id = overlap * proj[2] + j * proj[1] + i;
              blockPhi.get(id) = T{1};
            }
          }
        }
        // z = Nz-1 face
        if (maxZ > Lz_global - Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id = (nz - 1 - overlap) * proj[2] + j * proj[1] + i;
              blockPhi.get(id) = T{1};
            }
          }
        }
      }
    }
    // E2a: Re-communicate PHI after wall modification (ghost cells need updated phi)
    PFLattice.getField<PHI<T>>().Communicate();

    // E3: Gradients, normal, laplacian, chempot (Fortran computeOrderMacro)
    PFLattice.template ApplyCellDynamics<FFNormalSelector>(FlagFM);
    PFLattice.template ApplyCellDynamics<FFLaplacianSelector>(FlagFM);
    PFLattice.template ApplyCellDynamics<FFChemPotSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // E3a: Wall grad_phi special handling: copy normal component from first interior cell
    // Only apply to blocks that actually touch a global domain boundary
    {
      auto& gradField = PFLattice.getField<GRAD<T, LatSet::d>>();
      T Lx_global = T(Ni) * Cell_Len;
      T Ly_global = T(Nj) * Cell_Len;
      T Lz_global = T(Nz) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockGrad = gradField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int nz = block.getNz();
        int overlap = block.getOverlap();
        T minX = block.getMin()[0];
        T maxX = block.getMax()[0];
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        T minZ = block.getMin()[2];
        T maxZ = block.getMax()[2];
        // x = 0 face: copy grad[0]
        if (minX < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id0 = k * proj[2] + j * proj[1] + overlap;
              std::size_t id1 = k * proj[2] + j * proj[1] + (overlap + 1);
              blockGrad.get(id0)[0] = blockGrad.get(id1)[0];
            }
          }
        }
        // x = Ni-1 face
        if (maxX > Lx_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id0 = k * proj[2] + j * proj[1] + (nx - 1 - overlap);
              std::size_t id1 = k * proj[2] + j * proj[1] + (nx - 2 - overlap);
              blockGrad.get(id0)[0] = blockGrad.get(id1)[0];
            }
          }
        }
        // y = 0 face: copy grad[1]
        if (minY < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id0 = k * proj[2] + overlap * proj[1] + i;
              std::size_t id1 = k * proj[2] + (overlap + 1) * proj[1] + i;
              blockGrad.get(id0)[1] = blockGrad.get(id1)[1];
            }
          }
        }
        // y = Nj-1 face
        if (maxY > Ly_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id0 = k * proj[2] + (ny - 1 - overlap) * proj[1] + i;
              std::size_t id1 = k * proj[2] + (ny - 2 - overlap) * proj[1] + i;
              blockGrad.get(id0)[1] = blockGrad.get(id1)[1];
            }
          }
        }
        // z = 0 face: copy grad[2]
        if (minZ < Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id0 = overlap * proj[2] + j * proj[1] + i;
              std::size_t id1 = (overlap + 1) * proj[2] + j * proj[1] + i;
              blockGrad.get(id0)[2] = blockGrad.get(id1)[2];
            }
          }
        }
        // z = Nz-1 face
        if (maxZ > Lz_global - Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id0 = (nz - 1 - overlap) * proj[2] + j * proj[1] + i;
              std::size_t id1 = (nz - 2 - overlap) * proj[2] + j * proj[1] + i;
              blockGrad.get(id0)[2] = blockGrad.get(id1)[2];
            }
          }
        }
      }
    }

    // E4: Chempot extrapolation at all walls (Fortran setMacroOrderBC chpoten)
    // Only apply to blocks that actually touch a global domain boundary
    {
      auto& chpotenField = PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      T Lx_global = T(Ni) * Cell_Len;
      T Ly_global = T(Nj) * Cell_Len;
      T Lz_global = T(Nz) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockChpoten = chpotenField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int nz = block.getNz();
        int overlap = block.getOverlap();
        T minX = block.getMin()[0];
        T maxX = block.getMax()[0];
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        T minZ = block.getMin()[2];
        T maxZ = block.getMax()[2];
        // x = 0 face
        if (minX < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id1 = k * proj[2] + j * proj[1] + overlap;
              std::size_t id2 = k * proj[2] + j * proj[1] + (overlap + 1);
              std::size_t id3 = k * proj[2] + j * proj[1] + (overlap + 2);
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
        // x = Ni-1 face
        if (maxX > Lx_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              std::size_t id1 = k * proj[2] + j * proj[1] + (nx - 1 - overlap);
              std::size_t id2 = k * proj[2] + j * proj[1] + (nx - 2 - overlap);
              std::size_t id3 = k * proj[2] + j * proj[1] + (nx - 3 - overlap);
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
        // y = 0 face
        if (minY < Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id1 = k * proj[2] + overlap * proj[1] + i;
              std::size_t id2 = k * proj[2] + (overlap + 1) * proj[1] + i;
              std::size_t id3 = k * proj[2] + (overlap + 2) * proj[1] + i;
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
        // y = Nj-1 face
        if (maxY > Ly_global - Cell_Len * T(1.5)) {
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id1 = k * proj[2] + (ny - 1 - overlap) * proj[1] + i;
              std::size_t id2 = k * proj[2] + (ny - 2 - overlap) * proj[1] + i;
              std::size_t id3 = k * proj[2] + (ny - 3 - overlap) * proj[1] + i;
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
        // z = 0 face
        if (minZ < Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id1 = overlap * proj[2] + j * proj[1] + i;
              std::size_t id2 = (overlap + 1) * proj[2] + j * proj[1] + i;
              std::size_t id3 = (overlap + 2) * proj[2] + j * proj[1] + i;
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
        // z = Nz-1 face
        if (maxZ > Lz_global - Cell_Len * T(1.5)) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              std::size_t id1 = (nz - 1 - overlap) * proj[2] + j * proj[1] + i;
              std::size_t id2 = (nz - 2 - overlap) * proj[2] + j * proj[1] + i;
              std::size_t id3 = (nz - 3 - overlap) * proj[2] + j * proj[1] + i;
              blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
            }
          }
        }
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // E5: NS macro from streamed pops (Fortran computeMacro: p=Σf, u=Σe*f+0.5*F/ρ)
    {
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockLat = NSLattice.getBlockLat(blockid);
        auto& rhoField = NSLattice.getField<DENSITY<T>>();
        auto& presField = NSLattice.getField<PRESSURE<T>>();
        auto& velField = NSLattice.getField<VELOCITY<T, LatSet::d>>();
        auto& forceField = NSLattice.getField<FORCE<T, LatSet::d>>();
        auto& blockRho = rhoField.getBlockField(blockid);
        auto& blockPres = presField.getBlockField(blockid);
        auto& blockVel = velField.getBlockField(blockid);
        auto& blockForce = forceField.getBlockField(blockid);
        int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              NSCELL cell(id, blockLat);
              T pres = T{0};
              Vector<T, LatSet::d> u_raw{T{0}};
              for (unsigned int k_pop = 0; k_pop < LatSet::q; ++k_pop) {
                pres += cell[k_pop];
                for (unsigned int d = 0; d < LatSet::d; ++d) {
                  u_raw[d] += latset::c<LatSet>(k_pop)[d] * cell[k_pop];
                }
              }
              T rho = blockRho.get(id);
              const auto& F = blockForce.get(id);
              Vector<T, LatSet::d> u;
              for (unsigned int d = 0; d < LatSet::d; ++d) {
                u[d] = u_raw[d] + T{0.5} * F[d] / rho;
              }
              blockPres.get(id) = pres;
              blockVel.get(id) = u;
            }
          }
        }
      }
    }

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
      NSLattice.getField<VELOCITY<T, LatSet::d>>().Communicate();
      NSLattice.getField<PRESSURE<T>>().Communicate();
      NSLattice.getField<DENSITY<T>>().Communicate();
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
