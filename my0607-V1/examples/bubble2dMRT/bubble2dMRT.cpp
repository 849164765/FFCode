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
  iniReader param_reader("bubble2dMRT.ini");
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
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  Re = param_reader.getValue<T>("Two_Phase", "Re");
  U_g = param_reader.getValue<T>("Two_Phase", "U_g");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ===== Fortran-aligned hardcoded parameters (BubbleRising.f90 lines 269-289) =====
  // rho_low=0.001, rho_high=1.0, phi_low=0, phi_high=1
  // vis_low=0.6/3500, vis_high=0.6/35
  // surtencoefficent=1e-4*60/125 = 4.8e-5
  // gravity = 1e-4/60 (positive, applied as -rho*g)
  // tau_ddg = 3*mobility+0.5 = 0.53

  eta_l = T(0.6) / T(3500.0);                      // vis_low=0.6/3500 (dynamic)
  eta_h = T(0.6) / T(35.0);                         // vis_high=0.6/35 (dynamic)
  sigma = T(1.0e-4) * T(60.0) / T(125.0);           // surtencoefficent = 4.8e-5
  gravity = T(1.0e-4) / T(60.0);                    // g = 1.667e-6 (positive)

  // Derived: Eq.(14) β=12σ/W, κ=3Wσ/2
  Beta = T(12.0) * sigma / Interface_Width;          // = 1.44e-4
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5); // = 2.88e-4
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
    std::cout << "-------Bubble Rising Simulation (Fortran Aligned)-------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius
              << "  Center=(" << Bubble_Center[0] << "," << Bubble_Center[1] << ")\n";
    std::cout << "[Fortran]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " U_g=" << U_g << " Ma=" << Ma << "\n";
    std::cout << "[Fortran]: eta_l=" << eta_l << " eta_h=" << eta_h
              << " sigma=" << sigma << " g=" << gravity << "\n";
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
  constexpr std::uint8_t PeriodicFlag = std::uint8_t(8);

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
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> left(Vector<T, 2>(T(-Cell_Len), T(0)),
                  Vector<T, 2>(T(0), T(Nj * Cell_Len)));
  AABB<T, 2> right(Vector<T, 2>(T(Ni * Cell_Len), T(0)),
                   Vector<T, 2>(T((Ni + 1) * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(3,3);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.forEach(left, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(right, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // Write flag geometry for verification
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Bubble", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define NS lattice ------------------
  // Note: RHO<T> is cosmetic only (forcerhoU overwrites it with Σf≈1.0 each step)
  // Actual density variation enters via the FORCE field
  using NSFIELDS = TypePack<DENSITY<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, OMEGA<T>, PRESSURE<T>>;
  T omega_ns_init = T{1} / Tau_ns;
  ValuePack NSInitValues(T{1}, Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}}, omega_ns_init, T{});
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
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
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
  T W_phys = Interface_Width * Cell_Len;

  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    auto& blockLat = PFLattice.getBlockLat(blockid);
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

  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = PFLattice.getBlockLat(block_idx);
    auto& blockPhiField = phiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();
    
    for (int j = overlap; j < block.getNy() - overlap; ++j){ 
      for(int i = overlap; i < block.getNx() - overlap; ++i){ 
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhiField.get(id);
        
        for (unsigned int i = 0; i < LatSet::q; ++i) {
           T cu = 0;
           cell[i] = latset::w<LatSet>(i) * phi * (T(1) + LatSet::InvCs2 * cu);
           }
      }
    }
  }

  // Initialize PF interface width
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // Initialize NS populations to incompressible equilibrium
  // p_init = 0 (Fortran: pressure(i,j) = 0.d0)
  Vector<T, 2> u_zero{T{0}, T{0}};
  T p_init = T{0};
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(blockid);
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T u2 = T{0};
        NSCELL cell(id, blockLat);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
          T uc = u_zero * latset::c<LatSet>(k);
          // He-Luo incompressible: f_eq = w_k * (p + 3*(c.u) + 4.5*(c.u)^2 - 1.5*u^2)
          T feq = latset::w<LatSet>(k) * (p_init
                  + LatSet::InvCs2 * uc
                  + uc * uc * T{0.5} * LatSet::InvCs4
                  - LatSet::InvCs2 * u2 * T{0.5});
          cell[k] = feq;
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

  // Periodic BC in X direction
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using FM = BlockFieldManager<FLAG, T, LatSet::d>;

  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_Periodic("NS_Periodic", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_Periodic("PF_Periodic", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PF_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

#ifdef MPI_ENABLED
  NS_Periodic.SetupMPI(GeoHelper);
  PF_Periodic.SetupMPI(GeoHelper);
#endif

  // ------------------ define tasks ------------------
  // NS task: MRTForce with Guo force (moment-space equivalent of BGKForce)
  using NSBulkTask = tmp::Key_TypePair<BulkFlag, collision::MRTForce<NSCELL, FORCE<T, LatSet::d>>>;
  using NSPeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<NSCELL>>;
  using NSAllTasks = tmp::TupleWrapper<NSBulkTask, NSPeriodicTask>;
  using NSTaskSelector = tmp::TaskSelector<NSAllTasks, std::uint8_t, NSCELL>;

  // PF tasks: FF2D (∇φ + n with ε=0.005), FFLaplacian2D (∇²φ), FFChemPotential2D (λ)
  using FFNormalTask =
    tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLaplacianTask =
    tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemPotTask =
    tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
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
  using PFPeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<PFCELL>>;
  using PFAllTasks = tmp::TupleWrapper<PFCollisionTask, PFPeriodicTask>;
  using PFCollisionTaskSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  // ---- Coupling tasks (PF → NS) ----
  using STForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using STForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  using GravForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFGravityForce2D<PFCELL, NSCELL>>;
  using GravForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravForceTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  using RhoOmegaTask =
    tmp::Key_TypePair<BulkFlag, ff::FFRhoOmegaUpdate2D<PFCELL, NSCELL>>;
  using RhoOmegaTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RhoOmegaTask>;
  BlockLatManagerCoupling RhoOmegaCoupling(PFLattice, NSLattice);

  using PreForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFPreForce2D<PFCELL, NSCELL>>;
  using PreForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PreForceTask>;
  BlockLatManagerCoupling PreForceCoupling(PFLattice, NSLattice);

  using ViscoForceTask =
    tmp::Key_TypePair<BulkFlag, ff::FFViscoForce2D<PFCELL, NSCELL>>;
  using ViscoForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, ViscoForceTask>;
  BlockLatManagerCoupling ViscoForceCoupling(PFLattice, NSLattice);

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  //vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  //vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  //vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  //vtmo::ScalarWriter DensityWriter("Density", NSLattice.getField<DENSITY<T>>());
  //vtmo::ScalarWriter PressureWriter("Pressure", NSLattice.getField<PRESSURE<T>>());
  //vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubble2d", Geo);
  MainWriter.addWriterSet(PHIWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  
  PFLattice.NormalFullCommunicate();
  NSLattice.NormalFullCommunicate();
  NS_Periodic.Apply();
  PF_Periodic.Apply();

  // Compute initial phi gradients, normal, laplacian, chempot (Fortran initHydroMacroVars)
  PFLattice.template ApplyInnerCellDynamics<FFNormalSelector>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<FFLaplacianSelector>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<FFChemPotSelector>(FlagFM);
  PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
  PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);

  // Initial NS rho and omega from phi
  RhoOmegaCoupling.ApplyInnerCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);

  MainWriter.WriteBinary(MainLoopTimer());

  std::cout << "gravity = " << gravity << " kappa = " << Kappa << std::endl;
  Printer::Print_BigBanner(std::string("Start Calculation..."));

  while (MainLoopTimer() < MaxStep) {
    // ===== Fortran-aligned: force → PF collision → NS collision → stream → macro =====

    // ---- Phase A: Force setup (Fortran computeForce) ----
    // A1: Update per-cell rho and omega from phi
    RhoOmegaCoupling.ApplyInnerCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);

    // A2: Clear NS FORCE to zero
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    // A3: F_s = λ*∇φ (Fortran: surtenforce = chpoten*nablaphi)
    STCoupling.ApplyInnerCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A4: F_b = -ρ*g (Fortran: bodyforcey = -rho * 1e-4/60)
    GravCoupling.ApplyInnerCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A5: F_p = -(p/3)*Δρ*∇φ (Fortran: preforce from prev computeMacro2D)
    PreForceCoupling.ApplyInnerCellDynamics<PreForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A6: Communicate FORCE to ghost cells (aligned with bubble2d/shearflow2dMRT)
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // ---- Phase B: PF collision (Fortran collisionOrderDF2D) ----
    PFLattice.template ApplyInnerCellDynamics<PFCollisionTaskSelector>(FlagFM);
    PF_Periodic.Apply();
    // Communicate post-collision pops to ghost cells (needed for streaming)
    PFLattice.NormalFullCommunicate();

    // ---- Phase C: NS collision (Fortran collisionDenDF2D) ----
    // C1: Viscous force F_v from non-equilibrium moments
    ViscoForceCoupling.ApplyInnerCellDynamics<ViscoForceTaskSelector>(MainLoopTimer(), FlagFM);
    // C2: MRTForce with total accumulated FORCE (F_s+F_b+F_p+F_v)
    NSLattice.template ApplyInnerCellDynamics<NSTaskSelector>(FlagFM);
    NS_Periodic.Apply();
    // Communicate post-collision pops to ghost cells (needed for streaming)
    NSLattice.NormalFullCommunicate();

    // ---- Phase D: Streaming (Fortran streamOrderDF + streamDenDF) ----
    // D1: PF bounceback + stream
    PF_BB.Apply(MainLoopTimer());
    // D1a: Set PF Y ghost cell POPs to wall cell POPs
    // Y ghost cells at physical walls are not in any SharedComm and are
    // corrupted by pointer rotation (streaming). Their stale POPs are pulled
    // into wall cells during streaming, causing pressure/phi errors that
    // propagate to internal cells and eventually cause NaN at T4100.
    // Fix: copy wall cell POPs to ghost cell POPs (zero-gradient condition).
    {
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockLat = PFLattice.getBlockLat(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = 0 * proj[1] + i;
            std::size_t id_wall = overlap * proj[1] + i;
            PFCELL ghost_cell(id_ghost, blockLat);
            PFCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = (ny - 1) * proj[1] + i;
            std::size_t id_wall = (ny - 1 - overlap) * proj[1] + i;
            PFCELL ghost_cell(id_ghost, blockLat);
            PFCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
      }
    }
    PFLattice.Stream();
    PFLattice.NormalFullCommunicate();
    // D2: NS bounceback + stream
    NS_BB.Apply(MainLoopTimer());
    // D2a: Set NS Y ghost cell POPs to wall cell POPs (same rationale as D1a)
    {
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockLat = NSLattice.getBlockLat(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = 0 * proj[1] + i;
            std::size_t id_wall = overlap * proj[1] + i;
            NSCELL ghost_cell(id_ghost, blockLat);
            NSCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = (ny - 1) * proj[1] + i;
            std::size_t id_wall = (ny - 1 - overlap) * proj[1] + i;
            NSCELL ghost_cell(id_ghost, blockLat);
            NSCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
      }
    }
    NSLattice.Stream();
    NSLattice.NormalFullCommunicate();

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
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            PFCELL cell(id, blockLat);
            T phi_new = T(0);
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              phi_new += cell[k];
            }
            if (phi_new < T{0}) phi_new = T{0};
            if (phi_new > T{1}) phi_new = T{1};
            blockPhi.get(id) = phi_new;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // E2: phi=1 at top/bottom walls (Fortran setMacroOrderBC: phi(i,1)=1, phi(i,ny)=1)
    // Only apply to blocks that actually touch the physical wall (global j=0 or j=Nj-1)
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        int nx = block.getNx();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        // Bottom wall: only if block touches global j=0
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id_bot = overlap * proj[1] + i;
            blockPhi.get(id_bot) = T{1};
          }
        }
        // Top wall: only if block touches global j=Nj-1
        if (maxY > H_global - Cell_Len * T(1.5)) {
          int ny = block.getNy();
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id_top = (ny - 1 - overlap) * proj[1] + i;
            blockPhi.get(id_top) = T{1};
          }
        }
      }
    }
    // E2a: Re-communicate PHI after wall modification (ghost cells need updated phi)
    PFLattice.getField<PHI<T>>().Communicate();

    // E2b: Set PHI at Y ghost cells (physical wall ghost cells are NOT in any
    // SharedComm — they are "三不管" zones: block comm doesn't cover them,
    // BBLikeFixedBlockBdManager only processes internal cells, and no diagonal
    // neighbor can fill them. FF2D::apply reads PHI from all 8 neighbors
    // including these stale ghost cells, causing incorrect gradient at wall
    // cells and eventual NaN.)
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{1};
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{1};
            }
          }
        }
      }
    }

    // E3: Gradients, normal, laplacian, chempot (Fortran computeOrderMacro)
    PFLattice.template ApplyInnerCellDynamics<FFNormalSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFLaplacianSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFChemPotSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // E3a: Wall grad_phi special handling (Fortran: nablaphiy(i,1)=nablaphiy(i,2), nablaphiy(i,ny)=nablaphiy(i,ny-1))
    // Only apply to blocks that actually touch the physical wall
    {
      auto& gradField = PFLattice.getField<GRAD<T, LatSet::d>>();
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockGrad = gradField.getBlockField(blockid);
        int nx = block.getNx();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        // Bottom wall: only if block touches global j=0
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id_bot = overlap * proj[1] + i;
            std::size_t id_bot1 = (overlap + 1) * proj[1] + i;
            blockGrad.get(id_bot)[1] = blockGrad.get(id_bot1)[1];
          }
        }
        // Top wall: only if block touches global j=Nj-1
        if (maxY > H_global - Cell_Len * T(1.5)) {
          int ny = block.getNy();
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id_top = (ny - 1 - overlap) * proj[1] + i;
            std::size_t id_top1 = (ny - 2 - overlap) * proj[1] + i;
            blockGrad.get(id_top)[1] = blockGrad.get(id_top1)[1];
          }
        }
      }
    }

    // E4: Chempot extrapolation at walls (Fortran setMacroOrderBC chpoten)
    // Only apply to blocks that actually touch the physical wall
    {
      auto& chpotenField = PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockChpoten = chpotenField.getBlockField(blockid);
        int nx = block.getNx();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        // Bottom wall: only if block touches global j=0
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id1 = overlap * proj[1] + i;
            std::size_t id2 = (overlap + 1) * proj[1] + i;
            std::size_t id3 = (overlap + 2) * proj[1] + i;
            blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
          }
        }
        // Top wall: only if block touches global j=Nj-1
        if (maxY > H_global - Cell_Len * T(1.5)) {
          int ny = block.getNy();
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id1 = (ny - 1 - overlap) * proj[1] + i;
            std::size_t id2 = (ny - 2 - overlap) * proj[1] + i;
            std::size_t id3 = (ny - 3 - overlap) * proj[1] + i;
            blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
          }
        }
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // E4a: Extrapolate GRAD, NORMAL, CHEMPOT to Y ghost cells
    // (same "三不管" issue as E2b — these fields are not in SharedComm at
    // physical wall ghost cells, so they remain stale after Communicate)
    {
      T H_global = T(Nj) * Cell_Len;
      auto& gradField = PFLattice.getField<GRAD<T, LatSet::d>>();
      auto& normField = PFLattice.getField<NORMAL<T, LatSet::d>>();
      auto& chmField = PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockGrad = gradField.getBlockField(blockid);
        auto& blockNorm = normField.getBlockField(blockid);
        auto& blockChm = chmField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_wall = overlap * proj[1] + i;
              std::size_t id_ghost = j * proj[1] + i;
              blockGrad.get(id_ghost) = blockGrad.get(id_wall);
              blockNorm.get(id_ghost) = blockNorm.get(id_wall);
              blockChm.get(id_ghost) = blockChm.get(id_wall);
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_wall = (ny - 1 - overlap) * proj[1] + i;
              std::size_t id_ghost = j * proj[1] + i;
              blockGrad.get(id_ghost) = blockGrad.get(id_wall);
              blockNorm.get(id_ghost) = blockNorm.get(id_wall);
              blockChm.get(id_ghost) = blockChm.get(id_wall);
            }
          }
        }
      }
    }

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
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            NSCELL cell(id, blockLat);
            T pres = T{0};
            T ux_raw = T{0};
            T uy_raw = T{0};
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              pres += cell[k];
              ux_raw += latset::c<LatSet>(k)[0] * cell[k];
              uy_raw += latset::c<LatSet>(k)[1] * cell[k];
            }
            T rho = blockRho.get(id);
            const auto& F = blockForce.get(id);
            T ux = ux_raw + T{0.5} * F[0] / rho;
            T uy = uy_raw + T{0.5} * F[1] / rho;
            blockPres.get(id) = pres;
            blockVel.get(id) = Vector<T, 2>{ux, uy};
          }
        }
      }
    }
    // E5a: Communicate NS macro fields to ghost cells every step
    // (PF collision references NS velocity; PreForce references PRESSURE;
    //  ViscoForce/MRTForce reference DENSITY and OMEGA; stale ghost values
    //  cause feedback loop leading to exponential growth and NaN)
    NSLattice.getField<VELOCITY<T, LatSet::d>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    // E5b: Extrapolate NS macro fields to Y ghost cells
    // (same "三不管" issue — VELOCITY is read from neighbors by bounceback BC,
    // and all NS fields need correct ghost values for coupling calculations)
    {
      T H_global = T(Nj) * Cell_Len;
      auto& velField = NSLattice.getField<VELOCITY<T, LatSet::d>>();
      auto& presField = NSLattice.getField<PRESSURE<T>>();
      auto& rhoField = NSLattice.getField<DENSITY<T>>();
      auto& omegaField = NSLattice.getField<OMEGA<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockVel = velField.getBlockField(blockid);
        auto& blockPres = presField.getBlockField(blockid);
        auto& blockRho = rhoField.getBlockField(blockid);
        auto& blockOmega = omegaField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_wall = overlap * proj[1] + i;
              std::size_t id_ghost = j * proj[1] + i;
              blockVel.get(id_ghost) = blockVel.get(id_wall);
              blockPres.get(id_ghost) = blockPres.get(id_wall);
              blockRho.get(id_ghost) = blockRho.get(id_wall);
              blockOmega.get(id_ghost) = blockOmega.get(id_wall);
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_wall = (ny - 1 - overlap) * proj[1] + i;
              std::size_t id_ghost = j * proj[1] + i;
              blockVel.get(id_ghost) = blockVel.get(id_wall);
              blockPres.get(id_ghost) = blockPres.get(id_wall);
              blockRho.get(id_ghost) = blockRho.get(id_wall);
              blockOmega.get(id_ghost) = blockOmega.get(id_wall);
            }
          }
        }
      }
    }

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
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
