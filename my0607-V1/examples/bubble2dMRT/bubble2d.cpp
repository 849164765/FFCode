// bubble2d.cpp
// 2D bubble rising simulation - Phase field + NS flow coupling
// MRT version (Multi-Relaxation-Time)
// No magnetic field, no AMR, MPI parallel

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

using T = FLOAT;
using LatSet = D2Q9<T>;

// ===================================================================
// MRT collision with Guo force scheme for Navier-Stokes
// Equivalent to BGKForce in moment space
// ===================================================================
namespace collision {
template <typename CELL>
struct MRTForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    const T omega = cell.getOmega();
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();

    // 1. Compute rho and u_raw from populations
    T rho{};
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      rho += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
      }
    }
    for (unsigned int d = 0; d < LatSet::d; ++d) u[d] /= rho;

    // 2. Half-force correction: u_c = u + F/(2*rho)
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d) uc[d] += F[d] / (T{2} * rho);

    // Write macroscopic fields (corrected velocity)
    cell.template get<RHO<T>>() = rho;
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;

    // 3. NS relaxation time vector
    T rtvec[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) rtvec[i] = mrt::s<LatSet>(i);
    for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i)
      rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega;

    // 4. Moments from population: m_k = M * f
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];

    // 5. Equilibrium in population space, then to moment space
    std::array<T, LatSet::q> feq{};
    equilibrium::SecondOrder<CELL>::apply(feq, rho, uc);
    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];

    // 6. Force source in population space (uses corrected u, matching force::Force)
    // Fi_i = w_i * [F·c_i/cs^2 - u·F/cs^2 + (u·c_i)(F·c_i)/cs^4]
    T uF = T{};
    for (unsigned int d = 0; d < LatSet::d; ++d) uF += uc[d] * F[d];
    T fi_force[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cF = T{}, cu = T{};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cF += latset::c<LatSet>(i)[d] * F[d];
        cu += latset::c<LatSet>(i)[d] * uc[d];
      }
      fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                     + LatSet::InvCs4 * cu * cF);
    }

    // 7. Force source in moment space: F_k = M * fi
    T force_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j)
      for (unsigned int k = 0; k < LatSet::q; ++k)
        force_m[j] += mrt::M<LatSet>(j, k) * fi_force[k];

    // 8. Combined MRT collision + force source:
    //    f_i = f_i - M^{-1}_ik * s_k * (m_k - m_k^{eq})
    //         + M^{-1}_ik * (1 - s_k/2) * F_k
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll += invM_ij * rtvec[j] * (momenta[j] - momentaEq[j]);
        source += invM_ij * (T{1} - rtvec[j] / T{2}) * force_m[j];
      }
      cell[i] = cell[i] - coll + source;
    }
  }
};
} // namespace collision

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
T U_g;
T Tau_ns;

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
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  Re = param_reader.getValue<T>("Two_Phase", "Re");
  U_g = param_reader.getValue<T>("Two_Phase", "U_g");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  T D = T(2.0) * Bubble_Radius;
  if (U_g < T(0.01) || U_g > T(0.15)) {
    if (mpi().getRank() == 0) {
      std::cerr << "[Warning] U_g=" << U_g << " out of [0.01, 0.15] range\n";
    }
  }

  T g_abs = U_g * U_g / D;
  gravity = -g_abs;

  T nu = U_g * D / Re;
  eta_h = nu * rho_h;
  eta_l = eta_h / T(10);

  T DeltaRho = rho_h - rho_l;
  sigma = DeltaRho * g_abs * D * D / Eo;

  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  Beta = T(12.0) * sigma / Interface_Width;
  Tau_phi = T(0.5) + Mobility / LatSet::cs2;
  Omega_phi = T(1.0) / Tau_phi;

  Tau_ns = T(0.5) + nu / LatSet::cs2;

  T cs = std::sqrt(LatSet::cs2);
  T CFL = U_g + cs;
  bool cfl_ok = (CFL < T(1.2));

  MPI_RANK(0) {
    std::cout << "-------Bubble Rising Simulation (MRT)-------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius << " D=" << D
              << "  Center=(" << Bubble_Center[0] << "," << Bubble_Center[1] << ")\n";
    std::cout << "[Design] Step1: D = " << D << "\n";
    std::cout << "[Design] Step2: U_g = " << U_g << " (Ma=" << U_g/cs << ")\n";
    std::cout << "[Design] Step3: g = U_g^2/D = " << gravity << "\n";
    std::cout << "[Design] Step4: nu = U_g*D/Re = " << nu
              << "  eta_h = " << eta_h << "  eta_l = " << eta_l << "\n";
    std::cout << "[Design] Step5: sigma = DeltaRho*g*D^2/Eo = " << sigma << "\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " beta=" << Beta << " kappa=" << Kappa
              << " tau_phi=" << Tau_phi << "\n";
    std::cout << "[Flow]: tau_ns=" << Tau_ns << "  omega=" << (T(1)/Tau_ns)
              << "  Re=" << Re << "  Eo=" << Eo << "\n";
    std::cout << "[CFL]: (U_g+cs)*dt/dx = " << CFL
              << (cfl_ok ? "  OK" : "  WARNING: >1.2!") << "\n";
    std::cout << "[MRT]: NS=MRTForce  PF=MRTSource\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep << "  OutputStep=" << OutputStep << "\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\n";
#endif
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

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_BubbleMRT", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ define NS lattice ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ------------------ define PF lattice ------------------
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                            ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
                         T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PFBaseConv, &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  ff::BroadcastAllParams<T>(PFLattice,
                            rho_l, rho_h, eta_l, eta_h,
                            gravity, Beta, Kappa);

  // ---- Initialize PHI and POP fields ----
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

    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for(int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, blockLat);
        T phi = blockPhiField.get(id);
        for (unsigned int k = 0; k < LatSet::q; ++k) {
           T cu = 0;
           cell[k] = latset::w<LatSet>(k) * phi * (T(1) + LatSet::InvCs2 * cu);
        }
      }
    }
  }

  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // Initialize NS populations to equilibrium
  T ns_rho_init = BaseConv.getLatRhoInit();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockLat = NSLattice.getBlockLat(blockid);
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

  // ------------------ define BCs ------------------
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ define tasks ------------------
  // NS task: MRT collision with Guo force scheme
  using NSBulkTask = tmp::Key_TypePair<
    BulkFlag, collision::MRTForce<NSCELL>>;
  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask>;

  // PF tasks: FF2D, FFLaplacian2D, FFChemPotential2D
  using FFNormalTask =
    tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLaplacianTask =
    tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemPotTask =
    tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using PFFFTaskCollection = tmp::TupleWrapper<FFNormalTask, FFLaplacianTask, FFChemPotTask>;
  using PFFFTaskSelector = tmp::TaskSelector<PFFFTaskCollection, std::uint8_t, PFCELL>;

  // PF collision: MRT collision with source term (Allen-Cahn via NORMAL)
  using PFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::MRTSource<
      equilibrium::SecondOrder<PFCELL>,
      NORMAL<T, LatSet::d>,
      true>>;
  using PFCollisionTaskSelector = TaskSelector<std::uint8_t, PFCELL, PFCollisionTask>;

  // ---- Coupling tasks (PF -> NS) ----
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

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubble2dMRT", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VecWriter, RhoWriter, ForceWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  MainWriter.WriteBinary(MainLoopTimer());

  std::cout << "gravity = " << gravity << " kappa = " << Kappa << std::endl;
  Printer::Print_BigBanner(std::string("Start Bubble Rising MRT Calculation..."));

  while (MainLoopTimer() < MaxStep) {
    // ---- Phase field pre-processing ----
    // Step 1: Update PHI by summing populations
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

    // Step 2: Compute GRAD, NORMAL, LAPLACIAN, CHEMICALPOTENTIAL
    PFLattice.template ApplyCellDynamics<PFFFTaskSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // ---- NS force accumulation ----
    // Step 3: Clear NS FORCE
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    // Step 4: Surface tension force F_s = lambda*grad-phi -> NS FORCE
    STCoupling.ApplyCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);

    // Step 5: Gravity/buoyancy force F_b -> NS FORCE
    GravCoupling.ApplyCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);

    // ---- NS MRT collision + streaming ----
    // Step 6: NS MRTForce collision (reads FORCE from same cell)
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();
    // Step 7: NS BCs + Stream + Communicate
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---- PF MRT collision + streaming ----
    // Step 8: PF MRTSource collision (reads NORMAL, VELOCITY from NS ref)
    PFLattice.template ApplyCellDynamics<PFCollisionTaskSelector>(FlagFM);

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
      NSLattice.getField<RHO<T>>().Communicate();

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
