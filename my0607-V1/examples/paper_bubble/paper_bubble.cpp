// paper_bubble.cpp
// 2D bubble rising — Section 4.3 (Eo effect on interface dynamics)
// Pressure-based g-LBE + Allen-Cahn, spatially varying viscosity
// Uses fflbm clean operators

#include "freelb.h"
#include "freelb.hh"
#include "fflbm/macroscopic.hh"
#include "fflbm/equilibrium.hh"
#include "fflbm/collision.hh"
#include "fflbm/phasefield.hh"
#include "fflbm/coupling.hh"
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

T Bubble_Radius;
Vector<T, LatSet::d> Bubble_Center;

T Interface_Width;
T Mobility;
T Tau_phi;
T Omega_phi;
T Kappa;
T Beta;

T rho_l;
T rho_h;
T eta_l;
T eta_h;
T sigma;
T gravity;
T rho0;
T Eo;
T Re;
T U_g;
T Tau_ns;
T Omega_ns_init;

int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("paper_bubble.ini");
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

  try {
    eta_l = param_reader.getValue<T>("Two_Phase", "eta_l");
    eta_h = param_reader.getValue<T>("Two_Phase", "eta_h");
  } catch (const std::exception&) {
    T D_fb = T(2.0) * Bubble_Radius;
    T nu_fb = U_g * D_fb / Re;
    eta_h = nu_fb * rho_h;
    eta_l = eta_h / T(10.0);
  }

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ===== 5-step parameter design =====
  T D = T(2.0) * Bubble_Radius;
  rho0 = (rho_l + rho_h) / T{2.0};

  if (U_g < T(0.01) || U_g > T(0.15)) {
    MPI_RANK(0) std::cerr << "[Warning] U_g=" << U_g
                           << " out of [0.01,0.15]\n";
  }

  T g_abs = U_g * U_g / D;
  gravity = -g_abs;
  T nu = U_g * D / Re;
  sigma = (rho_h - rho_l) * g_abs * D * D / Eo;

  Kappa = T(1.5) * Interface_Width * sigma;
  Beta = T(12.0) * sigma / Interface_Width;
  Tau_phi = T(0.5) + Mobility / LatSet::cs2;
  Omega_phi = T(1.0) / Tau_phi;
  Tau_ns = T(0.5) + nu / LatSet::cs2;
  Omega_ns_init = T(1.0) / Tau_ns;

  T cs = std::sqrt(LatSet::cs2);
  T CFL = U_g + cs;
  bool cfl_ok = (CFL < T(1.2));

  MPI_RANK(0) {
    std::cout << "-------- Paper Bubble (g-LBE) ---------\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj
              << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Bubble]: R=" << Bubble_Radius << " D=" << D
              << "  Center=(" << Bubble_Center[0] << ","
              << Bubble_Center[1] << ")\n";
    std::cout << "[Design] rho0=" << rho0
              << "  U_g=" << U_g << " (Ma=" << U_g / cs << ")\n";
    std::cout << "[Design] g=" << gravity << "  nu=" << nu
              << "  sigma=" << sigma << "\n";
    std::cout << "[Visc] eta_l=" << eta_l << " eta_h=" << eta_h
              << "  tau_ns=" << Tau_ns << "  omega=" << Omega_ns_init << "\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " beta=" << Beta << " kappa=" << Kappa
              << " tau_phi=" << Tau_phi << "\n";
    std::cout << "[Flow]: Re=" << Re << "  Eo=" << Eo << "\n";
    std::cout << "[CFL]: " << CFL
              << (cfl_ok ? "  OK" : "  WARNING: >1.2!") << "\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep
              << "  OutputStep=" << OutputStep << "\n";
    std::cout << "---------------------------------------\n";
  }
  if (!cfl_ok) { MPI_RANK(0) std::cerr << "[ERROR] CFL>1.2!\n"; exit(1); }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT
  Printer::Print_BigBanner("Initializing Paper Bubble (g-LBE)...");
  readParam();

  // ------------------ converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);
  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_PaperBubble", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ------------------ NS lattice (g-LBE) ------------------
  using NSFIELDS = TypePack<RHO<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, OMEGA<T>>;
  ValuePack NSInitValues(BaseConv.getLatRhoInit(), Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}}, Omega_ns_init);
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
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
                         T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
      Geo, PFInitValues, PFBaseConv,
      &NSLattice.getField<VELOCITY<T, LatSet::d>>());
  ff::BroadcastAllParams<T>(PFLattice, rho_l, rho_h, eta_l, eta_h,
                            gravity, Beta, Kappa);

  // ---- Initialize PHI field (tanh profile) ----
  T R_phys = Bubble_Radius * Cell_Len;
  T xc_phys = Bubble_Center[0] * Cell_Len;
  T yc_phys = Bubble_Center[1] * Cell_Len;
  T W_phys = Interface_Width * Cell_Len;
  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    auto& bPhi = phiField.getBlockField(bid);
    T vxSize = block.getVoxelSize();
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        T x = block.getMin()[0] + T(i) * vxSize;
        T y = block.getMin()[1] + T(j) * vxSize;
        T dx = x - xc_phys, dy = y - yc_phys;
        T dist = std::sqrt(dx*dx + dy*dy);
        bPhi.get(id) = T(0.5) + T(0.5) * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
      }
    }
  }

  // Initialize PF populations: g_i = w_i * phi
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    auto& bLat = PFLattice.getBlockLat(bid);
    auto& bPhi = phiField.getBlockField(bid);
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        PFCELL cell(id, bLat);
        T phi = bPhi.get(id);
        for (unsigned int k = 0; k < LatSet::q; ++k)
          cell[k] = latset::w<LatSet>(k) * phi;
      }
    }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  // Initialize NS populations: g_i = w_i * p_init/cs²  (pressure-based, u=0)
  T p_init_cs2 = BaseConv.getLatRhoInit();
  for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
    auto& bLat = NSLattice.getBlockLat(bid);
    const auto& block = Geo.getBlock(bid);
    const auto& proj = block.getProjection();
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        std::size_t id = j * proj[1] + i;
        NSCELL cell(id, bLat);
        for (unsigned int k = 0; k < LatSet::q; ++k)
          cell[k] = latset::w<LatSet>(k) * p_init_cs2;
      }
    }
  }

  // Initialize per-cell omega to uniform base value
  NSLattice.getField<OMEGA<T>>().InitValue(Omega_ns_init);

  // set reference density for pressure-based moment (after NSCELL defined)
  fflbm::PressureMoment<NSCELL>::rho0 = rho0;

  // ------------------ BCs ------------------
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
      NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
      PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ fflbm tasks ------------------
  // PF pre-processing: ∇φ (GRAD, NORMAL=n), ∇²φ (LAPLACIAN), μ (CHEMPOT)
  using FFPhiGradTask  = tmp::Key_TypePair<BulkFlag, fflbm::PhiGradient<PFCELL>>;
  using FFLaplacianTask = tmp::Key_TypePair<BulkFlag, fflbm::PhiLaplacian<PFCELL>>;
  using FFChemPotTask  = tmp::Key_TypePair<BulkFlag, fflbm::ChemPotential<PFCELL>>;
  using PFPreCollection = tmp::TupleWrapper<FFPhiGradTask, FFLaplacianTask, FFChemPotTask>;
  using PFPreSelector = tmp::TaskSelector<PFPreCollection, std::uint8_t, PFCELL>;

  // PF collision: BGK + Allen-Cahn source (reads n from NORMAL)
  using PFCollTask = tmp::Key_TypePair<
      BulkFlag, fflbm::BGKSourceCollision<fflbm::FirstOrderEq<PFCELL>>>;
  using PFCollSel   = TaskSelector<std::uint8_t, PFCELL, PFCollTask>;

  // NS collision: pressure-based BGK + Guo force + per-cell omega
  using NSBulkTask = tmp::Key_TypePair<
      BulkFlag,
      fflbm::PressureBGKCollision<fflbm::PressureMoment<NSCELL>,
                                  fflbm::PressureEq<NSCELL>>>;
  using NSTaskSel   = TaskSelector<std::uint8_t, NSCELL, NSBulkTask>;

  // Coupling: surface tension (PF → NS)
  using STForceTask = tmp::Key_TypePair<BulkFlag, fflbm::SurfaceTension<PFCELL, NSCELL>>;
  using STForceSel  = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  // Coupling: gravity buoyancy (PF → NS)
  using GravForceTask = tmp::Key_TypePair<BulkFlag, fflbm::GravityBuoyancy<PFCELL, NSCELL>>;
  using GravForceSel  = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravForceTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  // Coupling: per-cell omega update (PF → NS)
  using OmegaTask = tmp::Key_TypePair<BulkFlag, fflbm::OmegaUpdate<PFCELL, NSCELL>>;
  using OmegaSel  = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, OmegaTask>;
  BlockLatManagerCoupling OmegaCoupling(PFLattice, NSLattice);

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());
  vtmo::ScalarWriter OmegaWriter("Omega", NSLattice.getField<OMEGA<T>>());
  vtmo::vtmWriter<T, 2> MainWriter("paper_bubble", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VecWriter, RhoWriter, ForceWriter, OmegaWriter);

  Timer MainLoopTimer;
  Timer OutputTimer;
  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  MainWriter.WriteBinary(MainLoopTimer());
  std::cout << "gravity=" << gravity << " kappa=" << Kappa
            << " beta=" << Beta << " rho0=" << rho0 << std::endl;
  Printer::Print_BigBanner("Start g-LBE Calculation...");

  // ===================== Main Loop =====================
  while (MainLoopTimer() < MaxStep) {
    // Step 1: Update PHI (φ = Σg_i)
    {
      auto& phiField2 = PFLattice.getField<PHI<T>>();
      for (int bid = 0; bid < Geo.getBlockNum(); ++bid) {
        const auto& block = Geo.getBlock(bid);
        const auto& proj = block.getProjection();
        auto& bPhi = phiField2.getBlockField(bid);
        auto& bLat = PFLattice.getBlockLat(bid);
        int ovl = block.getOverlap();
        for (int j = ovl; j < block.getNy() - ovl; ++j) {
          for (int i = ovl; i < block.getNx() - ovl; ++i) {
            std::size_t id = j * proj[1] + i;
            PFCELL cell(id, bLat);
            T s = T{0};
            for (unsigned int k = 0; k < LatSet::q; ++k) s += cell[k];
            if (s < T{0}) s = T{0};
            if (s > T{1}) s = T{1};
            bPhi.get(id) = s;
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Step 2: ∇φ (GRAD, NORMAL=n), ∇²φ (LAPLACIAN), μ (CHEMICALPOTENTIAL)
    PFLattice.template ApplyCellDynamics<PFPreSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // Step 3: Clear NS force
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    // Step 4: Surface tension  F_s =  μ · ∇φ
    STCoupling.ApplyCellDynamics<STForceSel>(MainLoopTimer(), FlagFM);

    // Step 5: Gravity buoyancy  F_b = (ρ-ρ_h)·g
    GravCoupling.ApplyCellDynamics<GravForceSel>(MainLoopTimer(), FlagFM);

    // Step 6: Update per-cell omega (spatially varying viscosity)
    OmegaCoupling.ApplyCellDynamics<OmegaSel>(MainLoopTimer(), FlagFM);

    // Step 7: NS collision (pressure-based BGK, per-cell ω, Guo force)
    NSLattice.template ApplyCellDynamics<NSTaskSel>(FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // Step 8: NS BCs + Stream + Communicate
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // Step 9: PF collision (Allen-Cahn source via ∇μ)
    PFLattice.template ApplyCellDynamics<PFCollSel>(FlagFM);

    // Step 10: PF BCs + Stream + Communicate
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();

    ++MainLoopTimer; ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
      PFLattice.getField<NORMAL<T, 2>>().Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();
      NSLattice.getField<RHO<T>>().Communicate();
      NSLattice.getField<OMEGA<T>>().Communicate();
      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner("Calculation Complete!");
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}
