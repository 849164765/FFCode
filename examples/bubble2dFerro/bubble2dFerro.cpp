// bubble2dFerro.cpp
// 2D bubble rising in ferrofluid under uniform magnetic field
// Phase field + NS flow + magnetic field (three-way coupling)
// D2Q5 lattice for magnetic potential, D2Q9 for PF/NS
// No AMR, MPI parallel
//
// Diagnostic flags:
//   DIAGNOSTIC_FREEZE_PHI — freeze phase field, output u and F_m in solvent
//   DIAGNOSTIC_H_LINE     — output H field (Hx,Hy) along centerline

// ========== DIAGNOSTIC FLAGS ==========
// #define DIAGNOSTIC_FREEZE_PHI
// #define DIAGNOSTIC_H_LINE
// ======================================

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

#ifdef DIAGNOSTIC_FREEZE_PHI
#include <fstream>
#include <iomanip>
#endif

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;

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

// magnetic field parameters
T mu_l;        // permeability of ferrofluid (light fluid)
T mu_h;        // permeability of non-magnetic fluid (heavy fluid)
T chi_l;       // magnetic susceptibility of ferrofluid
T chi_h;       // magnetic susceptibility of non-magnetic fluid (=0)
T H0_x;        // applied field x-component (lattice units)
T H0_y;        // applied field y-component (lattice units)
int mf_substeps;  // MF pseudo-time sub-iterations per NS step
T epsilon;        // pseudo-time relaxation parameter
T Tau_mf;         // MRT relaxation parameter (max relaxation for convergence)

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("bubble2dFerro.ini");
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

  mu_l = param_reader.getValue<T>("Magnetic_Field", "mu_l");
  mu_h = param_reader.getValue<T>("Magnetic_Field", "mu_h");
  chi_l = param_reader.getValue<T>("Magnetic_Field", "chi_l");
  chi_h = param_reader.getValue<T>("Magnetic_Field", "chi_h");
  H0_x = param_reader.getValue<T>("Magnetic_Field", "H0_x");
  H0_y = param_reader.getValue<T>("Magnetic_Field", "H0_y");
  mf_substeps = param_reader.getValue<int>("Magnetic_Field", "mf_substeps");
  epsilon = param_reader.getValue<T>("Magnetic_Field", "epsilon");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ===== 5-step parameter design =====
  T D = T(2.0) * Bubble_Radius;

  U_g = param_reader.getValue<T>("Two_Phase", "U_g");
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

  // MF relaxation: τ = 0.5 + μ_avg / (ε · cs²)
  // D2Q5 MRT diffusivity: D = cs²·(τ-0.5) = μ_avg/ε for magnetic potential eq.
  // With μ_avg=2, ε=1, cs²=1/3 → τ=6.5
  T mu_est = (mu_l + mu_h) * T(0.5);
  Tau_mf = T(0.5) + mu_est / (epsilon * MFLatSet::cs2);
  if (Tau_mf < T(0.51)) Tau_mf = T(0.51);  // clamp for stability

  T cs = std::sqrt(LatSet::cs2);
  T CFL = U_g + cs;
  bool cfl_ok = (CFL < T(1.2));

  MPI_RANK(0) {
    std::cout << "----------Bubble Rising in Ferrofluid Simulation----------\n";
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
    std::cout << "[Magnetic]: mu_l=" << mu_l << " mu_h=" << mu_h
              << " chi_l=" << chi_l << " chi_h=" << chi_h
              << " H0=(" << H0_x << "," << H0_y << ")\n";
    std::cout << "[MF Solver]: tau_mf=" << Tau_mf << " epsilon=" << epsilon
              << " sub-steps=" << mf_substeps << "\n";
    std::cout << "[CFL]: (U_g+cs)*dt/dx = " << CFL
              << (cfl_ok ? "  OK" : "  WARNING: >1.2!") << "\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep << "  OutputStep=" << OutputStep << "\n";
    std::cout << "--------------------------------------------------------\n";
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

  Printer::Print_BigBanner(std::string("Initializing Bubble Rising in Ferrofluid..."));

  readParam();

  // ------------------ define converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);

  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_mf);

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
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_BubbleFerro", Geo, 1);
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
  // Extended with MU_L, MU_H block constants and MU per-cell field for permeability
  using PFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                            NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                            ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                            ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                            ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>,
                            ff::MU_L<T>, ff::MU_H<T>, ff::MU<T>>;
  using PFFIELDREFS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFFIELDS, PFFIELDREFS>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
                         T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h,
                         mu_l, mu_h, T{});
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PFBaseConv, &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  ff::BroadcastAllParams<T>(PFLattice,
                            rho_l, rho_h, eta_l, eta_h,
                            gravity, Beta, Kappa);
  // Set magnetic block constants on PF lattice
  PFLattice.getField<ff::MU_L<T>>().InitValue(mu_l);
  PFLattice.getField<ff::MU_H<T>>().InitValue(mu_h);

  // ------------------ define MF lattice (D2Q5) ------------------
  using MFFIELDS = TypePack<ff::PSI<T>, POP<T, MFLatSet::q>,
                            ff::H_FIELD<T, MFLatSet::d>,
                            ff::H_SQ<T>, ff::CHI<T>,
                            ff::CHI_L<T>, ff::CHI_H<T>, OMEGA<T>,
                            ff::F_MAGNETIC<T, 2>>;
  T omega_mf_init = T{1} / Tau_mf;
  ValuePack MFInitValues(T{}, T{}, Vector<T, MFLatSet::d>{T{0}, T{0}}, T{}, T{}, chi_l, chi_h,
                         omega_mf_init,
                         Vector<T, 2>{T{0}, T{0}});
  using MFCELL = Cell<T, MFLatSet, MFFIELDS>;
  BlockLatticeManager<T, MFLatSet, MFFIELDS> MFLattice(Geo, MFInitValues, MFBaseConv);

  // Broadcast magnetic params to MF lattice blocks
  ff::BroadcastMagParams<T>(MFLattice, chi_l, chi_h);

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
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
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

#ifdef DIAGNOSTIC_FREEZE_PHI
  // Save snapshot of initial phi for diagnostic freeze
  IF_MPI_RANK(0) {
    std::cout << "[DIAGNOSTIC] Freezing phase field at initial state\n";
  }
  std::vector<std::vector<T>> frozenPhi(Geo.getBlockNum());
  {
    auto& phiField = PFLattice.getField<PHI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      const auto& blockPhi = phiField.getBlockField(blockid);
      int nx = block.getNx();
      int ny = block.getNy();
      frozenPhi[blockid].resize(static_cast<std::size_t>(nx) * ny);
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          std::size_t id = j * proj[1] + i;
          frozenPhi[blockid][id] = blockPhi.get(id);
        }
      }
    }
  }
#endif

  // Initialize NS populations to equilibrium
  Vector<T, 2> u_zero{T{0}, T{0}};
  T ns_rho_init = BaseConv.getLatRhoInit();
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
          T feq = latset::w<LatSet>(k) * ns_rho_init *
                  (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
                   LatSet::InvCs2 * u2 * T{0.5});
          cell[k] = feq;
        }
      }
    }
  }

  // Initialize MF populations to equilibrium (h_α^eq = w_eq_α * ψ)
  // Initial ψ = -H0_y * y (linear variation for uniform vertical field)
  // Use D2Q5 equilibrium weights w_eq = 0.2 (paper Eq. 40)
  constexpr T w_eq[MFLatSet::q] = {T(0.2), T(0.2), T(0.2), T(0.2), T(0.2)};
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockLat = MFLattice.getBlockLat(blockid);
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T y = block.getMin()[1] + static_cast<T>(j) * block.getVoxelSize();
        T psi_init = -H0_y * y;
        MFCELL cell(id, blockLat);
        for (unsigned int k = 0; k < MFLatSet::q; ++k) {
          cell[k] = w_eq[k] * psi_init;
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

  using MFBLKLAT = BlockLatticeManager<T, MFLatSet, MFFIELDS>;
  BBLikeFixedBlockBdManager<bounceback::normal<MFCELL>, MFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    MF_BB("MF_BB", MFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // ------------------ define tasks ------------------
  // NS task: BGKForce with Force vector field
  using NSBulkTask = tmp::Key_TypePair<
    BulkFlag,
    collision::BGKForce<
      moment::forcerhoU<NSCELL, force::Force<NSCELL>, true>,
      equilibrium::SecondOrder<NSCELL>,
      force::Force<NSCELL>>>;
  using NSTaskSelector = TaskSelector<std::uint8_t, NSCELL, NSBulkTask>;

  // PF tasks
  using FFNormalTask =
    tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLaplacianTask =
    tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemPotTask =
    tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using FFNormalSelector = TaskSelector<std::uint8_t, PFCELL, FFNormalTask>;
  using FFLaplacianSelector = TaskSelector<std::uint8_t, PFCELL, FFLaplacianTask>;
  using FFChemPotSelector = TaskSelector<std::uint8_t, PFCELL, FFChemPotTask>;

  using FFChemPotGradTask =
    tmp::Key_TypePair<BulkFlag, ff::FFChemPotentialGradient2D<PFCELL>>;
  using FFChemPotGradSel = TaskSelector<std::uint8_t, PFCELL, FFChemPotGradTask>;

  using PFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::BGKSource<
      equilibrium::SecondOrder<PFCELL>,
      NORMAL<T, LatSet::d>,
      true>>;
  using PFCollisionTaskSelector = TaskSelector<std::uint8_t, PFCELL, PFCollisionTask>;

  // PF permeability update task (μ)
  using MuUpdateTask =
    tmp::Key_TypePair<BulkFlag, ff::FFMuUpdate2D<PFCELL>>;
  using MuUpdateSelector = TaskSelector<std::uint8_t, PFCELL, MuUpdateTask>;

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

  // ---- Coupling task (MF → NS) ----
  using MFForceTask =
    tmp::Key_TypePair<BulkFlag, ff::MFForce2D<MFCELL, NSCELL>>;
  using MFForceTaskSelector =
    CoupledTaskSelector<std::uint8_t, MFCELL, NSCELL, MFForceTask>;
  BlockLatManagerCoupling MFCoupling(MFLattice, NSLattice);

  // MF internal tasks (on MF lattice)
  using MFGradTask =
    tmp::Key_TypePair<BulkFlag, ff::MFGradient2D<MFCELL>>;
  using MFGradSelector = TaskSelector<std::uint8_t, MFCELL, MFGradTask>;

  using MFHsqTask =
    tmp::Key_TypePair<BulkFlag, ff::MFHsq2D<MFCELL>>;
  using MFHsqSelector = TaskSelector<std::uint8_t, MFCELL, MFHsqTask>;

  // MF collision: MRTMag on D2Q5
  // M_5 moments: (ψ, J_x, J_y, e, p_xx)
  // rtvec = (1, ω_mag, ω_mag, 1, 1)
  // Higher-order moments fully relaxed (s=1), flux moments relaxed at diffusion rate ω_mag
  using MFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::MRTMag<MFCELL, true>>;
  using MFCollisionTaskSelector = TaskSelector<std::uint8_t, MFCELL, MFCollisionTask>;

  // ------------------ writers ------------------
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter GRADWriter("GRAD", PFLattice.getField<GRAD<T, 2>>());
  vtmo::VectorWriter NormalWriter("NORMAL", PFLattice.getField<NORMAL<T, 2>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::ScalarWriter RhoWriter("Rho", NSLattice.getField<RHO<T>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());
  vtmo::ScalarWriter PSIWriter("PSI", MFLattice.getField<ff::PSI<T>>());
  vtmo::VectorWriter HWriter("H", MFLattice.getField<ff::H_FIELD<T, 2>>());
  vtmo::ScalarWriter HSQWriter("H_SQ", MFLattice.getField<ff::H_SQ<T>>());
  vtmo::VectorWriter FmWriter("F_magnetic", MFLattice.getField<ff::F_MAGNETIC<T, 2>>());

  vtmo::vtmWriter<T, 2> MainWriter("bubble2dFerro", Geo);
  MainWriter.addWriterSet(PHIWriter, GRADWriter, NormalWriter,
                          VecWriter, RhoWriter, ForceWriter,
                          PSIWriter, HWriter, HSQWriter, FmWriter);

  // ------------------ timer ------------------
  Timer MainLoopTimer;
  Timer OutputTimer;

  PFLattice.NormalCommunicate();
  NSLattice.NormalCommunicate();
  MFLattice.NormalCommunicate();
  MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

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

#ifdef DIAGNOSTIC_FREEZE_PHI
    // Restore phi from frozen snapshot — locks the phase field
    {
      auto& phiField = PFLattice.getField<PHI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        for (int j = 0; j < block.getNy(); ++j) {
          for (int i = 0; i < block.getNx(); ++i) {
            std::size_t id = j * proj[1] + i;
            blockPhi.get(id) = frozenPhi[blockid][id];
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();
#endif

    // Step 2: Compute GRAD, NORMAL, LAPLACIAN, CHEMICALPOTENTIAL
    PFLattice.template ApplyCellDynamics<FFNormalSelector>(FlagFM);
    PFLattice.template ApplyCellDynamics<FFLaplacianSelector>(FlagFM);
    PFLattice.template ApplyCellDynamics<FFChemPotSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // Step 3: Update permeability μ = μ_l + (μ_h-μ_l)*h(φ) → PF  (quintic smooth-step)
    PFLattice.template ApplyCellDynamics<MuUpdateSelector>(FlagFM);

    // ---- PF→MF: compute CHI and OMEGA per cell ----
    // Use quintic smooth-step h(φ)=φ³(6φ²-15φ+10) for stable interpolation
    {
      auto& pfPhiField = PFLattice.getField<PHI<T>>();
      auto& mfChiField = MFLattice.getField<ff::CHI<T>>();
      auto& mfOmegaField = MFLattice.getField<OMEGA<T>>();
      T invCs2 = MFLatSet::cs2;  // 1/3
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& pfBlockPhi = pfPhiField.getBlockField(blockid);
        auto& mfBlockChi = mfChiField.getBlockField(blockid);
        auto& mfBlockOmega = mfOmegaField.getBlockField(blockid);
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            T phi = pfBlockPhi.get(id);
            T h = ff::quinticSmoothStep(phi);
            mfBlockChi.get(id) = chi_l + (chi_h - chi_l) * h;
            // Per-cell omega for variable-coefficient diffusion ∇·(μ∇ψ)
            T mu_local = mu_l + (mu_h - mu_l) * h;
            T tau_local = T(0.5) + mu_local / (epsilon * invCs2);
            mfBlockOmega.get(id) = T(1) / tau_local;
          }
        }
      }
    }
    MFLattice.getField<ff::CHI<T>>().Communicate();
    MFLattice.getField<OMEGA<T>>().Communicate();

    // ---- MF pseudo-time solver ----
    // Step 4: MF collision + stream sub-iterations
    // First update PSI from MF POP for initial condition
    {
      auto& psiField = MFLattice.getField<ff::PSI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPsi = psiField.getBlockField(blockid);
        auto& blockLat = MFLattice.getBlockLat(blockid);
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            MFCELL cell(id, blockLat);
            T psi = T(0);
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              psi += cell[k];
            }
            blockPsi.get(id) = psi;
          }
        }
      }
      psiField.Communicate();
    }

    for (int mf_iter = 0; mf_iter < mf_substeps; ++mf_iter) {
      // MRT collision on D2Q5 for magnetic potential equation
      MFLattice.template ApplyCellDynamics<MFCollisionTaskSelector>(FlagFM);

      // MF BCs + Stream + Communicate
      MF_BB.Apply(MainLoopTimer());
      MFLattice.Stream();
      MFLattice.NormalCommunicate();

      // Update PSI = Σ h_α after streaming
      auto& psiField = MFLattice.getField<ff::PSI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPsi = psiField.getBlockField(blockid);
        auto& blockLat = MFLattice.getBlockLat(blockid);
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            MFCELL cell(id, blockLat);
            T psi = T(0);
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              psi += cell[k];
            }
            blockPsi.get(id) = psi;
          }
        }
      }
      psiField.Communicate();

      // Fix PSI at domain boundary cells using Neumann BC.
      // BounceBack alone bounces populations but gives wrong ψ at walls.
      // Without this, H = -∇ψ produces spurious boundary gradients that
      // propagate via CommunicateHSQ into adjacent blocks.
      //
      // P0 fix: after correcting PSI, re-equilibrate MF populations to match.
      // Without this, the next MRT collision sees inconsistent populations
      // (sum ≠ fixed ψ) and generates spurious non-equilibrium moments that
      // create lattice-scale noise propagating upward through the domain.
      {
        auto& psiFixField = MFLattice.getField<ff::PSI<T>>();
        constexpr T w_eq = T{1} / T{MFLatSet::q};  // = 0.2 for D2Q5
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          T dx = block.getVoxelSize();
          auto& blockPsi = psiFixField.getBlockField(blockid);
          auto& flagBlock = FlagFM.getBlockField(blockid);
          auto& blockLat = MFLattice.getBlockLat(blockid);
          int overlap = block.getOverlap();
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = j * proj[1] + i;
              auto flag = flagBlock.get(id);
              if (flag == BulkFlag) continue;

              // Neumann: ∂ψ/∂n = -H_n. Sample interior psi, then extrapolate.
              T psi_set = blockPsi.get(id);
              bool fixed = false;

              // --- y-direction ---
              if (j + 1 < block.getNy() - overlap) {
                std::size_t id_up = (j + 1) * proj[1] + i;
                if (flagBlock.get(id_up) == BulkFlag) {
                  // Bottom boundary: interior above, (ψ_up - ψ_bd)/dx = -H0_y
                  // → ψ_bd = ψ_up + H0_y·dx
                  psi_set = blockPsi.get(id_up) + H0_y * dx;
                  fixed = true;
                }
              }
              if (!fixed && j > overlap) {
                std::size_t id_dn = (j - 1) * proj[1] + i;
                if (flagBlock.get(id_dn) == BulkFlag) {
                  // Top boundary: interior below, (ψ_bd - ψ_dn)/dx = -H0_y
                  // → ψ_bd = ψ_dn - H0_y·dx
                  psi_set = blockPsi.get(id_dn) - H0_y * dx;
                  fixed = true;
                }
              }

              // --- x-direction (only if y-direction didn't fix) ---
              if (!fixed && i + 1 < block.getNx() - overlap) {
                std::size_t id_rt = j * proj[1] + (i + 1);
                if (flagBlock.get(id_rt) == BulkFlag) {
                  // Left boundary: interior right, (ψ_rt - ψ_bd)/dx = -H0_x
                  // → ψ_bd = ψ_rt + H0_x·dx
                  psi_set = blockPsi.get(id_rt) + H0_x * dx;
                  fixed = true;
                }
              }
              if (!fixed && i > overlap) {
                std::size_t id_lt = j * proj[1] + (i - 1);
                if (flagBlock.get(id_lt) == BulkFlag) {
                  // Right boundary: interior left, (ψ_bd - ψ_lt)/dx = -H0_x
                  // → ψ_bd = ψ_lt - H0_x·dx
                  psi_set = blockPsi.get(id_lt) - H0_x * dx;
                  fixed = true;
                }
              }

              if (fixed) {
                blockPsi.get(id) = psi_set;
                // Re-equilibrate MF populations to match the corrected ψ.
                // Uses MagEquilibrium: h_k^eq = w_eq * ψ with w_eq = 1/q = 0.2
                MFCELL cell(id, blockLat);
                for (unsigned int k = 0; k < MFLatSet::q; ++k) {
                  cell[k] = w_eq * psi_set;
                }
              }
            }
          }
        }
      }
    }

    // ---- Final MF sync after sub-iterations (multi-process fix) ----
    // The last sub-iteration's Neumann BC fix modifies boundary PSI and POP
    // AFTER the last psiField.Communicate().  Without this final sync, ghost
    // cells on neighbour processes carry stale PSI values, causing H = -∇ψ
    // to compute wrong gradients at process boundaries.  The error is largest
    // at the bottom boundary where H0_y * dx correction is maximized, which
    // matches the observed bottom dent in the ferrofluid droplet shape.
    {
      auto& psiField = MFLattice.getField<ff::PSI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPsi = psiField.getBlockField(blockid);
        auto& blockLat = MFLattice.getBlockLat(blockid);
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            MFCELL cell(id, blockLat);
            T psi = T(0);
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              psi += cell[k];
            }
            blockPsi.get(id) = psi;
          }
        }
      }
      psiField.Communicate();
    }

    // Step 5: Compute H = -∇ψ and |H|²
    MFLattice.template ApplyCellDynamics<MFGradSelector>(FlagFM);
    MFLattice.getField<ff::H_FIELD<T, MFLatSet::d>>().Communicate();
    MFLattice.template ApplyCellDynamics<MFHsqSelector>(FlagFM);
    ff::CommunicateHSQ<T>(MFLattice);

    // ---- NS force accumulation ----
    // Step 6: Clear NS FORCE
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    // Step 7: Surface tension force F_s = λ*∇φ → NS FORCE
    STCoupling.ApplyCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);

    // Step 8: Gravity/buoyancy force F_b → NS FORCE
    GravCoupling.ApplyCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);

    // Step 9: Magnetic force F_m = (χ/2)∇(|H|²) → NS FORCE
    MFCoupling.ApplyCellDynamics<MFForceTaskSelector>(MainLoopTimer(), FlagFM);

    // ---- NS collision and streaming ----
    NSLattice.template ApplyCellDynamics<NSTaskSelector>(FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // ---- PF collision and streaming ----
#ifndef DIAGNOSTIC_FREEZE_PHI
    PFLattice.template ApplyCellDynamics<PFCollisionTaskSelector>(FlagFM);
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();
#endif

    ++MainLoopTimer;
    ++OutputTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      PFLattice.getField<GRAD<T, 2>>().Communicate();
      PFLattice.getField<PHI<T>>().Communicate();
      NSLattice.getField<VELOCITY<T, 2>>().Communicate();
      NSLattice.getField<RHO<T>>().Communicate();
      MFLattice.getField<ff::PSI<T>>().Communicate();
      MFLattice.getField<ff::H_FIELD<T, 2>>().Communicate();
      MFLattice.getField<ff::H_SQ<T>>().Communicate();
      MFLattice.getField<ff::F_MAGNETIC<T, 2>>().Communicate();

      OutputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MainWriter.WriteBinary(MainLoopTimer());

#ifdef DIAGNOSTIC_FREEZE_PHI
      // Write diagnostic text output: centerline cut of u and F_m
      {
        auto& velField = NSLattice.getField<VELOCITY<T, 2>>();
        auto& forceField = NSLattice.getField<FORCE<T, 2>>();
        auto& phiField = PFLattice.getField<PHI<T>>();
        auto& chiField = MFLattice.getField<ff::CHI<T>>();

        // Find which block contains the bubble center y-line
        T y_target = Bubble_Center[1] * Cell_Len;

        std::ostringstream fname;
        fname << "diagnose_force_step" << std::setw(6) << std::setfill('0')
              << MainLoopTimer() << ".csv";
        std::ofstream diagFile(fname.str());
        diagFile << "# Diagnostic: frozen phi, u and F_m along y=" << y_target << "\n";
        diagFile << "# step,x,y,ux,uy,Fmx,Fmy,phi,chi\n";

        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          T voxelSize = block.getVoxelSize();
          const auto& blockVel = velField.getBlockField(blockid);
          const auto& blockForce = forceField.getBlockField(blockid);
          const auto& blockPhi = phiField.getBlockField(blockid);
          const auto& blockChi = chiField.getBlockField(blockid);
          int overlap = block.getOverlap();
          T y_min = block.getMin()[1];

          // Find j index closest to y_target
          int j_target = -1;
          T best_dist = T(999);
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            T y = y_min + static_cast<T>(j) * voxelSize;
            T dist = std::abs(y - y_target);
            if (dist < best_dist) {
              best_dist = dist;
              j_target = j;
            }
          }

          if (j_target >= overlap && j_target < block.getNy() - overlap) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = j_target * proj[1] + i;
              T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
              T y = y_min + static_cast<T>(j_target) * voxelSize;
              const auto& vel = blockVel.get(id);
              const auto& force = blockForce.get(id);
              T phi = blockPhi.get(id);
              T chi = blockChi.get(id);
              diagFile << MainLoopTimer() << ","
                       << x << "," << y << ","
                       << vel[0] << "," << vel[1] << ","
                       << force[0] << "," << force[1] << ","
                       << phi << "," << chi << "\n";
            }
          }
        }
        diagFile.close();
        IF_MPI_RANK(0) {
          std::cout << "[DIAGNOSTIC] Wrote " << fname.str() << "\n";
        }
      }
#endif

#ifdef DIAGNOSTIC_H_LINE
      // Write H-field centerline: Hx, Hy along y = bubble_center_y
      {
        auto& hField = MFLattice.getField<ff::H_FIELD<T, 2>>();
        auto& phiField = PFLattice.getField<PHI<T>>();
        auto& psiField = MFLattice.getField<ff::PSI<T>>();

        T y_target = Bubble_Center[1] * Cell_Len;

        std::ostringstream fnameH;
        fnameH << "diagnose_H_line_step" << std::setw(6) << std::setfill('0')
               << MainLoopTimer() << ".csv";
        std::ofstream hFile(fnameH.str());
        hFile << "# Diagnostic: H = -grad(psi) along y=" << y_target << "\n";
        hFile << "# step,x,y,Hx,Hy,psi,phi\n";

        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          T voxelSize = block.getVoxelSize();
          const auto& blockH = hField.getBlockField(blockid);
          const auto& blockPsi = psiField.getBlockField(blockid);
          const auto& blockPhi = phiField.getBlockField(blockid);
          int overlap = block.getOverlap();
          T y_min = block.getMin()[1];

          int j_target = -1;
          T best_dist = T(999);
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            T y = y_min + static_cast<T>(j) * voxelSize;
            T dist = std::abs(y - y_target);
            if (dist < best_dist) {
              best_dist = dist;
              j_target = j;
            }
          }

          if (j_target >= overlap && j_target < block.getNy() - overlap) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = j_target * proj[1] + i;
              T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
              T y = y_min + static_cast<T>(j_target) * voxelSize;
              const auto& H = blockH.get(id);
              T psi = blockPsi.get(id);
              T phi = blockPhi.get(id);
              hFile << MainLoopTimer() << ","
                    << x << "," << y << ","
                    << H[0] << "," << H[1] << ","
                    << psi << "," << phi << "\n";
            }
          }
        }
        hFile.close();
        IF_MPI_RANK(0) {
          std::cout << "[DIAGNOSTIC-H] Wrote " << fnameH.str() << "\n";
        }
      }
#endif

    }
  }

#ifdef DIAGNOSTIC_FREEZE_PHI
  // Print summary of spurious forces in solvent (chi≈0) region
  {
    auto& forceField = NSLattice.getField<FORCE<T, 2>>();
    auto& velField = NSLattice.getField<VELOCITY<T, 2>>();
    auto& phiField = PFLattice.getField<PHI<T>>();

    T max_F_solvent = T{0}, max_u_solvent = T{0};
    T max_F_all = T{0}, max_u_all = T{0};
    int solvent_count = 0;
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      const auto& blockForce = forceField.getBlockField(blockid);
      const auto& blockVel = velField.getBlockField(blockid);
      const auto& blockPhi = phiField.getBlockField(blockid);
      int overlap = block.getOverlap();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j * proj[1] + i;
          T phi = blockPhi.get(id);
          const auto& F = blockForce.get(id);
          const auto& u = blockVel.get(id);
          T Fmag = std::sqrt(F[0]*F[0] + F[1]*F[1]);
          T umag = std::sqrt(u[0]*u[0] + u[1]*u[1]);
          if (Fmag > max_F_all) max_F_all = Fmag;
          if (umag > max_u_all) max_u_all = umag;
          // Solvent = non-magnetic region (phi < 0.05, chi ≈ 0)
          if (phi < T{0.05}) {
            ++solvent_count;
            if (Fmag > max_F_solvent) max_F_solvent = Fmag;
            if (umag > max_u_solvent) max_u_solvent = umag;
          }
        }
      }
    }
    IF_MPI_RANK(0) {
      std::cout << "\n[DIAGNOSTIC SUMMARY] Final state (phi frozen, MF active):\n";
      std::cout << "  Solvent region (phi<0.05, " << solvent_count << " cells):\n";
      std::cout << "    max |F_m| = " << max_F_solvent << "\n";
      std::cout << "    max |u|   = " << max_u_solvent << "\n";
      std::cout << "  Full domain:\n";
      std::cout << "    max |F_m| = " << max_F_all << "\n";
      std::cout << "    max |u|   = " << max_u_all << "\n";
      if (max_F_solvent > T{1e-10}) {
        std::cout << "  => SPURIOUS FORCES DETECTED in solvent region!\n";
        std::cout << "     Suspect 1 (MST divergence error) is active.\n";
      } else {
        std::cout << "  => No spurious forces in solvent region.\n";
      }
      std::cout << std::endl;
    }
  }
#endif

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  MainLoopTimer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  Printer::Print("Total PhysTime", BaseConv.getPhysTime(MainLoopTimer()));
  Printer::Endl();
  return 0;
}
