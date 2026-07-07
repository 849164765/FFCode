// ferroBubbleRising2d.cpp
// 2D bubble rising in ferrofluid — NS + PF + Magnetic Field coupling
// Guo et al. (2025) Phys. Fluids 37, 022148 — Case C
// Three coupled lattices: NS(D2Q9) + PF(D2Q9) + MF(D2Q5)
// periodic X, no-slip Y, uniform vertical H₀

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

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

// bubble
T Bubble_Radius;
T Bubble_Diam;
Vector<T, LatSet::d> Bubble_Center;

// phase field
T Interface_Width;
T Mobility;
T Tau_phi;
T Omega_phi;
T Kappa;
T Beta;

// two-phase
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

// magnetic field — Guo et al. (2025), Hu & Li (2018)
T Mu_l;         // μ_l  — light fluid permeability
T Mu_h;         // μ_h  — heavy fluid permeability
T Chi_l;        // χ_l  — light fluid susceptibility
T Chi_h;        // χ_h  — heavy fluid susceptibility
T H0;           // H₀   — applied field strength (vertical)
T Epsilon;       // ε   — free parameter for magnetic pseudo-time
T Tau_mag_ref;   // τ_mag reference (Eq.42-43: default 4.0)
T Bom;          // magnetic Bond number: μ₀χ₀H₀²D/(2σ)

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("ferroBubbleRising2d.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Bubble_Diam = param_reader.getValue<T>("Bubble", "Diameter");
  Bubble_Radius = Bubble_Diam / T{2};
  Bubble_Center[0] = param_reader.getValue<T>("Bubble", "CenterX");
  Bubble_Center[1] = param_reader.getValue<T>("Bubble", "CenterY");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");
  Eo = param_reader.getValue<T>("Two_Phase", "Eo");
  Re = param_reader.getValue<T>("Two_Phase", "Re");

  // Magnetic parameters
  Mu_l = param_reader.getValue<T>("Magnetic", "Mu_l");
  Mu_h = param_reader.getValue<T>("Magnetic", "Mu_h");
  Chi_l = param_reader.getValue<T>("Magnetic", "Chi_l");
  Chi_h = param_reader.getValue<T>("Magnetic", "Chi_h");
  H0 = param_reader.getValue<T>("Magnetic", "H0");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ---- Guo et al. (2025) Case C hardcoded lattice parameters ----
  // ρ_l=1, ρ_h=1000, η_h/η_l=100, Re=40, Eo=20
  // τ_ns = 0.6 → ν_h = (τ-0.5)/3 = 0.03333, η_h = ν_h·ρ_h = 33.33
  // Re = sqrt(g·ρ_h²·D³)/η_h  →  g = (Re·η_h)²/(ρ_h²·D³)
  // Eo = g·ρ_h·D²/σ           →  σ = g·ρ_h·D²/Eo
  T D = Bubble_Diam;
  T tau_ns_target = T{0.6};
  T nu_h = (tau_ns_target - T{0.5}) * LatSet::cs2;
  eta_h = nu_h * rho_h;
  eta_l = eta_h / T{100};

  T g = (Re * eta_h) * (Re * eta_h) / (rho_h * rho_h * D * D * D);
  gravity = g;

  sigma = gravity * rho_h * D * D / Eo;
  U_g = std::sqrt(gravity * D);

  Beta = T{12} * sigma / Interface_Width;
  Kappa = T{3} * Interface_Width * sigma * T{0.5};
  Tau_phi = T{3} * Mobility + T{0.5};
  Omega_phi = T{1} / Tau_phi;
  Tau_ns = tau_ns_target;

  // Magnetic diffusion: τ(x) = 0.5 + ε·μ(x)/c_s²  (Eq.42-43)
  Tau_mag_ref = T{4.0};
  T mu_avg = (Mu_l + Mu_h) / T{2};
  Epsilon = (Tau_mag_ref - T{0.5}) * MFLatSet::cs2 / mu_avg;

  // Magnetic Bond number: Bom = μ₀H₀²D/(2σ), μ₀=1  (Eq.65)
  Bom = H0 * H0 * D / (T{2} * sigma);

  T cs = std::sqrt(LatSet::cs2);
  T Ma = U_g / cs;

  MPI_RANK(0) {
    std::cout << "-----Ferrofluid Bubble Rising (Guo2025 Case C)-----\\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj
              << "  BlockCellLen=" << BlockCellLen << "\\n";
    std::cout << "[Bubble]: D=" << Bubble_Diam
              << "  Center=(" << Bubble_Center[0] << ","
              << Bubble_Center[1] << ")\\n";
    std::cout << "[Fluid]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " eta_l=" << eta_l << " eta_h=" << eta_h
              << "  Re=" << Re << " Eo=" << Eo << "\\n";
    std::cout << "[Derived]: g=" << gravity << " sigma=" << sigma
              << " U_g=" << U_g << " Ma=" << Ma << "\\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " tau_phi=" << Tau_phi
              << " beta=" << Beta << " kappa=" << Kappa << "\\n";
    std::cout << "[Magnetic]: mu_l=" << Mu_l << " mu_h=" << Mu_h
              << " chi_l=" << Chi_l << " chi_h=" << Chi_h
              << " H0=" << H0 << " Bom=" << Bom
              << " eps=" << Epsilon << " tau_ref=" << Tau_mag_ref << "\\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep
              << "  OutputStep=" << OutputStep << "\\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\\n";
#endif
    std::cout << "-----------------------------------------------------\\n";
  }

  if (Ma > T(0.2)) {
    MPI_RANK(0) {
      std::cerr << "[Warning] Ma=" << Ma << " > 0.2\\n";
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

  Printer::Print_BigBanner(std::string("Initializing Bubble Rising (Guo2025 Case C)..."));

  readParam();

  // ------------------ converters ------------------
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_ns);

  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_phi);

  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni, T(0.01), Tau_mag_ref);

  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // ------------------ geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>(T(0), T(0)),
                    Vector<T, 2>(T(Ni * Cell_Len), T(Nj * Cell_Len)));
  AABB<T, 2> left(Vector<T, 2>(T(-Cell_Len), T(0)),
                  Vector<T, 2>(T(0), T(Nj * Cell_Len)));
  AABB<T, 2> right(Vector<T, 2>(T(Ni * Cell_Len), T(0)),
                   Vector<T, 2>(T((Ni + 1) * Cell_Len), T(Nj * Cell_Len)));

  BlockGeometryHelper2D<T> GeoHelper(Ni, Nj, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(4, 4);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.forEach(left, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(right, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  vtmo::ScalarWriter FlagWriter("flag", FlagFM);
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_BubbleRising", Geo, 1);
  GeoWriter.addWriterSet(FlagWriter);
  GeoWriter.WriteBinary();

  // ================== NS Lattice (D2Q9) ==================
  using NSFIELDS = TypePack<DENSITY<T>, VELOCITY<T, 2>, POP<T, LatSet::q>,
                            FORCE<T, LatSet::d>, OMEGA<T>, PRESSURE<T>>;
  T omega_ns_init = T{1} / Tau_ns;
  ValuePack NSInitValues(T{1}, Vector<T, 2>{T{0}, T{0}},
                         T{}, Vector<T, 2>{T{0}, T{0}}, omega_ns_init, T{});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSInitValues, BaseConv);

  // ================== PF Lattice (D2Q9) ==================
  using PFSELFFIELDS = TypePack<PHI<T>, POP<T, LatSet::q>, GRAD<T, LatSet::d>,
                                NORMAL<T, LatSet::d>, INTERFACEWIDTH<T>,
                                ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>,
                                ff::GRAVITY<T>, ff::BETA<T>, ff::KAPPA<T>,
                                ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>,
                                ff::DELTARHO<T>,
                                ff::MU_L<T>, ff::MU_H<T>,
                                ff::MU<T>>;
  using PFREFFIELDS = TypePack<VELOCITY<T, LatSet::d>>;
  using PFFIELDPACK = TypePack<PFSELFFIELDS, PFREFFIELDS>;
  ValuePack PFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}},
                         Vector<T, 2>{T{0}, T{0}}, Interface_Width,
                         T{}, T{},
                         gravity, Beta, Kappa, rho_l, rho_h, eta_l, eta_h,
                         rho_h - rho_l,
                         Mu_l, Mu_h,
                         T{});
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFFIELDPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFFIELDPACK> PFLattice(
    Geo, PFInitValues, PFBaseConv, &NSLattice.getField<VELOCITY<T, LatSet::d>>());

  ff::BroadcastAllParams<T>(PFLattice, rho_l, rho_h, eta_l, eta_h, gravity, Beta, Kappa);

  {
#ifdef MPI_ENABLED
    mpi().bCast(Mu_l, 0);
    mpi().bCast(Mu_h, 0);
#endif
    PFLattice.template getField<ff::MU_L<T>>().InitValue(Mu_l);
    PFLattice.template getField<ff::MU_H<T>>().InitValue(Mu_h);
  }

  // ================== MF Lattice (D2Q5) ==================
  using MFSELFFIELDS = TypePack<ff::PSI<T>, POP<T, MFLatSet::q>,
                                ff::H_FIELD<T, 2>, ff::H_SQ<T>,
                                ff::CHI<T>,
                                ff::CHI_L<T>, ff::CHI_H<T>,
                                ff::MU_L<T>, ff::MU_H<T>,
                                ff::MAGOMEGA<T>>;
  ValuePack MFInitValues(T{}, T{}, Vector<T, 2>{T{0}, T{0}}, T{},
                         T{}, Chi_l, Chi_h, Mu_l, Mu_h,
                         T{1} / Tau_mag_ref);
  using MFCELL = Cell<T, MFLatSet, MFSELFFIELDS>;
  BlockLatticeManager<T, MFLatSet, MFSELFFIELDS> MFLattice(Geo, MFInitValues, MFBaseConv);

  ff::BroadcastMagParams<T>(MFLattice, Mu_l, Mu_h, Chi_l, Chi_h);

  // ---- Init PHI: φ=1 outside (ferrofluid), φ=0 inside (bubble) ----
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

  // Init PF populations
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
  PFLattice.getField<ff::DELTARHO<T>>().InitValue(rho_h - rho_l);

  // Init MF populations (h_α = w_eq * ψ₀  with ψ₀ = -H₀*y)
  auto& psiField = MFLattice.getField<ff::PSI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPsi = psiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = 0;
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
        T psi = -H0 * y;
        blockPsi.get(id) = psi;
      }
    }
  }

  for (int block_idx = 0; block_idx < Geo.getBlockNum(); ++block_idx) {
    auto& blockLat = MFLattice.getBlockLat(block_idx);
    auto& blockPsiField = psiField.getBlockField(block_idx);
    const auto& block = Geo.getBlock(block_idx);
    int overlap = 0;
    const auto& proj = block.getProjection();
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        std::size_t id = j * proj[1] + i;
        MFCELL cell(id, blockLat);
        T psi = blockPsiField.get(id);
        constexpr T w_eq = T{1} / T{MFLatSet::q};
        for (unsigned int k = 0; k < MFLatSet::q; ++k) {
          cell[k] = w_eq * psi;
        }
      }
    }
  }

  // Init CHI from phi
  {
    auto& chiField = MFLattice.getField<ff::CHI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& blockPhi = phiField.getBlockField(blockid);
      auto& blockChi = chiField.getBlockField(blockid);
      int overlap = 0;
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j * proj[1] + i;
          T phi = blockPhi.get(id);
          blockChi.get(id) = Chi_l + phi * (Chi_h - Chi_l);
        }
      }
    }
    chiField.Communicate();
  }

  // Init OMEGA on MF lattice from chi (Eq.42-43)
  {
    auto& omegaField = MFLattice.getField<ff::MAGOMEGA<T>>();
    auto& chiField = MFLattice.getField<ff::CHI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& blockOmega = omegaField.getBlockField(blockid);
      auto& blockChi = chiField.getBlockField(blockid);
      int overlap = 0;
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = j * proj[1] + i;
          T chi = blockChi.get(id);
          T mu = T{1} + chi;
          T tau = T{0.5} + Epsilon * mu * MFLatSet::InvCs2;
          T omega = T{1} / tau;
          if (omega > T{1.9}) omega = T{1.9};
          if (omega < T{0.05}) omega = T{0.05};
          blockOmega.get(id) = omega;
        }
      }
    }
    omegaField.Communicate();
  }

  // Init NS populations to incompressible equilibrium
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
          T feq = latset::w<LatSet>(k) * (p_init
                  + LatSet::InvCs2 * uc
                  + uc * uc * T{0.5} * LatSet::InvCs4
                  - LatSet::InvCs2 * u2 * T{0.5});
          cell[k] = feq;
        }
      }
    }
  }

  // ================== Boundary Conditions ==================
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>,
                            BlockLatticeManager<T, LatSet, NSFIELDS>,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);

  using PFBLKLAT = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, PFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  using MFBLKLAT = BlockLatticeManager<T, MFLatSet, MFSELFFIELDS>;
  BBLikeFixedBlockBdManager<bounceback::normal<MFCELL>, MFBLKLAT,
                            BlockFieldManager<FLAG, T, LatSet::d>>
    MF_BB("MF_BB", MFLattice, FlagFM, BouncebackFlag, VoidFlag);

  // NS/PF/MF: periodic BC in X direction (MPI parallel decomposition)
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using LM_MF = BlockLatticeManager<T, MFLatSet, MFSELFFIELDS>;
  using FM = BlockFieldManager<FLAG, T, LatSet::d>;

  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_Periodic("NS_Periodic", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_Periodic("PF_Periodic", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PF_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

  FixedPeriodicBoundaryManager<LM_MF, FM>
      MF_Periodic("MF_Periodic", MFLattice, FlagFM, PeriodicFlag, VoidFlag);
  MF_Periodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  MF_Periodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);

#ifdef MPI_ENABLED
  NS_Periodic.SetupMPI(GeoHelper);
  PF_Periodic.SetupMPI(GeoHelper);
  MF_Periodic.SetupMPI(GeoHelper);
#endif

  // ================== Tasks ==================
  // --- NS tasks ---
  using NSBulkTask = tmp::Key_TypePair<BulkFlag,
    collision::MRTForce<NSCELL, FORCE<T, LatSet::d>>>;
  using NSPeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<NSCELL>>;
  using NSAllTasks = tmp::TupleWrapper<NSBulkTask, NSPeriodicTask>;
  using NSTaskSelector = tmp::TaskSelector<NSAllTasks, std::uint8_t, NSCELL>;

  // --- PF tasks ---
  using FFNormalTask = tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using FFLaplacianTask = tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using FFChemPotTask = tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using FFNormalSelector = TaskSelector<std::uint8_t, PFCELL, FFNormalTask>;
  using FFLaplacianSelector = TaskSelector<std::uint8_t, PFCELL, FFLaplacianTask>;
  using FFChemPotSelector = TaskSelector<std::uint8_t, PFCELL, FFChemPotTask>;

  using PFCollisionTask = tmp::Key_TypePair<
    BulkFlag,
    collision::MRTSource<equilibrium::FirstOrder<PFCELL>,
                         NORMAL<T, LatSet::d>, true, true>>;
  using PFPeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<PFCELL>>;
  using PFAllTasks = tmp::TupleWrapper<PFCollisionTask, PFPeriodicTask>;
  using PFCollisionTaskSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  // --- MF tasks ---
  using MFCollisionTask = tmp::Key_TypePair<BulkFlag, collision::MRTMag<MFCELL>>;
  using MFPeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<MFCELL>>;
  using MFBulkGradTask = tmp::Key_TypePair<BulkFlag, ff::MFGradient2D<MFCELL>>;
  using MFBulkHsqTask = tmp::Key_TypePair<BulkFlag, ff::MFHsq2D<MFCELL>>;
  using MFAllCollTasks = tmp::TupleWrapper<MFCollisionTask, MFPeriodicTask>;
  using MFAllGradTasks = tmp::TupleWrapper<MFBulkGradTask>;
  using MFAllHsqTasks = tmp::TupleWrapper<MFBulkHsqTask>;
  using MFCollisionSelector = tmp::TaskSelector<MFAllCollTasks, std::uint8_t, MFCELL>;
  using MFGradientSelector = tmp::TaskSelector<MFAllGradTasks, std::uint8_t, MFCELL>;
  using MFHsqSelector = tmp::TaskSelector<MFAllHsqTasks, std::uint8_t, MFCELL>;

  // --- PF → NS coupling tasks ---
  using STForceTask = tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using STForceTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STForceTask>;
  BlockLatManagerCoupling STCoupling(PFLattice, NSLattice);

  using GravForceTask = tmp::Key_TypePair<BulkFlag, ff::FFGravityForce2D<PFCELL, NSCELL>>;
  using GravForceTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravForceTask>;
  BlockLatManagerCoupling GravCoupling(PFLattice, NSLattice);

  using RhoOmegaTask = tmp::Key_TypePair<BulkFlag, ff::FFRhoOmegaUpdate2D<PFCELL, NSCELL>>;
  using RhoOmegaTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RhoOmegaTask>;
  BlockLatManagerCoupling RhoOmegaCoupling(PFLattice, NSLattice);

  using PreForceTask = tmp::Key_TypePair<BulkFlag, ff::FFPreForce2D<PFCELL, NSCELL>>;
  using PreForceTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PreForceTask>;
  BlockLatManagerCoupling PreForceCoupling(PFLattice, NSLattice);

  using ViscoForceTask = tmp::Key_TypePair<BulkFlag, ff::FFViscoForce2D<PFCELL, NSCELL>>;
  using ViscoForceTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, ViscoForceTask>;
  BlockLatManagerCoupling ViscoForceCoupling(PFLattice, NSLattice);

  using MagForceTask = tmp::Key_TypePair<BulkFlag, ff::MFForce2D<MFCELL, NSCELL>>;
  using MagForceTaskSelector = CoupledTaskSelector<std::uint8_t, MFCELL, NSCELL, MagForceTask>;
  BlockLatManagerCoupling MagForceCoupling(MFLattice, NSLattice);

  using ChiUpdateTask = tmp::Key_TypePair<BulkFlag, ff::FFCchiUpdate2D<PFCELL, MFCELL>>;
  using ChiUpdateTaskSelector = CoupledTaskSelector<std::uint8_t, PFCELL, MFCELL, ChiUpdateTask>;
  BlockLatManagerCoupling ChiUpdateCoupling(PFLattice, MFLattice);

  using MuUpdateTask = tmp::Key_TypePair<BulkFlag, ff::FFMuUpdate2D<PFCELL>>;
  using MuUpdateSelector = TaskSelector<std::uint8_t, PFCELL, MuUpdateTask>;

  // ================== Writers ==================
  vtmo::ScalarWriter PHIWriter("PHI", PFLattice.getField<PHI<T>>());
  vtmo::VectorWriter VecWriter("Velocity", NSLattice.getField<VELOCITY<T, 2>>());
  vtmo::VectorWriter ForceWriter("Force", NSLattice.getField<FORCE<T, 2>>());
  vtmo::ScalarWriter DensityWriter("Density", NSLattice.getField<DENSITY<T>>());
  vtmo::ScalarWriter PressureWriter("Pressure", NSLattice.getField<PRESSURE<T>>());
  vtmo::ScalarWriter PsiWriter("Psi", MFLattice.getField<ff::PSI<T>>());
  vtmo::VectorWriter HWriter("HField", MFLattice.getField<ff::H_FIELD<T, 2>>());
  vtmo::ScalarWriter HSqWriter("HSq", MFLattice.getField<ff::H_SQ<T>>());
  vtmo::ScalarWriter ChiWriter("Chi", MFLattice.getField<ff::CHI<T>>());

  vtmo::vtmWriter<T, 2> MainWriter("ferroBubbleRising2d", Geo);
  MainWriter.addWriterSet(PHIWriter, VecWriter, ForceWriter,
                          DensityWriter, PressureWriter,
                          PsiWriter, HWriter, HSqWriter, ChiWriter);

  // ================== Timer ==================
  Timer MainLoopTimer;
  Timer OutputTimer;

  // Initial communication
  PFLattice.NormalFullCommunicate();
  NSLattice.NormalFullCommunicate();
  MFLattice.NormalFullCommunicate();
  NS_Periodic.Apply();
  PF_Periodic.Apply();

  // Initialize ψ field = -H₀·y (nice initial guess for the SOR solver)
  {
    auto& psiField = MFLattice.getField<ff::PSI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& bp = psiField.getBlockField(blockid);
      T minY = block.getMin()[1];
      T dy = block.getVoxelSize();
      for (int j = 0; j < block.getNy(); ++j) {
        T psi_bg = -H0 * (minY + j * dy);
        for (int i = 0; i < block.getNx(); ++i)
          bp.get(j * proj[1] + i) = psi_bg;
      }
    }
    MFLattice.NormalFullCommunicate();
  }

  // Initial PF gradients
  PFLattice.template ApplyInnerCellDynamics<FFNormalSelector>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<FFLaplacianSelector>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<FFChemPotSelector>(FlagFM);
  PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
  PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);

  // Initial NS rho & omega
  RhoOmegaCoupling.ApplyInnerCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);

  // Initial MF gradient and |H|²
  MFLattice.template ApplyInnerCellDynamics<MFGradientSelector>(FlagFM);
  MFLattice.template ApplyInnerCellDynamics<MFHsqSelector>(FlagFM);
  ff::CommunicateMagFields<T>(MFLattice);

  MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  // ================== Main Loop ==================
  while (MainLoopTimer() < MaxStep) {

    // ---- Phase A: Force accumulation ----
    RhoOmegaCoupling.ApplyInnerCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    STCoupling.ApplyInnerCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);
    GravCoupling.ApplyInnerCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);
    PreForceCoupling.ApplyInnerCellDynamics<PreForceTaskSelector>(MainLoopTimer(), FlagFM);

    // ---- Phase A_mag: Magnetic force ----
    // CHI update
    {
      auto& chiField = MFLattice.getField<ff::CHI<T>>();
      auto& phiField = PFLattice.getField<PHI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPhi = phiField.getBlockField(blockid);
        auto& blockChi = chiField.getBlockField(blockid);
        int overlap = 0;
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            T phi = blockPhi.get(id);
            blockChi.get(id) = Chi_l + phi * (Chi_h - Chi_l);
          }
        }
      }
      chiField.Communicate();
    }

    // MAGOMEGA update (required for MFLattice NormalFullCommunicate consistency)
    {
      auto& omegaField = MFLattice.getField<ff::MAGOMEGA<T>>();
      auto& chiField = MFLattice.getField<ff::CHI<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& bo = omegaField.getBlockField(blockid);
        auto& bc = chiField.getBlockField(blockid);
        for (int j = 0; j < block.getNy(); ++j)
          for (int i = 0; i < block.getNx(); ++i) {
            std::size_t id = j * proj[1] + i;
            T mu = T{1} + bc.get(id);
            T tau = T{0.5} + Epsilon * mu * MFLatSet::InvCs2;
            T omega = T{1} / tau;
            if (omega > T{1.9}) omega = T{1.9};
            if (omega < T{0.05}) omega = T{0.05};
            bo.get(id) = omega;
          }
      }
      omegaField.Communicate();
    }

    // ---- MF LBM: Fortran-style collision + ghost POP BC ----
    // Collision: MRTMag on per-cell omega from mu
    // Ghost BC: bounceback-like with flux 2*wt*H₀*μ/cs²
    // After streaming: ψ=Σpop, ghost ψ extrapolation
    // One step per NS step (no sub-iterations) — steady state reached
    // within ~O(R²/D) ≈ 100-200 NS steps
    {
      // MF collision
      MFLattice.template ApplyInnerCellDynamics<MFCollisionSelector>(FlagFM);
      MF_Periodic.Apply();

      // Ghost POP BC: Fortran-style (collisionPsiDF2D lines 472-488)
      // ghost[incoming] = wall[opposite] ∓ 2*wt*H₀*μ/cs²
      {
        auto& chiField  = MFLattice.getField<ff::CHI<T>>();
        T invCs2 = MFLatSet::InvCs2;
        T H_global = T(Nj) * Cell_Len;
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& blockLat = MFLattice.getBlockLat(blockid);
          auto& bc = chiField.getBlockField(blockid);
          int nx = block.getNx(), ny = block.getNy();
          int ov = block.getOverlap();
          T minY = block.getMin()[1], maxY = block.getMax()[1];
          T dy = block.getVoxelSize();
          if (minY < Cell_Len * T{1.5}) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_g = 0 * proj[1] + i;
              std::size_t id_w = ov * proj[1] + i;
              MFCELL gc(id_g, blockLat);
              MFCELL wc(id_w, blockLat);
              T mu = T{1} + bc.get(id_w);
              for (unsigned int k = 1; k < MFLatSet::q; ++k) {
                int cky = latset::c<MFLatSet>(k)[1];
                if (cky > 0) {  // direction pointing INTO domain from bottom
                  int opp = latset::opp<MFLatSet>(k);
                  T flux = T{2} * latset::w<MFLatSet>(k) * H0 * mu * invCs2;
                  gc[k] = wc[opp] - flux * dy;
                }
              }
            }
          }
          if (maxY > H_global - Cell_Len * T{1.5}) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_g = (ny - 1) * proj[1] + i;
              std::size_t id_w = (ny - 1 - ov) * proj[1] + i;
              MFCELL gc(id_g, blockLat);
              MFCELL wc(id_w, blockLat);
              T mu = T{1} + bc.get(id_w);
              for (unsigned int k = 1; k < MFLatSet::q; ++k) {
                int cky = latset::c<MFLatSet>(k)[1];
                if (cky < 0) {  // direction pointing INTO domain from top
                  int opp = latset::opp<MFLatSet>(k);
                  T flux = T{2} * latset::w<MFLatSet>(k) * H0 * mu * invCs2;
                  gc[k] = wc[opp] + flux * dy;
                }
              }
            }
          }
        }
      }

      // Bounceback + streaming
      MF_BB.Apply(MainLoopTimer());
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

      // Psi = Σ pop
      {
        auto& psiField = MFLattice.getField<ff::PSI<T>>();
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& bp = psiField.getBlockField(blockid);
          auto& blockLat = MFLattice.getBlockLat(blockid);
          int ov = block.getOverlap();
          for (int j = ov; j < block.getNy() - ov; ++j)
            for (int i = ov; i < block.getNx() - ov; ++i) {
              MFCELL cell(j * proj[1] + i, blockLat);
              T psi = T{};
              for (unsigned int k = 0; k < MFLatSet::q; ++k) psi += cell[k];
              bp.get(j * proj[1] + i) = psi;
            }
        }
        psiField.Communicate();
      }

      // Psi ghost extrapolation (2nd order, like Fortran computeOrderMacro2D)
      {
        auto& psiField = MFLattice.getField<ff::PSI<T>>();
        T H_global = T(Nj) * Cell_Len;
        T dy = Cell_Len;  // voxel size
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& bp = psiField.getBlockField(blockid);
          int nx = block.getNx(), ny = block.getNy();
          int ov = block.getOverlap();
          T minY = block.getMin()[1], maxY = block.getMax()[1];
          if (minY < Cell_Len * T{1.5}) {
            // Only 1 ghost layer supported; use 2nd-order Neumann extrapolation
            for (int i = 0; i < nx; ++i) {
              std::size_t ig = 0 * proj[1] + i;
              std::size_t i1 = ov * proj[1] + i;
              std::size_t i2 = (ov + 1) * proj[1] + i;
              bp.get(ig) = (T{4} * bp.get(i1) - bp.get(i2) - T{2} * H0 * dy) / T{3};
            }
          }
          if (maxY > H_global - Cell_Len * T{1.5}) {
            for (int i = 0; i < nx; ++i) {
              std::size_t ig = (ny - 1) * proj[1] + i;
              std::size_t i1 = (ny - 1 - ov) * proj[1] + i;
              std::size_t i2 = (ny - 2 - ov) * proj[1] + i;
              bp.get(ig) = (T{4} * bp.get(i1) - bp.get(i2) + T{2} * H0 * dy) / T{3};
            }
          }
        }
        psiField.Communicate();
      }
    }

    // MF: H = -∇ψ, |H|²
    MFLattice.template ApplyInnerCellDynamics<MFGradientSelector>(FlagFM);
    MFLattice.template ApplyInnerCellDynamics<MFHsqSelector>(FlagFM);
    ff::CommunicateMagFields<T>(MFLattice);

    // MF → NS: F_m = (χ/2)·∇(|H|²)
    MagForceCoupling.ApplyInnerCellDynamics<MagForceTaskSelector>(MainLoopTimer(), FlagFM);

    // Communicate FORCE to ghosts
    NSLattice.getField<FORCE<T, LatSet::d>>().Communicate();

    // ---- Phase B: PF collision ----
    PFLattice.template ApplyInnerCellDynamics<PFCollisionTaskSelector>(FlagFM);
    PF_Periodic.Apply();
    PFLattice.NormalFullCommunicate();

    // ---- Phase C: NS collision ----
    ViscoForceCoupling.ApplyInnerCellDynamics<ViscoForceTaskSelector>(MainLoopTimer(), FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSTaskSelector>(FlagFM);
    NS_Periodic.Apply();
    NSLattice.NormalFullCommunicate();

    // ---- Phase D: Streaming ----
    // PF stream
    PF_BB.Apply(MainLoopTimer());
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

    // NS stream
    NS_BB.Apply(MainLoopTimer());
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

    // ---- Phase E: Macro update ----
    // PF phi from pops
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

    // φ=1 at walls + ghost cell fix
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
          for (int i = overlap; i < nx - overlap; ++i) {
            blockPhi.get(overlap * proj[1] + i) = T{1};
          }
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{1};
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            blockPhi.get((ny - 1 - overlap) * proj[1] + i) = T{1};
          }
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{1};
            }
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // PF gradients, normal, laplacian, chempot
    PFLattice.template ApplyInnerCellDynamics<FFNormalSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFLaplacianSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFChemPotSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // Wall grad_phi and chempot extrapolation
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
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id_bot = overlap * proj[1] + i;
            std::size_t id_bot1 = (overlap + 1) * proj[1] + i;
            blockGrad.get(id_bot)[1] = blockGrad.get(id_bot1)[1];
          }
        }
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
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            std::size_t id1 = overlap * proj[1] + i;
            std::size_t id2 = (overlap + 1) * proj[1] + i;
            std::size_t id3 = (overlap + 2) * proj[1] + i;
            blockChpoten.get(id1) = (T{4} * blockChpoten.get(id2) - blockChpoten.get(id3)) / T{3};
          }
        }
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

    // Ghost extrapolation
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

    // NS macro from streamed pops
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
            T pres = T{0}, ux_raw = T{0}, uy_raw = T{0};
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              pres += cell[k];
              ux_raw += latset::c<LatSet>(k)[0] * cell[k];
              uy_raw += latset::c<LatSet>(k)[1] * cell[k];
            }
            T rho = blockRho.get(id);
            const auto& F = blockForce.get(id);
            blockPres.get(id) = pres;
            blockVel.get(id) = Vector<T, 2>{ux_raw + T{0.5} * F[0] / rho,
                                            uy_raw + T{0.5} * F[1] / rho};
          }
        }
      }
    }
    NSLattice.getField<VELOCITY<T, LatSet::d>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    // NS ghost extrapolation
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
      MFLattice.getField<ff::PSI<T>>().Communicate();
      MFLattice.getField<ff::H_FIELD<T, 2>>().Communicate();
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
