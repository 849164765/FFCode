// rosensweig2d.cpp - clean Rosensweig instability benchmark
// Phase field + incompressible flow + scalar-potential magnetic field.

#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

#include <cmath>
#include <cstdlib>
#include <type_traits>
#include <utility>
#include <vector>

using T = FLOAT;
using LatSet = D2Q9<T>;
using MFLatSet = D2Q5<T>;
using namespace mfield;

int Thread_Num = 1;
T RosenPressureScale = T{0.6};

struct RosenParams {
  int Ni{}, Nj{}, BlockCellLen{}, ThreadNum{};
  T CellLen{};
  T InterfaceY{}, SeedAmplitude{}, SeedPeriods{}, InterfaceWidth{}, Mobility{};
  T rho_l{}, rho_h{}, eta_l{}, eta_h{}, sigma{}, DeltaRho{}, gravity{};
  T chi_l{}, chi_h{}, mu_l{}, mu_h{}, H0{}, H0_kAm{}, Hc_kAm{};
  T HcLat{}, lambdaCLat{};
  T LxMm{}, LyMm{}, sigmaMm{}, lambdaCMm{};
  T PhiFloor{}, PhiCeiling{};
  int PsiSolverIter{};
  T PsiSolverK{};
  int PsiWallLayers{}, MagneticForceWallBand{};
  int MaxStep{}, OutputStep{};

  static constexpr T Pi = T{3.1415926535897932384626433832795};

  void read(int argc, char** argv) {
    std::string iniName = "rosensweig2d.ini";
    T h0Override = T{-1};
    if (argc > 1) iniName = argv[1];
    if (argc > 2) h0Override = std::atof(argv[2]);

    iniReader r(iniName);
    ThreadNum = r.getValue<int>("parallel", "thread_num");
    Ni = r.getValue<int>("Mesh", "Ni");
    Nj = r.getValue<int>("Mesh", "Nj");
    CellLen = r.getValue<T>("Mesh", "Cell_Len");
    BlockCellLen = r.getValue<int>("Mesh", "BlockCellLen");

    InterfaceY = r.getValue<T>("Interface", "Y0_cells");
    SeedAmplitude = r.getValue<T>("Interface", "SeedAmplitude_cells");
    SeedPeriods = r.getValue<T>("Interface", "SeedPeriods");

    InterfaceWidth = r.getValue<T>("Phase_Field", "InterfaceWidth");
    Mobility = r.getValue<T>("Phase_Field", "Mobility");

    rho_l = r.getValue<T>("Two_Phase", "rho_l");
    rho_h = r.getValue<T>("Two_Phase", "rho_h");
    eta_l = r.getValue<T>("Two_Phase", "eta_l");
    eta_h = r.getValue<T>("Two_Phase", "eta_h");
    sigma = r.getValue<T>("Two_Phase", "sigma");

    LxMm = r.getValue<T>("Physical", "Lx_mm");
    LyMm = r.getValue<T>("Physical", "Ly_mm");
    sigmaMm = r.getValue<T>("Physical", "sigma_mNm");
    lambdaCMm = r.getValue<T>("Physical", "lambda_c_mm");

    chi_l = r.getValue<T>("Magnetic_Field", "chi_l");
    chi_h = r.getValue<T>("Magnetic_Field", "chi_h");
    mu_l = r.getValue<T>("Magnetic_Field", "mu_l");
    mu_h = r.getValue<T>("Magnetic_Field", "mu_h");
    H0_kAm = r.getValue<T>("Magnetic_Field", "H0_kAm");
    Hc_kAm = r.getValue<T>("Magnetic_Field", "Hc_kAm");
    PsiSolverIter = r.getValue<int>("Magnetic_Field", "PsiSolverIter");
    PsiSolverK = r.getValue<T>("Magnetic_Field", "PsiSolverK");

    RosenPressureScale = r.getValue<T>("Numerics", "PressureCorrectionScale");
    PhiFloor = r.getValue<T>("Numerics", "PhiFloor");
    PhiCeiling = r.getValue<T>("Numerics", "PhiCeiling");
    PsiWallLayers = r.getValue<int>("Numerics", "PsiWallLayers");
    MagneticForceWallBand = r.getValue<int>("Numerics", "MagneticForceWallBand");

    if (h0Override >= T{0}) H0_kAm = h0Override;
    MaxStep = r.getValue<int>("Simulation_Settings", "TotalStep");
    if (argc > 3) MaxStep = std::atoi(argv[3]);
    OutputStep = r.getValue<int>("Simulation_Settings", "OutputStep");

    DeltaRho = rho_h - rho_l;
    lambdaCLat = lambdaCMm * T(Ni) / LxMm;
    gravity = sigma /
      (DeltaRho * std::pow(lambdaCLat / (T{2} * Pi), T{2}));

    const T mu0 = T{4} * Pi * T{1e-7};
    const T sigmaPhysical = sigmaMm * T{1e-3};
    const T lambdaPhysical = lambdaCMm * T{1e-3};
    const T HcPhysical = Hc_kAm * T{1e3};
    HcLat = std::sqrt(mu0 * HcPhysical * HcPhysical * lambdaPhysical * sigma /
                      (sigmaPhysical * lambdaCLat));
    H0 = H0_kAm * HcLat / Hc_kAm;

    MPI_RANK(0) {
      printf("---- Clean Rosensweig benchmark ----\n");
      printf("mesh=%dx%d blocks=%d cell=%.3f\n", Ni, Nj, BlockCellLen, CellLen);
      printf("interface_y=%.3f seed_amp=%.3f seed_periods=%.0f W=%.3f\n",
             InterfaceY, SeedAmplitude, SeedPeriods, InterfaceWidth);
      printf("rho=(%.5f,%.5f) eta=(%.5f,%.5f) sigma=%.6g gravity=%.6g\n",
             rho_l, rho_h, eta_l, eta_h, sigma, gravity);
      printf("mu=(%.5f,%.5f) chi=(%.5f,%.5f) H0=%.6g HcLat=%.6g\n",
             mu_l, mu_h, chi_l, chi_h, H0, HcLat);
      printf("H0 physical=%.4f kA/m, Hc=%.4f kA/m, lambda_c=%.4f cells\n",
             H0_kAm, Hc_kAm, lambdaCLat);
      printf("psi_iter=%d psi_K=%.6g pressure_scale=%.6g\n",
             PsiSolverIter, PsiSolverK, RosenPressureScale);
      printf("phi_floor=%.6g phi_ceiling=%.6g force_wall_band=%d\n",
             PhiFloor, PhiCeiling, MagneticForceWallBand);
      printf("steps=%d output=%d\n", MaxStep, OutputStep);
      printf("--------------------------------------\n");
    }
  }
};

template <typename PFCELL, typename NSCELL>
struct RosenPressureForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;

  __any__ static void apply(PFCELL& pfCell, NSCELL& nsCell) {
    const T pressure = nsCell.template get<PRESSURE<T>>();
    const T deltaRho = pfCell.template get<ff::DELTARHO<T>>();
    const auto& gradPhi = pfCell.template get<GRAD<T, LatSet::d>>();
    const T coefficient = -pressure * deltaRho * RosenPressureScale / T{3};
    auto& force = nsCell.template get<FORCE<T, LatSet::d>>();
    force[0] += coefficient * gradPhi[0];
    force[1] += coefficient * gradPhi[1];
  }
};

int main(int argc, char** argv) {
  constexpr std::uint8_t VoidFlag = 1;
  constexpr std::uint8_t BulkFlag = 2;
  constexpr std::uint8_t BouncebackFlag = 4;
  constexpr std::uint8_t PeriodicFlag = 8;

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  RosenParams p;
  p.read(argc, argv);
  Thread_Num = p.ThreadNum;

  BaseConverter<T> NSConv(LatSet::cs2);
  const T tauNS = T{0.5} + p.eta_h / p.rho_h / LatSet::cs2;
  NSConv.SimplifiedConverterFromRT(p.Ni, T{0.01}, tauNS);
  BaseConverter<T> PFConv(LatSet::cs2);
  const T tauPhi = T{0.5} + T{3} * p.Mobility;
  PFConv.SimplifiedConverterFromRT(p.Ni, T{0.01}, tauPhi);
  BaseConverter<T> MFConv(MFLatSet::cs2);
  MFConv.SimplifiedConverterFromRT(p.Ni, T{0.01}, T{1});

  AABB<T, 2> domain({0, 0}, {T(p.Ni) * p.CellLen, T(p.Nj) * p.CellLen});
  AABB<T, 2> left({-p.CellLen, 0}, {0, T(p.Nj) * p.CellLen});
  AABB<T, 2> right({T(p.Ni) * p.CellLen, 0},
                   {T(p.Ni + 1) * p.CellLen, T(p.Nj) * p.CellLen});
  BlockGeometryHelper2D<T> geoHelper(
    p.Ni, p.Nj, domain, p.CellLen, p.BlockCellLen);
  geoHelper.CreateBlocks(8, 16);
  geoHelper.AdaptiveOptimization(mpi().getSize());
  geoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> geo(geoHelper);

  BlockFieldManager<FLAG, T, 2> flags(geo, VoidFlag);
  flags.forEach(domain, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, BulkFlag);
  });
  flags.forEach(left, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  flags.forEach(right, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  flags.template SetupBoundary<LatSet>(domain, BouncebackFlag);
  AABB<T, 2> leftColumn({0, p.CellLen},
                        {p.CellLen, T(p.Nj - 1) * p.CellLen});
  AABB<T, 2> rightColumn({T(p.Ni - 1) * p.CellLen, p.CellLen},
                         {T(p.Ni) * p.CellLen,
                          T(p.Nj - 1) * p.CellLen});
  flags.forEach(leftColumn, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, BulkFlag);
  });
  flags.forEach(rightColumn, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, BulkFlag);
  });

  using NSFIELDS = TypePack<
    DENSITY<T>, VELOCITY<T, 2>, POP<T, LatSet::q>, FORCE<T, 2>,
    OMEGA<T>, PRESSURE<T>>;
  ValuePack NSInit(
    T{1}, Vector<T, 2>{0, 0}, T{}, Vector<T, 2>{0, 0}, T{1} / tauNS, T{});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NS(geo, NSInit, NSConv);

  using PFFIELDS = TypePack<
    PHI<T>, POP<T, LatSet::q>, GRAD<T, 2>, NORMAL<T, 2>, INTERFACEWIDTH<T>,
    ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>, ff::GRAVITY<T>, ff::BETA<T>,
    ff::KAPPA<T>, ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>,
    ff::DELTARHO<T>>;
  using PFREF = TypePack<VELOCITY<T, 2>>;
  using PFPACK = TypePack<PFFIELDS, PFREF>;
  T beta = T{12} * p.sigma / p.InterfaceWidth;
  T kappa = T{3} * p.InterfaceWidth * p.sigma / T{2};
  ValuePack PFInit(
    T{}, T{}, Vector<T, 2>{0, 0}, Vector<T, 2>{0, 0}, p.InterfaceWidth,
    T{}, T{}, p.gravity, beta, kappa, p.rho_l, p.rho_h,
    p.eta_l, p.eta_h, p.DeltaRho);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFPACK> PF(
    geo, PFInit, PFConv, &NS.getField<VELOCITY<T, 2>>());
  ff::BroadcastAllParams<T>(
    PF, p.rho_l, p.rho_h, p.eta_l, p.eta_h, p.gravity, beta, kappa);
  PF.getField<ff::DELTARHO<T>>().InitValue(p.DeltaRho);

  using MFFIELDS = TypePack<
    PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>, HX<T>, HY<T>,
    HMAG<T>, POP<T, MFLatSet::q>, MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>,
    H_0<T>, PSI_K<T>>;
  using MFREF = TypePack<PHI<T>>;
  using MFPACK = TypePack<MFFIELDS, MFREF>;
  ValuePack MFInit(
    T{}, T{1}, T{p.mu_l}, T{p.chi_l}, T{}, T{}, T{}, T{},
    p.mu_l, p.mu_h, p.chi_l, p.chi_h, p.H0, p.PsiSolverK);
  using MFCELL = Cell<T, MFLatSet, ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T, MFLatSet, MFPACK> MF(
    geo, MFInit, MFConv, &PF.getField<PHI<T>>());
  BroadcastAllMFParams<T>(
    MF, p.mu_l, p.mu_h, p.chi_l, p.chi_h, p.H0, p.PsiSolverK);
  MF.getField<OMEGA_PSI<T>>().InitValue(T{1});

  const T y0 = p.InterfaceY * p.CellLen;
  const T width = p.InterfaceWidth * p.CellLen;
  const T domainHeight = T(p.Nj) * p.CellLen;
  const T domainWidth = T(p.Ni) * p.CellLen;
  const T waveNumber = T{2} * RosenParams::Pi * p.SeedPeriods / domainWidth;

  auto& phiField = PF.getField<PHI<T>>();
  for (int b = 0; b < geo.getBlockNum(); ++b) {
    const auto& block = geo.getBlock(b);
    const auto& projection = block.getProjection();
    auto& blockPhi = phiField.getBlockField(b);
    const T voxel = block.getVoxelSize();
    const T minX = block.getMin()[0];
    const T minY = block.getMin()[1];
    for (int j = 0; j < block.getNy(); ++j) {
      const T y = minY + T(j) * voxel;
      for (int i = 0; i < block.getNx(); ++i) {
        const T x = minX + T(i) * voxel;
        const T interfaceY = y0 + p.SeedAmplitude * p.CellLen *
          std::cos(waveNumber * x);
        blockPhi.get(j * projection[1] + i) =
          T{0.5} - T{0.5} * std::tanh(T{2} * (y - interfaceY) / width);
      }
    }
  }

  for (int b = 0; b < geo.getBlockNum(); ++b) {
    auto& blockLat = PF.getBlockLat(b);
    auto& blockPhi = phiField.getBlockField(b);
    const auto& block = geo.getBlock(b);
    const auto& projection = block.getProjection();
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        const std::size_t id = j * projection[1] + i;
        PFCELL cell(id, blockLat);
        const T phi = blockPhi.get(id);
        for (unsigned k = 0; k < LatSet::q; ++k) {
          cell[k] = latset::w<LatSet>(k) * phi;
        }
      }
    }
  }
  PF.getField<INTERFACEWIDTH<T>>().InitValue(p.InterfaceWidth);

  for (int b = 0; b < geo.getBlockNum(); ++b) {
    auto& blockLat = NS.getBlockLat(b);
    const auto& block = geo.getBlock(b);
    const auto& projection = block.getProjection();
    for (int j = 0; j < block.getNy(); ++j) {
      for (int i = 0; i < block.getNx(); ++i) {
        NSCELL cell(j * projection[1] + i, blockLat);
        for (unsigned k = 0; k < LatSet::q; ++k) {
          // 必须从零压初始化（与旧版 pz=0 一致）：w*(-3) 给出负密度/负压力，
          // 压力修正力 -p*DeltaRho/3*grad_phi 会被负 p 反向放大，界面整体
          // 上漂直至淹没顶壁，导致峰/谷提取全部为 0。
          cell[k] = T{0};
        }
      }
    }
  }

  // Seed the exact flat-interface solution plus the evanescent field of the
  // initial cosine perturbation. This removes a long magnetic transient while
  // keeping the boundary data tied to the declared initial interface.
  const T hBelow = p.H0 * p.mu_l / p.mu_h;
  const T psiAtY0 = -hBelow * y0;
  const T muRatio = p.mu_l / p.mu_h;
  const T ca = p.SeedAmplitude * p.CellLen * p.H0 * (T{1} - muRatio) *
    p.mu_h / (p.mu_h + p.mu_l);
  const T cb = -muRatio * ca;
  auto seedPsi = [&](T x, T y) {
    if (y >= y0) {
      return psiAtY0 - p.H0 * (y - y0) +
        ca * std::cos(waveNumber * x) * std::exp(-waveNumber * (y - y0));
    }
    return -hBelow * y +
      cb * std::cos(waveNumber * x) * std::exp(waveNumber * (y - y0));
  };
  auto& psiField = MF.getField<PSI<T>>();
  for (int b = 0; b < geo.getBlockNum(); ++b) {
    auto& blockLat = MF.getBlockLat(b);
    const auto& block = geo.getBlock(b);
    const auto& projection = block.getProjection();
    const T voxel = block.getVoxelSize();
    const T minX = block.getMin()[0];
    const T minY = block.getMin()[1];
    auto& blockPsi = psiField.getBlockField(b);
    for (int j = 0; j < block.getNy(); ++j) {
      const T y = minY + T(j) * voxel;
      for (int i = 0; i < block.getNx(); ++i) {
        const T x = minX + T(i) * voxel;
        const T psi = seedPsi(x, y);
        const std::size_t id = j * projection[1] + i;
        MFCELL cell(id, blockLat);
        for (unsigned k = 0; k < MFLatSet::q; ++k) {
          cell[k] = latset::w<MFLatSet>(k) * psi;
        }
        blockPsi.get(id) = psi;
      }
    }
  }

  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFPACK>;
  using LM_MF = BlockLatticeManager<T, MFLatSet, MFPACK>;
  using FM = BlockFieldManager<FLAG, T, 2>;
  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>, LM_NS, FM>
    NSBounce("NSBounce", NS, flags, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, LM_PF, FM>
    PFBounce("PFBounce", PF, flags, BouncebackFlag, VoidFlag);

  FixedPeriodicBoundaryManager<LM_NS, FM> NSPeriodic(
    "NSPeriodic", NS, flags, PeriodicFlag, VoidFlag);
  FixedPeriodicBoundaryManager<LM_PF, FM> PFPeriodic(
    "PFPeriodic", PF, flags, PeriodicFlag, VoidFlag);
  FixedPeriodicBoundaryManager<LM_MF, FM> MFPeriodic(
    "MFPeriodic", MF, flags, PeriodicFlag, VoidFlag);
  NSPeriodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NSPeriodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
  PFPeriodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PFPeriodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
  MFPeriodic.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  MFPeriodic.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
#ifdef MPI_ENABLED
  NSPeriodic.SetupMPI(geoHelper);
  PFPeriodic.SetupMPI(geoHelper);
  MFPeriodic.SetupMPI(geoHelper);
#endif

  using PFNormalTask = tmp::Key_TypePair<BulkFlag, ff::FF2D<PFCELL>>;
  using PFLaplacianTask = tmp::Key_TypePair<BulkFlag, ff::FFLaplacian2D<PFCELL>>;
  using PFChemicalTask = tmp::Key_TypePair<BulkFlag, ff::FFChemPotential2D<PFCELL>>;
  using PFNormalSelector = TaskSelector<std::uint8_t, PFCELL, PFNormalTask>;
  using PFLaplacianSelector = TaskSelector<std::uint8_t, PFCELL, PFLaplacianTask>;
  using PFChemicalSelector = TaskSelector<std::uint8_t, PFCELL, PFChemicalTask>;
  using PFCollisionTask = tmp::Key_TypePair<BulkFlag,
    collision::MRTSource<equilibrium::FirstOrder<PFCELL>, NORMAL<T, 2>, true, true>>;
  using PFPeriodicTask = tmp::Key_TypePair<PeriodicFlag,
    collision::PeriodicBoundary<PFCELL>>;
  using PFAllTasks = tmp::TupleWrapper<PFCollisionTask, PFPeriodicTask>;
  using PFSelector = tmp::TaskSelector<PFAllTasks, std::uint8_t, PFCELL>;

  using NSCollisionTask = tmp::Key_TypePair<BulkFlag,
    collision::MRTForce<NSCELL, FORCE<T, 2>>>;
  using NSPeriodicTask = tmp::Key_TypePair<PeriodicFlag,
    collision::PeriodicBoundary<NSCELL>>;
  using NSAllTasks = tmp::TupleWrapper<NSCollisionTask, NSPeriodicTask>;
  using NSSelector = tmp::TaskSelector<NSAllTasks, std::uint8_t, NSCELL>;

  using SurfaceTask = tmp::Key_TypePair<BulkFlag,
    ff::FFSurfaceTension2D<PFCELL, NSCELL>>;
  using GravityTask = tmp::Key_TypePair<BulkFlag,
    ff::FFGravityForce2D<PFCELL, NSCELL>>;
  using PressureTask = tmp::Key_TypePair<BulkFlag,
    RosenPressureForce2D<PFCELL, NSCELL>>;
  using RhoOmegaTask = tmp::Key_TypePair<BulkFlag,
    ff::FFRhoOmegaUpdate2D<PFCELL, NSCELL>>;
  BlockLatManagerCoupling surfaceCoupling(PF, NS);
  BlockLatManagerCoupling gravityCoupling(PF, NS);
  BlockLatManagerCoupling pressureCoupling(PF, NS);
  BlockLatManagerCoupling rhoOmegaCoupling(PF, NS);

  using MagneticCoefficientsTask = tmp::Key_TypePair<BulkFlag,
    MFUpdateCoeffs2D<PFCELL, MFCELL>>;
  using MagneticCoefficientsSelector = CoupledTaskSelector<
    std::uint8_t, PFCELL, MFCELL, MagneticCoefficientsTask>;
  BlockLatManagerCoupling magneticCoefficients(PF, MF);

  vtmo::ScalarWriter phiWriter("PHI", PF.getField<PHI<T>>());
  vtmo::ScalarWriter hmagWriter("HMAG", MF.getField<HMAG<T>>());
  vtmo::vtmWriter<T, 2> writer("rosensweig2d", geo);
  writer.addWriterSet(phiWriter, hmagWriter);

  // The framework periodic manager synchronizes populations, but per-cell
  // diagnostic fields need the same seam treatment. This helper handles both
  // same-rank and cross-rank edge blocks for arbitrary scalar/vector fields.
  auto syncPeriodicField = [&](auto& field) {
    using FieldBlock = std::decay_t<decltype(field.getBlockField(0))>;
    using Value = typename FieldBlock::value_type;
    constexpr int valueComponents =
      std::is_same_v<Value, Vector<T, 2>> ? 2 : 1;
    const T widthX = T(p.Ni) * p.CellLen;
    const int rank = mpi().getRank();
    const auto& globalGeometry = geoHelper.getBlockGeometry();
#ifdef MPI_ENABLED
    std::vector<std::vector<T>> sendBuffers;
    std::vector<std::vector<T>> receiveBuffers;
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> receiveRequests;
    std::vector<std::pair<int, int>> receiveJobs;
#endif

    auto writeValue = [&](std::vector<T>& buffer, int row, const Value& value) {
      if constexpr (valueComponents == 2) {
        buffer[row * valueComponents] = value[0];
        buffer[row * valueComponents + 1] = value[1];
      } else {
        buffer[row] = value;
      }
    };
    auto readValue = [&](const std::vector<T>& buffer, int row) {
      if constexpr (valueComponents == 2) {
        return Value{buffer[row * valueComponents],
                     buffer[row * valueComponents + 1]};
      } else {
        return Value{buffer[row]};
      }
    };

    for (int b = 0; b < geo.getBlockNum(); ++b) {
      const auto& block = geo.getBlock(b);
      const bool atLeft = block.getMin()[0] < p.CellLen * T{0.5};
      const bool atRight = block.getMax()[0] > widthX - p.CellLen * T{0.5};
      if (!atLeft && !atRight) continue;

      const int nx = block.getNx();
      const int ny = block.getNy();
      const int overlap = block.getOverlap();
      const auto& projection = block.getProjection();
      const T yMid = (block.getMin()[1] + block.getMax()[1]) / T{2};
      const T xProbe = atLeft ? widthX - p.CellLen * T{0.5}
                              : p.CellLen * T{0.5};
      int sourceGlobalId = -1;
      for (std::size_t candidate = 0;
           candidate < globalGeometry.getBlockNum(); ++candidate) {
        const auto& source = globalGeometry.getBlock(candidate);
        if (source.getSelfBlock().isInside(Vector<T, 2>{xProbe, yMid})) {
          sourceGlobalId = source.getBlockId();
          break;
        }
      }
      if (sourceGlobalId < 0) continue;

      const int ghostColumn = atLeft ? 0 : nx - overlap;
      const int sourceColumn = atLeft ? nx - 1 - overlap : overlap;
      const int sourceRank = geoHelper.whichRank(sourceGlobalId);
      if (sourceRank == rank) {
        const int sourceBlock = geo.findBlockIndex(sourceGlobalId);
        auto& targetField = field.getBlockField(b);
        const auto& sourceField = field.getBlockField(sourceBlock);
        for (int j = 0; j < ny; ++j) {
          const Value value = sourceField.get(j * projection[1] + sourceColumn);
          for (int c = 0; c < overlap; ++c) {
            targetField.get(j * projection[1] + ghostColumn + c) = value;
          }
        }
      } else {
#ifdef MPI_ENABLED
        const int physicalColumn = atLeft ? overlap : nx - 1 - overlap;
        std::vector<T> send(ny * valueComponents);
        auto& targetField = field.getBlockField(b);
        for (int j = 0; j < ny; ++j) {
          const Value value = targetField.get(j * projection[1] + physicalColumn);
          writeValue(send, j, value);
        }
        sendBuffers.emplace_back(std::move(send));
        MPI_Request sendRequest;
        mpi().iSend(sendBuffers.back().data(),
                    static_cast<int>(sendBuffers.back().size()), sourceRank,
                    &sendRequest, 9600 + block.getBlockId());
        sendRequests.push_back(sendRequest);

        receiveBuffers.emplace_back(ny * valueComponents);
        receiveJobs.emplace_back(b, ghostColumn);
        MPI_Request receiveRequest;
        mpi().iRecv(receiveBuffers.back().data(),
                    static_cast<int>(receiveBuffers.back().size()), sourceRank,
                    &receiveRequest, 9600 + sourceGlobalId);
        receiveRequests.push_back(receiveRequest);
#endif
      }
    }

#ifdef MPI_ENABLED
    MPI_Waitall(static_cast<int>(sendRequests.size()), sendRequests.data(),
                MPI_STATUSES_IGNORE);
    for (std::size_t i = 0; i < receiveRequests.size(); ++i) {
      MPI_Wait(&receiveRequests[i], MPI_STATUS_IGNORE);
      auto& targetField = field.getBlockField(receiveJobs[i].first);
      const auto& block = geo.getBlock(receiveJobs[i].first);
      const auto& projection = block.getProjection();
      const int overlap = block.getOverlap();
      const int ghostColumn = receiveJobs[i].second;
      const auto& buffer = receiveBuffers[i];
      for (int j = 0; j < block.getNy(); ++j) {
        const Value value = readValue(buffer, j);
        for (int c = 0; c < overlap; ++c) {
          targetField.get(j * projection[1] + ghostColumn + c) = value;
        }
      }
    }
#endif
  };

  auto setPsiWalls = [&]() {
    auto& field = MF.getField<PSI<T>>();
    const T wallDistance = T(p.PsiWallLayers) * p.CellLen;
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = MF.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockField = field.getBlockField(b);
      const int overlap = block.getOverlap();
      const T minX = block.getMin()[0];
      const T minY = block.getMin()[1];
      const T voxel = block.getVoxelSize();
      for (int j = 0; j < block.getNy(); ++j) {
        const T y = minY + T(j - overlap) * voxel;
        const bool nearBottom = y <= wallDistance && y >= -wallDistance;
        const bool nearTop = y <= domainHeight + wallDistance &&
          y >= domainHeight - wallDistance;
        if (!nearBottom && !nearTop) continue;
        for (int i = 0; i < block.getNx(); ++i) {
          const T x = minX + T(i) * voxel;
          const T psi = seedPsi(x, y);
          const std::size_t id = j * projection[1] + i;
          MFCELL cell(id, blockLat);
          for (unsigned k = 0; k < MFLatSet::q; ++k) {
            cell[k] = latset::w<MFLatSet>(k) * psi;
          }
          blockField.get(id) = psi;
        }
      }
    }
  };

  auto updatePsiFromPops = [&]() {
    auto& field = MF.getField<PSI<T>>();
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = MF.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockField = field.getBlockField(b);
      const int overlap = block.getOverlap();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          MFCELL cell(j * projection[1] + i, blockLat);
          T psi = T{0};
          for (unsigned k = 0; k < MFLatSet::q; ++k) psi += cell[k];
          blockField.get(j * projection[1] + i) = psi;
        }
      }
    }
  };

  auto computeMagneticField = [&]() {
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = MF.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      const int overlap = block.getOverlap();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          MFCELL cell(j * projection[1] + i, blockLat);
          MFComputeH2D<MFCELL>::apply(cell);
        }
      }
    }
    CommunicateAllMFFields<T>(MF);
    syncPeriodicField(MF.getField<HX<T>>());
    syncPeriodicField(MF.getField<HY<T>>());
    syncPeriodicField(MF.getField<HMAG<T>>());
  };

  auto copyPFYGhosts = [&]() {
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = PF.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      const int overlap = block.getOverlap();
      const int nx = block.getNx();
      const int ny = block.getNy();
      const T minY = block.getMin()[1];
      const T maxY = block.getMax()[1];
      if (minY < p.CellLen * T{1.5}) {
        for (int i = 0; i < nx; ++i) {
          PFCELL ghost(i, blockLat);
          PFCELL wall(overlap * projection[1] + i, blockLat);
          for (unsigned k = 0; k < LatSet::q; ++k) ghost[k] = wall[k];
        }
      }
      if (maxY > domainHeight - p.CellLen * T{1.5}) {
        for (int i = 0; i < nx; ++i) {
          PFCELL ghost((ny - 1) * projection[1] + i, blockLat);
          PFCELL wall((ny - 1 - overlap) * projection[1] + i, blockLat);
          for (unsigned k = 0; k < LatSet::q; ++k) ghost[k] = wall[k];
        }
      }
    }
  };

  auto copyNSYGhosts = [&]() {
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = NS.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      const int overlap = block.getOverlap();
      const int nx = block.getNx();
      const int ny = block.getNy();
      const T minY = block.getMin()[1];
      const T maxY = block.getMax()[1];
      if (minY < p.CellLen * T{1.5}) {
        for (int i = 0; i < nx; ++i) {
          NSCELL ghost(i, blockLat);
          NSCELL wall(overlap * projection[1] + i, blockLat);
          for (unsigned k = 0; k < LatSet::q; ++k) ghost[k] = wall[k];
        }
      }
      if (maxY > domainHeight - p.CellLen * T{1.5}) {
        for (int i = 0; i < nx; ++i) {
          NSCELL ghost((ny - 1) * projection[1] + i, blockLat);
          NSCELL wall((ny - 1 - overlap) * projection[1] + i, blockLat);
          for (unsigned k = 0; k < LatSet::q; ++k) ghost[k] = wall[k];
        }
      }
    }
  };

  auto setPhaseWalls = [&]() {
    auto& field = PF.getField<PHI<T>>();
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockField = field.getBlockField(b);
      const int overlap = block.getOverlap();
      const int nx = block.getNx();
      const int ny = block.getNy();
      const T minY = block.getMin()[1];
      const T maxY = block.getMax()[1];
      if (minY < p.CellLen * T{1.5}) {
        for (int i = 0; i < nx; ++i) {
          blockField.get(overlap * projection[1] + i) = T{1};
        }
      }
      if (maxY > domainHeight - p.CellLen * T{1.5}) {
        const int row = ny - 1 - overlap;
        for (int i = 0; i < nx; ++i) {
          blockField.get(row * projection[1] + i) = T{0};
        }
      }
      if (minY < p.CellLen * T{1.5}) {
        for (int j = 0; j < overlap; ++j) {
          for (int i = 0; i < nx; ++i) {
            blockField.get(j * projection[1] + i) = T{1};
          }
        }
      }
      if (maxY > domainHeight - p.CellLen * T{1.5}) {
        for (int j = ny - overlap; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            blockField.get(j * projection[1] + i) = T{0};
          }
        }
      }
    }
  };

  auto updatePhaseAndGradients = [&]() {
    auto& field = PF.getField<PHI<T>>();
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = PF.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockField = field.getBlockField(b);
      const int overlap = block.getOverlap();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          PFCELL cell(j * projection[1] + i, blockLat);
          T phi = T{0};
          for (unsigned k = 0; k < LatSet::q; ++k) phi += cell[k];
          phi = std::max(p.PhiFloor, std::min(p.PhiCeiling, phi));
          blockField.get(j * projection[1] + i) = phi;
        }
      }
    }
    field.Communicate();
    setPhaseWalls();
    field.Communicate();
    syncPeriodicField(field);

    PF.template ApplyInnerCellDynamics<PFNormalSelector>(flags);
    PF.template ApplyInnerCellDynamics<PFLaplacianSelector>(flags);
    PF.template ApplyInnerCellDynamics<PFChemicalSelector>(flags);
    PF.getField<NORMAL<T, 2>>().Communicate();
    PF.getField<GRAD<T, 2>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PF);

    auto& grad = PF.getField<GRAD<T, 2>>();
    auto& chemical = PF.getField<ff::CHEMICALPOTENTIAL<T>>();
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockGrad = grad.getBlockField(b);
      auto& blockChemical = chemical.getBlockField(b);
      const int overlap = block.getOverlap();
      const int nx = block.getNx();
      const T minY = block.getMin()[1];
      const T maxY = block.getMax()[1];
      if (minY < p.CellLen * T{1.5}) {
        for (int i = overlap; i < nx - overlap; ++i) {
          blockGrad.get(overlap * projection[1] + i)[1] =
            blockGrad.get((overlap + 1) * projection[1] + i)[1];
          blockChemical.get(overlap * projection[1] + i) =
            (T{4} * blockChemical.get((overlap + 1) * projection[1] + i) -
             blockChemical.get((overlap + 2) * projection[1] + i)) / T{3};
        }
      }
      if (maxY > domainHeight - p.CellLen * T{1.5}) {
        const int row = block.getNy() - 1 - overlap;
        for (int i = overlap; i < nx - overlap; ++i) {
          blockGrad.get(row * projection[1] + i)[1] =
            blockGrad.get((row - 1) * projection[1] + i)[1];
          blockChemical.get(row * projection[1] + i) =
            (T{4} * blockChemical.get((row - 1) * projection[1] + i) -
             blockChemical.get((row - 2) * projection[1] + i)) / T{3};
        }
      }
    }
    ff::CommunicateAllSelfFields<T>(PF);
  };

  auto updateNSMacros = [&]() {
    auto& density = NS.getField<DENSITY<T>>();
    auto& pressure = NS.getField<PRESSURE<T>>();
    auto& velocity = NS.getField<VELOCITY<T, 2>>();
    auto& force = NS.getField<FORCE<T, 2>>();
    for (int b = 0; b < geo.getBlockNum(); ++b) {
      auto& blockLat = NS.getBlockLat(b);
      const auto& block = geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockDensity = density.getBlockField(b);
      auto& blockPressure = pressure.getBlockField(b);
      auto& blockVelocity = velocity.getBlockField(b);
      auto& blockForce = force.getBlockField(b);
      const int overlap = block.getOverlap();
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          const std::size_t id = j * projection[1] + i;
          NSCELL cell(id, blockLat);
          T pressureValue = T{0};
          T ux = T{0};
          T uy = T{0};
          for (unsigned k = 0; k < LatSet::q; ++k) {
            pressureValue += cell[k];
            ux += latset::c<LatSet>(k)[0] * cell[k];
            uy += latset::c<LatSet>(k)[1] * cell[k];
          }
          const T rho = blockDensity.get(id);
          const auto bodyForce = blockForce.get(id);
          blockPressure.get(id) = pressureValue;
          blockVelocity.get(id) = Vector<T, 2>{
            ux + bodyForce[0] / (T{2} * rho),
            uy + bodyForce[1] / (T{2} * rho)};
        }
      }
    }
    density.Communicate();
    pressure.Communicate();
    velocity.Communicate();
    force.Communicate();
    syncPeriodicField(velocity);
    syncPeriodicField(force);
  };

  PF.NormalFullCommunicate();
  NS.NormalFullCommunicate();
  MF.NormalFullCommunicate();
  NSPeriodic.Apply();
  PFPeriodic.Apply();
  MFPeriodic.Apply();
  PF.template ApplyInnerCellDynamics<PFNormalSelector>(flags);
  PF.template ApplyInnerCellDynamics<PFLaplacianSelector>(flags);
  PF.template ApplyInnerCellDynamics<PFChemicalSelector>(flags);
  PF.getField<NORMAL<T, 2>>().Communicate();
  PF.getField<GRAD<T, 2>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PF);
  magneticCoefficients.ApplyInnerCellDynamics<MagneticCoefficientsSelector>(0, flags);
  CommunicateOMEGAPSI<T>(MF);
  syncPeriodicField(MF.getField<PSI<T>>());
  computeMagneticField();
  updateNSMacros();
  writer.WriteBinary(0);

  Printer::Print_BigBanner(std::string("Start Rosensweig calculation..."));
  Timer timer;
  Timer outputTimer;
  const T psiWallHeight = domainHeight;

  while (timer() < p.MaxStep) {
    magneticCoefficients.ApplyInnerCellDynamics<MagneticCoefficientsSelector>(
      timer(), flags);
    CommunicateOMEGAPSI<T>(MF);

    for (int sub = 0; sub < p.PsiSolverIter; ++sub) {
      setPsiWalls();
      for (int b = 0; b < geo.getBlockNum(); ++b) {
        auto& blockLat = MF.getBlockLat(b);
        const auto& block = geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            MFCELL cell(j * projection[1] + i, blockLat);
            collision::MRTDiffusion<MFCELL, OMEGA_PSI<T>>::apply(cell);
          }
        }
      }
      MFPeriodic.Apply();
      MF.NormalFullCommunicate();
      MF.Stream();
      MF.NormalFullCommunicate();
      updatePsiFromPops();
      CommunicatePSI<T>(MF);
      setPsiWalls();
    }
    syncPeriodicField(MF.getField<PSI<T>>());
    computeMagneticField();

    rhoOmegaCoupling.ApplyInnerCellDynamics<
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RhoOmegaTask>>(
        timer(), flags);
    NS.getField<FORCE<T, 2>>().InitValue(Vector<T, 2>{0, 0});
    surfaceCoupling.ApplyInnerCellDynamics<
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, SurfaceTask>>(
        timer(), flags);

    if (p.H0 > T{0}) {
      for (int b = 0; b < geo.getBlockNum(); ++b) {
        auto& pfBlock = PF.getBlockLat(b);
        auto& mfBlock = MF.getBlockLat(b);
        auto& nsBlock = NS.getBlockLat(b);
        const auto& block = geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int overlap = block.getOverlap();
        const T minY = block.getMin()[1];
        const T voxel = block.getVoxelSize();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          const T y = minY + T(j - overlap) * voxel;
          if (y <= T(p.MagneticForceWallBand) * p.CellLen ||
              y >= psiWallHeight - T(p.MagneticForceWallBand) * p.CellLen) {
            continue;
          }
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            const std::size_t id = j * projection[1] + i;
            PFCELL pfCell(id, pfBlock);
            MFCELL mfCell(id, mfBlock);
            NSCELL nsCell(id, nsBlock);
            MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
              pfCell, mfCell, nsCell);
          }
        }
      }
    }

    gravityCoupling.ApplyInnerCellDynamics<
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GravityTask>>(
        timer(), flags);
    pressureCoupling.ApplyInnerCellDynamics<
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PressureTask>>(
        timer(), flags);
    NS.getField<FORCE<T, 2>>().Communicate();

    PF.template ApplyInnerCellDynamics<PFSelector>(flags);
    PFPeriodic.Apply();
    PF.NormalFullCommunicate();
    NS.template ApplyInnerCellDynamics<NSSelector>(flags);
    NSPeriodic.Apply();
    NS.NormalFullCommunicate();

    PFBounce.Apply(timer());
    copyPFYGhosts();
    PF.Stream();
    PF.NormalFullCommunicate();

    NSBounce.Apply(timer());
    copyNSYGhosts();
    NS.Stream();
    NS.NormalFullCommunicate();

    updatePhaseAndGradients();
    updateNSMacros();

    ++timer;
    ++outputTimer;
    if (timer() % p.OutputStep == 0) {
      outputTimer.Print_InnerLoopPerformance(geo.getTotalCellNum(), p.OutputStep);
      Printer::Endl();
      writer.WriteBinary(timer());
    }
  }

  Printer::Print_BigBanner(std::string("Rosensweig calculation complete"));
  timer.Print_MainLoopPerformance(geo.getTotalCellNum());
  return 0;
}
