// rosen3d.cpp - 3D Rosensweig instability (phase field + NS + magnetic)
#include "freelb.h"
#include "freelb.hh"
#include "ff/ff2d.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

using T = FLOAT;
using LatSet = D3Q19<T>;
using MFLatSet = D3Q7<T>;
using namespace mfield;

// ---- Simulation parameters ----
int Lattice_L{};
int Ni{}, Nj{}, Nk{};
T Cell_Len{};
int BlockCellLen{}, Thread_Num{};

// Rosensweig geometry
T Interface_Width{}, Mobility{};
T Wavelength{}, InterfaceZ{}, PerturbationAmplitude{};

// Phase field and two-phase parameters
T Tau_phi{}, Omega_phi{}, Kappa{}, Beta{};
T rho_l{}, rho_h{}, eta_l{}, eta_h{}, sigma{}, gravity{}, DeltaRho{}, Tau_ns{};

// Magnetic field parameters
T mu_l{}, mu_h{}, chi_l{}, chi_h{};
T H0{}, Hc{}, H0_over_Hc{};
T PsiSolver_K{};
int PsiSolver_Iter{};
T PsiSolver_Residual{};

int MaxStep{}, OutputStep{};
std::string work_dir;

void readParam(const std::string& filename) {
  iniReader r(filename);
  work_dir = r.getValue<std::string>("workdir", "workdir_");
  Thread_Num = r.getValue<int>("parallel", "thread_num");

  Lattice_L = r.getValue<int>("Mesh", "L");
  Cell_Len = r.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = r.getValue<int>("Mesh", "BlockCellLen");
  Ni = 3 * Lattice_L;
  Nj = 3 * Lattice_L;
  Nk = Lattice_L;

  Interface_Width = r.getValue<T>("Phase_Field", "W");
  Mobility = r.getValue<T>("Phase_Field", "M");

  rho_l = r.getValue<T>("Two_Phase", "rho_l");
  rho_h = r.getValue<T>("Two_Phase", "rho_h");
  eta_l = r.getValue<T>("Two_Phase", "eta_l");
  eta_h = r.getValue<T>("Two_Phase", "eta_h");
  sigma = r.getValue<T>("Two_Phase", "sigma");

  mu_l = r.getValue<T>("Magnetic_Field", "mu_l");
  mu_h = r.getValue<T>("Magnetic_Field", "mu_h");
  chi_l = r.getValue<T>("Magnetic_Field", "chi_l");
  chi_h = r.getValue<T>("Magnetic_Field", "chi_h");
  H0_over_Hc = r.getValue<T>("Magnetic_Field", "H0_over_Hc");
  PsiSolver_Iter = r.getValue<int>("Magnetic_Field", "PsiSolver_Iter");
  PsiSolver_K = r.getValue<T>("Magnetic_Field", "PsiSolver_K");
  PsiSolver_Residual = r.getValue<T>("Magnetic_Field", "PsiSolver_Residual");

  MaxStep = r.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = r.getValue<int>("Simulation_Settings", "OutputStep");

  const T pi = std::acos(T{-1});
  Wavelength = T{0.5} * T(Lattice_L) * Cell_Len;
  InterfaceZ = T{1.0} / T{3.0} * T(Lattice_L) * Cell_Len;
  PerturbationAmplitude = T{0.01} * Wavelength;
  DeltaRho = rho_h - rho_l;

  // lambda_c = 2*pi*sqrt(sigma/(g*DeltaRho)).
  const T wavelength_lattice = T{0.5} * T(Lattice_L);
  gravity = T{4} * pi * pi * sigma / (DeltaRho * wavelength_lattice * wavelength_lattice);

  // Cowley-Rosensweig critical field, paper Eq. (71), with mu0=1 in LBM units.
  const T mu0 = T{1};
  const T mu_ratio = mu0 / mu_h;
  const T magnetic_factor = (T{2} / mu0) * (mu_ratio + T{1}) /
                           ((mu_ratio - T{1}) * (mu_ratio - T{1}));
  Hc = std::sqrt(magnetic_factor) * std::pow(sigma * gravity * DeltaRho, T{0.25});
  H0 = H0_over_Hc * Hc;

  Beta = T{12} * sigma / Interface_Width;
  Kappa = T{1.5} * Interface_Width * sigma;
  Tau_phi = T{3} * Mobility + T{0.5};
  Omega_phi = T{1} / Tau_phi;
  Tau_ns = T{0.5} + eta_h / rho_h / LatSet::cs2;

  MPI_RANK(0) {
    std::printf("---- 3D Rosensweig Instability ----\n");
    std::printf("Mesh: %dx%dx%d (L=%d) BlockCellLen=%d\n",
                Ni, Nj, Nk, Lattice_L, BlockCellLen);
    std::printf("Interface: z0=%.4f lambda=%.4f A=%.4f W=%.4f\n",
                InterfaceZ, Wavelength, PerturbationAmplitude, Interface_Width);
    std::printf("rho: l=%.6f h=%.6f  eta: l=%.6f h=%.6f  sigma=%.6f g=%.6e\n",
                rho_l, rho_h, eta_l, eta_h, sigma, gravity);
    std::printf("M=%.6f tau_phi=%.6f  DeltaRho=%.6f\n",
                Mobility, Tau_phi, DeltaRho);
    std::printf("mu=(%.6f,%.6f) chi=(%.6f,%.6f) Hc=%.6e H0/Hc=%.6f H0=%.6e\n",
                mu_l, mu_h, chi_l, chi_h, Hc, H0_over_Hc, H0);
    std::printf("PsiSolver: K=%.6f max_iter=%d residual=%.3e\n",
                PsiSolver_K, PsiSolver_Iter, PsiSolver_Residual);
    std::printf("-----------------------------------\n");
  }
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = 1;
  constexpr std::uint8_t BulkFlag = 2;
  constexpr std::uint8_t BouncebackFlag = 4;
  constexpr std::uint8_t PeriodicFlag = 8;

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT
  Printer::Print_BigBanner(std::string("Initializing 3D Rosensweig Instability..."));
  readParam(argc > 1 ? std::string(argv[1]) : std::string("rosen3d.ini"));

  // -- converters --
  BaseConverter<T> BaseConv(LatSet::cs2);
  BaseConv.SimplifiedConverterFromRT(Ni, T{0.01}, Tau_ns);
  BaseConverter<T> PFBaseConv(LatSet::cs2);
  PFBaseConv.SimplifiedConverterFromRT(Ni, T{0.01}, Tau_phi);
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(Ni, T{0.01}, T{1.0});
  UnitConvManager<T> ConvManager(&BaseConv);
  ConvManager.Check_and_Print();

  // -- geometry --
  // The physical domain is 3L x 3L x L. X/Y are periodic and Z is solid.
  AABB<T, 3> domain({0, 0, 0},
                    {T(Ni) * Cell_Len, T(Nj) * Cell_Len, T(Nk) * Cell_Len});
  AABB<T, 3> left({-Cell_Len, 0, 0},
                  {0, T(Nj) * Cell_Len, T(Nk) * Cell_Len});
  AABB<T, 3> right({T(Ni) * Cell_Len, 0, 0},
                   {T(Ni + 1) * Cell_Len, T(Nj) * Cell_Len, T(Nk) * Cell_Len});
  AABB<T, 3> front({0, -Cell_Len, 0},
                   {T(Ni) * Cell_Len, 0, T(Nk) * Cell_Len});
  AABB<T, 3> back({0, T(Nj) * Cell_Len, 0},
                  {T(Ni) * Cell_Len, T(Nj + 1) * Cell_Len, T(Nk) * Cell_Len});

  BlockGeometryHelper3D<T> GeoHelper(
      Ni, Nj, Nk, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(2, 2, 4);
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry3D<T> Geo(GeoHelper);

  // The periodic magnetic scalar fields need an explicit edge-plane exchange.
  // The normal block communicator handles block neighbors, while this table
  // also identifies the x/y partner when it is owned by another MPI rank.
  std::vector<double> GlobalBlockTable;
  {
    const int nranks = mpi().getSize();
    const T xLeft = T{0};
    const T xRight = T(Ni) * Cell_Len;
    const T yFront = T{0};
    const T yBack = T(Nj) * Cell_Len;
    std::vector<double> local;
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& block = Geo.getBlock(b);
      local.push_back(double(mpi().getRank()));
      local.push_back(block.getMin()[0]);
      local.push_back(block.getMin()[1]);
      local.push_back(block.getMin()[2]);
      local.push_back(block.getNx());
      local.push_back(block.getNy());
      local.push_back(block.getNz());
      local.push_back(block.getMin()[0] < xLeft + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(block.getMax()[0] > xRight - Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(block.getMin()[1] < yFront + Cell_Len * T{0.5} ? 1.0 : 0.0);
      local.push_back(block.getMax()[1] > yBack - Cell_Len * T{0.5} ? 1.0 : 0.0);
    }

    std::vector<int> counts(nranks, 0);
    {
      double send_size[1]{};
      double recv_size[1]{};
      send_size[0] = double(local.size());
      if (mpi().getRank() == 0) {
        for (int rank = 0; rank < nranks; ++rank) {
          mpi().sendRecv(send_size, recv_size, 1, rank, rank, 1000);
          counts[rank] = int(recv_size[0]);
        }
      } else {
        mpi().sendRecv(send_size, recv_size, 1, 0, 0, 1000);
      }
    }
    mpi().bCast(counts.data(), nranks, 0);

    int total = 0;
    for (int count : counts) total += count;
    std::vector<double> all(total, 0.0);
    {
      std::vector<double> recv;
      if (mpi().getRank() == 0) {
        std::copy(local.begin(), local.end(), all.begin());
        int offset = counts[0];
        for (int rank = 1; rank < nranks; ++rank) {
          recv.assign(counts[rank], 0.0);
          mpi().sendRecv(recv.data(), all.data() + offset,
                         counts[rank], rank, rank, 1001);
          offset += counts[rank];
        }
      } else {
        recv.assign(counts[mpi().getRank()], 0.0);
        mpi().sendRecv(local.data(), recv.data(), int(local.size()), 0, 0, 1001);
      }
    }
    mpi().bCast(all.data(), total, 0);
    GlobalBlockTable = std::move(all);
  }

  // -- flags --
  BlockFieldManager<FLAG, T, 3> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, BulkFlag);
  });
  FlagFM.forEach(left, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  FlagFM.forEach(right, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  FlagFM.forEach(front, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  FlagFM.forEach(back, [&](FLAG& flag, std::size_t id) {
    flag.SetField(id, PeriodicFlag);
  });
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // -- D3Q19 NS lattice --
  using NSFIELDS = TypePack<DENSITY<T>, VELOCITY<T, 3>, POP<T, LatSet::q>,
                            FORCE<T, 3>, OMEGA<T>, PRESSURE<T>>;
  const T omega_ns = T{1} / Tau_ns;
  ValuePack NSI(T{1}, Vector<T, 3>{0, 0, 0}, T{}, Vector<T, 3>{0, 0, 0},
                omega_ns, T{});
  using NSCELL = Cell<T, LatSet, NSFIELDS>;
  BlockLatticeManager<T, LatSet, NSFIELDS> NSLattice(Geo, NSI, BaseConv);

  // -- D3Q19 phase-field lattice --
  using PFFIELDS = TypePack<
      PHI<T>, POP<T, LatSet::q>, GRAD<T, 3>, NORMAL<T, 3>, INTERFACEWIDTH<T>,
      ff::LAPLACIAN<T>, ff::CHEMICALPOTENTIAL<T>, ff::GRAVITY<T>, ff::BETA<T>,
      ff::KAPPA<T>, ff::RHO_L<T>, ff::RHO_H<T>, ff::ETA_L<T>, ff::ETA_H<T>,
      ff::DELTARHO<T>>;
  using PFREF = TypePack<VELOCITY<T, 3>>;
  using PFPACK = TypePack<PFFIELDS, PFREF>;
  ValuePack PFI(T{}, T{}, Vector<T, 3>{0, 0, 0}, Vector<T, 3>{0, 0, 0},
                Interface_Width, T{}, T{}, gravity, Beta, Kappa, rho_l, rho_h,
                eta_l, eta_h, DeltaRho);
  using PFCELL = Cell<T, LatSet, ExtractFieldPack<PFPACK>::mergedpack>;
  BlockLatticeManager<T, LatSet, PFPACK> PFLattice(
      Geo, PFI, PFBaseConv, &NSLattice.getField<VELOCITY<T, 3>>());
  ff::BroadcastAllParams<T>(PFLattice, rho_l, rho_h, eta_l, eta_h, gravity,
                            Beta, Kappa);
  PFLattice.getField<ff::DELTARHO<T>>().InitValue(DeltaRho);

  // -- D3Q7 magnetic lattice --
  using MFFIELDS = TypePack<
      PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>, HX<T>, HY<T>, HZ<T>,
      HMAG<T>, POP<T, MFLatSet::q>, MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>,
      H_0<T>, PSI_K<T>>;
  using MFREF = TypePack<PHI<T>>;
  using MFPACK = TypePack<MFFIELDS, MFREF>;
  ValuePack MFI(T{}, T{1.0}, T{mu_l}, T{chi_l}, T{}, T{}, T{}, T{}, T{},
                mu_l, mu_h, chi_l, chi_h, H0, PsiSolver_K);
  using MFCELL = Cell<T, MFLatSet, ExtractFieldPack<MFPACK>::mergedpack>;
  BlockLatticeManager<T, MFLatSet, MFPACK> MFLattice(
      Geo, MFI, MFBaseConv, &PFLattice.getField<PHI<T>>());
  BroadcastAllMFParams3D<T>(MFLattice, mu_l, mu_h, chi_l, chi_h, H0,
                            PsiSolver_K);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(T{1.0});

  // -- initial phase field: ferrofluid below a perturbed interface --
  const T two_pi = T{2} * std::acos(T{-1});
  auto& phiField = PFLattice.getField<PHI<T>>();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    const auto& block = Geo.getBlock(b);
    const auto& projection = block.getProjection();
    auto& blockPhi = phiField.getBlockField(b);
    const T voxelSize = block.getVoxelSize();
    const T minX = block.getMin()[0];
    const T minY = block.getMin()[1];
    const T minZ = block.getMin()[2];
    for (int k = 0; k < block.getNz(); ++k) {
      const T z = minZ + T(k) * voxelSize;
      for (int j = 0; j < block.getNy(); ++j) {
        const T y = minY + T(j) * voxelSize;
        for (int i = 0; i < block.getNx(); ++i) {
          const T x = minX + T(i) * voxelSize;
          const T interfaceZ = InterfaceZ + PerturbationAmplitude *
              std::cos(two_pi * x / Wavelength) *
              std::cos(two_pi * y / Wavelength);
          const T signedDistance = interfaceZ - z;
          const T phi = T{0.5} + T{0.5} *
              std::tanh(T{2} * signedDistance / (Interface_Width * Cell_Len));
          blockPhi.get(k * projection[2] + j * projection[1] + i) = phi;
        }
      }
    }
  }

  // Initialize phase-field populations from phi and NS populations at rest.
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& blockLat = PFLattice.getBlockLat(b);
    auto& blockPhi = phiField.getBlockField(b);
    const auto& block = Geo.getBlock(b);
    const auto& projection = block.getProjection();
    for (int k = 0; k < block.getNz(); ++k) {
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          const std::size_t id = k * projection[2] + j * projection[1] + i;
          PFCELL cell(id, blockLat);
          const T phi = blockPhi.get(id);
          for (unsigned int q = 0; q < LatSet::q; ++q) {
            cell[q] = latset::w<LatSet>(q) * phi;
          }
        }
      }
    }
  }
  PFLattice.getField<INTERFACEWIDTH<T>>().InitValue(Interface_Width);

  const Vector<T, 3> zeroVelocity{0, 0, 0};
  const T zeroPressure = T{0};
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& blockLat = NSLattice.getBlockLat(b);
    const auto& block = Geo.getBlock(b);
    const auto& projection = block.getProjection();
    for (int k = 0; k < block.getNz(); ++k) {
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          const std::size_t id = k * projection[2] + j * projection[1] + i;
          NSCELL cell(id, blockLat);
          for (unsigned int q = 0; q < LatSet::q; ++q) {
            const T cu = zeroVelocity * latset::c<LatSet>(q);
            cell[q] = latset::w<LatSet>(q) *
                      (zeroPressure + LatSet::InvCs2 * cu +
                       T{0.5} * LatSet::InvCs4 * cu * cu - LatSet::InvCs2);
          }
        }
      }
    }
  }

  // Initial magnetic potential for a uniform +z field: H=-grad(psi)=H0 ez.
  auto& psiField = MFLattice.getField<PSI<T>>();
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& blockLat = MFLattice.getBlockLat(b);
    const auto& block = Geo.getBlock(b);
    const auto& projection = block.getProjection();
    auto& blockPsi = psiField.getBlockField(b);
    const T voxelSize = block.getVoxelSize();
    const T minZ = block.getMin()[2];
    for (int k = 0; k < block.getNz(); ++k) {
      const T z = minZ + T(k) * voxelSize;
      const T psi = -H0 * z;
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          const std::size_t id = k * projection[2] + j * projection[1] + i;
          MFCELL cell(id, blockLat);
          for (unsigned int q = 0; q < MFLatSet::q; ++q) {
            cell[q] = latset::w<MFLatSet>(q) * psi;
          }
          blockPsi.get(id) = psi;
        }
      }
    }
  }

  // -- boundary managers --
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFPACK>;
  using LM_MF = BlockLatticeManager<T, MFLatSet, MFPACK>;
  using FM = BlockFieldManager<FLAG, T, 3>;

  BBLikeFixedBlockBdManager<bounceback::normal<NSCELL>, LM_NS, FM>
      NS_BB("NS_BB", NSLattice, FlagFM, BouncebackFlag, VoidFlag);
  BBLikeFixedBlockBdManager<bounceback::normal<PFCELL>, LM_PF, FM>
      PF_BB("PF_BB", PFLattice, FlagFM, BouncebackFlag, VoidFlag);

  FixedPeriodicBoundaryManager<LM_NS, FM> NS_Per(
      "NS_Per", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_Per.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  NS_Per.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
  NS_Per.Setup(front, NbrDirection::YN, back, NbrDirection::YP);
  NS_Per.Setup(back, NbrDirection::YP, front, NbrDirection::YN);

  FixedPeriodicBoundaryManager<LM_PF, FM> PF_Per(
      "PF_Per", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_Per.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  PF_Per.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
  PF_Per.Setup(front, NbrDirection::YN, back, NbrDirection::YP);
  PF_Per.Setup(back, NbrDirection::YP, front, NbrDirection::YN);

  FixedPeriodicBoundaryManager<LM_MF, FM> MF_Per(
      "MF_Per", MFLattice, FlagFM, PeriodicFlag, VoidFlag);
  MF_Per.Setup(left, NbrDirection::XN, right, NbrDirection::XP);
  MF_Per.Setup(right, NbrDirection::XP, left, NbrDirection::XN);
  MF_Per.Setup(front, NbrDirection::YN, back, NbrDirection::YP);
  MF_Per.Setup(back, NbrDirection::YP, front, NbrDirection::YN);
#ifdef MPI_ENABLED
  NS_Per.SetupMPI(GeoHelper);
  PF_Per.SetupMPI(GeoHelper);
  MF_Per.SetupMPI(GeoHelper);
#endif

  // -- phase-field, NS, and coupling tasks --
  using PFNT = tmp::Key_TypePair<BulkFlag, ff::FF3D<PFCELL>>;
  using PFLT = tmp::Key_TypePair<BulkFlag, ff::FFLaplacian3D<PFCELL>>;
  using PFCT = tmp::Key_TypePair<BulkFlag, ff::FFChemPotential3D<PFCELL>>;
  using PFSelN = TaskSelector<std::uint8_t, PFCELL, PFNT>;
  using PFSelL = TaskSelector<std::uint8_t, PFCELL, PFLT>;
  using PFSelC = TaskSelector<std::uint8_t, PFCELL, PFCT>;
  using PFColT = tmp::Key_TypePair<
      BulkFlag,
      collision::MRTSource<equilibrium::FirstOrder<PFCELL>, NORMAL<T, 3>, true, true>>;
  using PFPerT = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<PFCELL>>;
  using PFAll = tmp::TupleWrapper<PFColT, PFPerT>;
  using PFSel = tmp::TaskSelector<PFAll, std::uint8_t, PFCELL>;

  using NSMT = tmp::Key_TypePair<BulkFlag, collision::MRTForce<NSCELL, FORCE<T, 3>>>;
  using NSPT = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<NSCELL>>;
  using NSAll = tmp::TupleWrapper<NSMT, NSPT>;
  using NSSel = tmp::TaskSelector<NSAll, std::uint8_t, NSCELL>;

  using STT = tmp::Key_TypePair<BulkFlag, ff::FFSurfaceTension3D<PFCELL, NSCELL>>;
  BlockLatManagerCoupling STC(PFLattice, NSLattice);
  using GrT = tmp::Key_TypePair<BulkFlag, ff::FFGravityForce3D<PFCELL, NSCELL>>;
  BlockLatManagerCoupling GrC(PFLattice, NSLattice);
  using PrT = tmp::Key_TypePair<BulkFlag, ff::FFPreForce3D<PFCELL, NSCELL>>;
  BlockLatManagerCoupling PrC(PFLattice, NSLattice);
  using ViT = tmp::Key_TypePair<BulkFlag, ff::FFViscoForce3DM<PFCELL, NSCELL>>;
  BlockLatManagerCoupling ViC(PFLattice, NSLattice);
  using RoT = tmp::Key_TypePair<BulkFlag, ff::FFRhoOmegaUpdate3D<PFCELL, NSCELL>>;
  BlockLatManagerCoupling RoC(PFLattice, NSLattice);

  using MCT = tmp::Key_TypePair<BulkFlag, MFUpdateCoeffs3D<PFCELL, MFCELL>>;
  using MCSel = CoupledTaskSelector<std::uint8_t, PFCELL, MFCELL, MCT>;
  BlockLatManagerCoupling MCC(PFLattice, MFLattice);

  // -- output --
  vtmo::ScalarWriter PW("PHI", PFLattice.getField<PHI<T>>());
  vtmo::ScalarWriter PS("PSI", MFLattice.getField<PSI<T>>());
  vtmo::VectorWriter VW("Velocity", NSLattice.getField<VELOCITY<T, 3>>());
  vtmo::ScalarWriter Dw("Density", NSLattice.getField<DENSITY<T>>());
  vtmo::VectorWriter Fw("Force", NSLattice.getField<FORCE<T, 3>>());
  vtmo::ScalarWriter Hxw("HX", MFLattice.getField<HX<T>>());
  vtmo::ScalarWriter Hyw("HY", MFLattice.getField<HY<T>>());
  vtmo::ScalarWriter Hzw("HZ", MFLattice.getField<HZ<T>>());
  vtmo::ScalarWriter Hmw("HMAG", MFLattice.getField<HMAG<T>>());
  vtmo::vtmWriter<T, 3> MW("rosen3d", Geo);
  MW.addWriterSet(PW, PS, VW, Dw, Fw, Hxw, Hyw, Hzw, Hmw);

  // -- initial communications and fields --
  PFLattice.NormalFullCommunicate();
  NSLattice.NormalFullCommunicate();
  MFLattice.NormalFullCommunicate();
  NS_Per.Apply();
  PF_Per.Apply();
  MF_Per.Apply();

  PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
  PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
  PFLattice.getField<NORMAL<T, 3>>().Communicate();
  PFLattice.getField<GRAD<T, 3>>().Communicate();
  ff::CommunicateAllSelfFields<T>(PFLattice);
  RoC.ApplyInnerCellDynamics<
      CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RoT>>(0, FlagFM);

  // Populate the initial magnetic output fields before writing T0.  The
  // first physical step performs the full periodic ghost synchronization;
  // this initial pass is sufficient for the interior diagnostic snapshot.
  for (int b = 0; b < Geo.getBlockNum(); ++b) {
    auto& blockLat = MFLattice.getBlockLat(b);
    const auto& block = Geo.getBlock(b);
    const auto& projection = block.getProjection();
    const int overlap = block.getOverlap();
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          MFCELL cell(k * projection[2] + j * projection[1] + i, blockLat);
          MFComputeH3D<MFCELL>::apply(cell);
        }
      }
    }
  }
  CommunicateAllMFFields3D<T>(MFLattice);
  MW.WriteBinary(0);

  // Pinning the scalar potential on the first few z rows supplies the uniform
  // applied field and keeps the D3Q7 wall treatment consistent with the 3D
  // magnetic solver. The interior populations are warm-started between steps.
  const T H_global = T(Nk) * Cell_Len;
  auto PinPsiWalls = [&]() {
    auto& field = MFLattice.getField<PSI<T>>();
    const T wallBand = Cell_Len * T{3};
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& blockLat = MFLattice.getBlockLat(b);
      const auto& block = Geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& blockPsi = field.getBlockField(b);
      const int nx = block.getNx();
      const int ny = block.getNy();
      const int nz = block.getNz();
      const int overlap = block.getOverlap();
      const T minZ = block.getMin()[2];
      const T voxelSize = block.getVoxelSize();
      for (int k = 0; k < nz; ++k) {
        const T z = minZ + T(k - overlap) * voxelSize;
        if (!((z <= wallBand && z >= -wallBand) ||
              (z <= H_global + wallBand && z >= H_global - wallBand))) {
          continue;
        }
        const T wallPsi = -H0 * z;
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            const std::size_t id = k * projection[2] + j * projection[1] + i;
            MFCELL cell(id, blockLat);
            for (unsigned int q = 0; q < MFLatSet::q; ++q) {
              cell[q] = latset::w<MFLatSet>(q) * wallPsi;
            }
            blockPsi.get(id) = wallPsi;
          }
        }
      }
    }
  };

  auto SyncMFPeriodicGhosts = [&](auto& field, int fieldIndex) {
    const T xLeft = T{0};
    const T xRight = T(Ni) * Cell_Len;
    const T yFront = T{0};
    const T yBack = T(Nj) * Cell_Len;
    const int entryCount = int(GlobalBlockTable.size() / 11);
    const int myRank = mpi().getRank();

    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      const auto& block = Geo.getBlock(b);
      const bool atLeft = block.getMin()[0] < xLeft + Cell_Len * T{0.5};
      const bool atRight = block.getMax()[0] > xRight - Cell_Len * T{0.5};
      const bool atFront = block.getMin()[1] < yFront + Cell_Len * T{0.5};
      const bool atBack = block.getMax()[1] > yBack - Cell_Len * T{0.5};
      if (!atLeft && !atRight && !atFront && !atBack) continue;

      const auto& projection = block.getProjection();
      const int nx = block.getNx();
      const int ny = block.getNy();
      const int nz = block.getNz();
      auto& blockField = field.getBlockField(b);

      if (atLeft || atRight) {
        int localEntry = -1;
        int partnerEntry = -1;
        for (int e = 0; e < entryCount; ++e) {
          const double* entry = &GlobalBlockTable[e * 11];
          if (entry[2] != block.getMin()[1] || entry[3] != block.getMin()[2]) continue;
          if (atLeft && entry[8] > 0.5) partnerEntry = e;
          if (atRight && entry[7] > 0.5) partnerEntry = e;
          if (int(entry[0]) == myRank && entry[1] == block.getMin()[0] &&
              entry[2] == block.getMin()[1] && entry[3] == block.getMin()[2]) {
            localEntry = e;
          }
        }
        if (localEntry >= 0 && partnerEntry >= 0) {
          const double* partner = &GlobalBlockTable[partnerEntry * 11];
          const int partnerRank = int(partner[0]);
          const int plane = ny * nz;
          if (partnerRank == myRank) {
            int partnerBlock = -1;
            for (int bp = 0; bp < Geo.getBlockNum(); ++bp) {
              const auto& candidate = Geo.getBlock(bp);
              if (candidate.getMin()[0] == partner[1] &&
                  candidate.getMin()[1] == partner[2] &&
                  candidate.getMin()[2] == partner[3]) {
                partnerBlock = bp;
                break;
              }
            }
            if (partnerBlock >= 0) {
              auto& partnerField = field.getBlockField(partnerBlock);
              const int partnerNx = Geo.getBlock(partnerBlock).getNx();
              if (atLeft) {
                for (int k = 0; k < nz; ++k) {
                  for (int j = 0; j < ny; ++j) {
                    blockField.get(k * projection[2] + j * projection[1]) =
                        partnerField.get(k * projection[2] + j * projection[1] + partnerNx - 2);
                  }
                }
              }
              if (atRight) {
                for (int k = 0; k < nz; ++k) {
                  for (int j = 0; j < ny; ++j) {
                    blockField.get(k * projection[2] + j * projection[1] + nx - 1) =
                        partnerField.get(k * projection[2] + j * projection[1] + 1);
                  }
                }
              }
            }
          } else {
            const int tag = fieldIndex * 1000 + std::min(localEntry, partnerEntry) * 2;
            const int sourceI = atLeft ? 1 : nx - 2;
            const int targetI = atLeft ? 0 : nx - 1;
            std::vector<T> sendBuffer(plane);
            std::vector<T> receiveBuffer(plane);
            int index = 0;
            for (int k = 0; k < nz; ++k) {
              for (int j = 0; j < ny; ++j) {
                sendBuffer[index++] = blockField.get(k * projection[2] +
                                                    j * projection[1] + sourceI);
              }
            }
            mpi().sendRecv(sendBuffer.data(), receiveBuffer.data(), plane,
                           partnerRank, partnerRank, tag);
            index = 0;
            for (int k = 0; k < nz; ++k) {
              for (int j = 0; j < ny; ++j) {
                blockField.get(k * projection[2] + j * projection[1] + targetI) =
                    receiveBuffer[index++];
              }
            }
          }
        }
      }

      if (atFront || atBack) {
        int localEntry = -1;
        int partnerEntry = -1;
        for (int e = 0; e < entryCount; ++e) {
          const double* entry = &GlobalBlockTable[e * 11];
          if (entry[1] != block.getMin()[0] || entry[3] != block.getMin()[2]) continue;
          if (atFront && entry[10] > 0.5) partnerEntry = e;
          if (atBack && entry[9] > 0.5) partnerEntry = e;
          if (int(entry[0]) == myRank && entry[1] == block.getMin()[0] &&
              entry[2] == block.getMin()[1] && entry[3] == block.getMin()[2]) {
            localEntry = e;
          }
        }
        if (localEntry >= 0 && partnerEntry >= 0) {
          const double* partner = &GlobalBlockTable[partnerEntry * 11];
          const int partnerRank = int(partner[0]);
          const int plane = nx * nz;
          if (partnerRank == myRank) {
            int partnerBlock = -1;
            for (int bp = 0; bp < Geo.getBlockNum(); ++bp) {
              const auto& candidate = Geo.getBlock(bp);
              if (candidate.getMin()[0] == partner[1] &&
                  candidate.getMin()[1] == partner[2] &&
                  candidate.getMin()[2] == partner[3]) {
                partnerBlock = bp;
                break;
              }
            }
            if (partnerBlock >= 0) {
              auto& partnerField = field.getBlockField(partnerBlock);
              const int partnerNy = Geo.getBlock(partnerBlock).getNy();
              if (atFront) {
                for (int k = 0; k < nz; ++k) {
                  for (int i = 0; i < nx; ++i) {
                    blockField.get(k * projection[2] + i) =
                        partnerField.get(k * projection[2] + (partnerNy - 2) * projection[1] + i);
                  }
                }
              }
              if (atBack) {
                for (int k = 0; k < nz; ++k) {
                  for (int i = 0; i < nx; ++i) {
                    blockField.get(k * projection[2] + (ny - 1) * projection[1] + i) =
                        partnerField.get(k * projection[2] + projection[1] + i);
                  }
                }
              }
            }
          } else {
            const int tag = fieldIndex * 1000 + std::min(localEntry, partnerEntry) * 2 + 1;
            const int sourceJ = atFront ? 1 : ny - 2;
            const int targetJ = atFront ? 0 : ny - 1;
            std::vector<T> sendBuffer(plane);
            std::vector<T> receiveBuffer(plane);
            int index = 0;
            for (int k = 0; k < nz; ++k) {
              for (int i = 0; i < nx; ++i) {
                sendBuffer[index++] = blockField.get(k * projection[2] +
                                                    sourceJ * projection[1] + i);
              }
            }
            mpi().sendRecv(sendBuffer.data(), receiveBuffer.data(), plane,
                           partnerRank, partnerRank, tag);
            index = 0;
            for (int k = 0; k < nz; ++k) {
              for (int i = 0; i < nx; ++i) {
                blockField.get(k * projection[2] + targetJ * projection[1] + i) =
                    receiveBuffer[index++];
              }
            }
          }
        }
      }
    }
  };

  Printer::Print_BigBanner(std::string("Start Calculation..."));
  Timer timer;
  Timer outputTimer;

  while (timer() < MaxStep) {
    // ===== Phase 0: warm-started magnetic field solve =====
    MCC.ApplyInnerCellDynamics<MCSel>(timer(), FlagFM);
    CommunicateOMEGAPSI<T>(MFLattice);
    PinPsiWalls();

    T psiResidual = T{};
    for (int sub = 0; sub < PsiSolver_Iter; ++sub) {
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& blockLat = MFLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              MFCELL cell(k * projection[2] + j * projection[1] + i, blockLat);
              collision::MRTDiffusion<MFCELL, OMEGA_PSI<T>>::apply(cell);
            }
          }
        }
      }
      MF_Per.Apply();

      MFLattice.NormalFullCommunicate();
      MFLattice.Stream();
      MFLattice.NormalFullCommunicate();

      T localResidual = T{};
      auto& blockPsiField = MFLattice.getField<PSI<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& blockLat = MFLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockPsi = blockPsiField.getBlockField(b);
        const int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              const std::size_t id = k * projection[2] + j * projection[1] + i;
              MFCELL cell(id, blockLat);
              T psi = T{};
              for (unsigned int q = 0; q < MFLatSet::q; ++q) psi += cell[q];
              localResidual = std::max(localResidual,
                                       std::abs(psi - blockPsi.get(id)));
              blockPsi.get(id) = psi;
            }
          }
        }
      }
      CommunicatePSI<T>(MFLattice);
      PinPsiWalls();

#ifdef MPI_ENABLED
      T globalResidual = T{};
      mpi().reduce(localResidual, globalResidual, MPI_MAX);
      mpi().bCast(globalResidual, 0);
      psiResidual = globalResidual;
#else
      psiResidual = localResidual;
#endif
      if (psiResidual <= PsiSolver_Residual) break;
    }

    // Periodic ghost planes of scalar fields are not filled by MF_Per.Apply.
    SyncMFPeriodicGhosts(MFLattice.getField<PSI<T>>(), 0);

    // 0f: H = -grad(psi), then synchronize all magnetic output fields.
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& blockLat = MFLattice.getBlockLat(b);
      const auto& block = Geo.getBlock(b);
      const auto& projection = block.getProjection();
      const int overlap = block.getOverlap();
      for (int k = overlap; k < block.getNz() - overlap; ++k) {
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            MFCELL cell(k * projection[2] + j * projection[1] + i, blockLat);
            MFComputeH3D<MFCELL>::apply(cell);
          }
        }
      }
    }
    CommunicateAllMFFields3D<T>(MFLattice);
    SyncMFPeriodicGhosts(MFLattice.getField<HX<T>>(), 1);
    SyncMFPeriodicGhosts(MFLattice.getField<HY<T>>(), 2);
    SyncMFPeriodicGhosts(MFLattice.getField<HZ<T>>(), 3);
    SyncMFPeriodicGhosts(MFLattice.getField<HMAG<T>>(), 4);

    // ===== Phase A: force setup =====
    RoC.ApplyInnerCellDynamics<
        CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, RoT>>(timer(), FlagFM);
    NSLattice.getField<FORCE<T, 3>>().InitValue(Vector<T, 3>{0, 0, 0});
    STC.ApplyInnerCellDynamics<
        CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, STT>>(timer(), FlagFM);

    // Fm = (mu0*chi/2)*grad(|H|^2) = chi*|H|*grad(|H|), mu0=1.
    if (H0 > T{0} && chi_h > T{0}) {
      const T wallBand = T{16} * Cell_Len;
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& pfBlock = PFLattice.getBlockLat(b);
        auto& mfBlock = MFLattice.getBlockLat(b);
        auto& nsBlock = NSLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T voxelSize = block.getVoxelSize();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          const T z = minZ + T(k - overlap) * voxelSize;
          if (z <= wallBand || z >= H_global - wallBand) continue;
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              const std::size_t id = k * projection[2] + j * projection[1] + i;
              PFCELL pfCell(id, pfBlock);
              MFCELL mfCell(id, mfBlock);
              NSCELL nsCell(id, nsBlock);
              MFMagneticForce3D<PFCELL, MFCELL, NSCELL>::apply(
                  pfCell, mfCell, nsCell);
            }
          }
        }
      }
    }
    GrC.ApplyInnerCellDynamics<
        CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, GrT>>(timer(), FlagFM);
    PrC.ApplyInnerCellDynamics<
        CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, PrT>>(timer(), FlagFM);
    NSLattice.getField<FORCE<T, 3>>().Communicate();

    // ===== Phase B-C: phase-field and NS collision =====
    PFLattice.template ApplyInnerCellDynamics<PFSel>(FlagFM);
    PF_Per.Apply();
    PFLattice.NormalFullCommunicate();

    ViC.ApplyInnerCellDynamics<
        CoupledTaskSelector<std::uint8_t, PFCELL, NSCELL, ViT>>(timer(), FlagFM);
    NSLattice.template ApplyInnerCellDynamics<NSSel>(FlagFM);
    NS_Per.Apply();
    NSLattice.NormalFullCommunicate();

    // ===== Phase D: streaming and z no-slip boundary handling =====
    PF_BB.Apply(timer());
    {
      const T domainHeight = T(Nk) * Cell_Len;
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& blockLat = PFLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int nz = block.getNz();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              PFCELL ghost(j * projection[1] + i, blockLat);
              PFCELL wall(overlap * projection[2] + j * projection[1] + i,
                          blockLat);
              for (unsigned int q = 0; q < LatSet::q; ++q) ghost[q] = wall[q];
            }
          }
        }
        if (maxZ > domainHeight - Cell_Len * T{1.5}) {
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              PFCELL ghost((nz - 1) * projection[2] + j * projection[1] + i,
                           blockLat);
              PFCELL wall((nz - 1 - overlap) * projection[2] +
                              j * projection[1] + i,
                          blockLat);
              for (unsigned int q = 0; q < LatSet::q; ++q) ghost[q] = wall[q];
            }
          }
        }
      }
    }
    PFLattice.Stream();
    PFLattice.NormalFullCommunicate();

    NS_BB.Apply(timer());
    {
      const T domainHeight = T(Nk) * Cell_Len;
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& blockLat = NSLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int nz = block.getNz();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              NSCELL ghost(j * projection[1] + i, blockLat);
              NSCELL wall(overlap * projection[2] + j * projection[1] + i,
                          blockLat);
              for (unsigned int q = 0; q < LatSet::q; ++q) ghost[q] = wall[q];
            }
          }
        }
        if (maxZ > domainHeight - Cell_Len * T{1.5}) {
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              NSCELL ghost((nz - 1) * projection[2] + j * projection[1] + i,
                           blockLat);
              NSCELL wall((nz - 1 - overlap) * projection[2] +
                              j * projection[1] + i,
                          blockLat);
              for (unsigned int q = 0; q < LatSet::q; ++q) ghost[q] = wall[q];
            }
          }
        }
      }
    }
    NSLattice.Stream();
    NSLattice.NormalFullCommunicate();

    // ===== Phase E: phase-field macro and wall values =====
    {
      auto& field = PFLattice.getField<PHI<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        auto& blockLat = PFLattice.getBlockLat(b);
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockPhi = field.getBlockField(b);
        const int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              const std::size_t id = k * projection[2] + j * projection[1] + i;
              PFCELL cell(id, blockLat);
              T phi = T{};
              for (unsigned int q = 0; q < LatSet::q; ++q) phi += cell[q];
              phi = std::max(T{0}, std::min(T{1}, phi));
              blockPhi.get(id) = phi;
            }
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // The bottom wall is ferrofluid and the top wall is organic fluid.
    {
      auto& field = PFLattice.getField<PHI<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockPhi = field.getBlockField(b);
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int nz = block.getNz();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(overlap * projection[2] + j * projection[1] + i) = T{1};
            }
          }
        }
        if (maxZ > H_global - Cell_Len * T{1.5}) {
          const int top = nz - 1 - overlap;
          for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(top * projection[2] + j * projection[1] + i) = T{0};
            }
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();

    // Keep the phase ghost layers consistent with their adjacent wall phase.
    {
      auto& field = PFLattice.getField<PHI<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockPhi = field.getBlockField(b);
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int nz = block.getNz();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int k = 0; k < overlap; ++k) {
            for (int j = 0; j < ny; ++j) {
              for (int i = 0; i < nx; ++i) {
                blockPhi.get(k * projection[2] + j * projection[1] + i) = T{1};
              }
            }
          }
        }
        if (maxZ > H_global - Cell_Len * T{1.5}) {
          for (int k = nz - overlap; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
              for (int i = 0; i < nx; ++i) {
                blockPhi.get(k * projection[2] + j * projection[1] + i) = T{0};
              }
            }
          }
        }
      }
    }

    // Gradients and chemical potential use the updated wall phase values.
    PFLattice.template ApplyInnerCellDynamics<PFSelN>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelL>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<PFSelC>(FlagFM);
    PFLattice.getField<NORMAL<T, 3>>().Communicate();
    PFLattice.getField<GRAD<T, 3>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    {
      auto& field = PFLattice.getField<GRAD<T, 3>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockGrad = field.getBlockField(b);
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              blockGrad.get(overlap * projection[2] + j * projection[1] + i)[2] =
                  blockGrad.get((overlap + 1) * projection[2] +
                                j * projection[1] + i)[2];
            }
          }
        }
        if (maxZ > H_global - Cell_Len * T{1.5}) {
          const int top = block.getNz() - 1 - overlap;
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              blockGrad.get(top * projection[2] + j * projection[1] + i)[2] =
                  blockGrad.get((top - 1) * projection[2] +
                                j * projection[1] + i)[2];
            }
          }
        }
      }
    }
    {
      auto& field = PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      for (int b = 0; b < Geo.getBlockNum(); ++b) {
        const auto& block = Geo.getBlock(b);
        const auto& projection = block.getProjection();
        auto& blockChem = field.getBlockField(b);
        const int nx = block.getNx();
        const int ny = block.getNy();
        const int overlap = block.getOverlap();
        const T minZ = block.getMin()[2];
        const T maxZ = block.getMax()[2];
        if (minZ < Cell_Len * T{1.5}) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              blockChem.get(overlap * projection[2] + j * projection[1] + i) =
                  (T{4} * blockChem.get((overlap + 1) * projection[2] +
                                        j * projection[1] + i) -
                   blockChem.get((overlap + 2) * projection[2] +
                                 j * projection[1] + i)) / T{3};
            }
          }
        }
        if (maxZ > H_global - Cell_Len * T{1.5}) {
          const int top = block.getNz() - 1 - overlap;
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              blockChem.get(top * projection[2] + j * projection[1] + i) =
                  (T{4} * blockChem.get((top - 1) * projection[2] +
                                        j * projection[1] + i) -
                   blockChem.get((top - 2) * projection[2] +
                                 j * projection[1] + i)) / T{3};
            }
          }
        }
      }
    }
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // NS macroscopic pressure, velocity, and density fields.
    for (int b = 0; b < Geo.getBlockNum(); ++b) {
      auto& blockLat = NSLattice.getBlockLat(b);
      const auto& block = Geo.getBlock(b);
      const auto& projection = block.getProjection();
      auto& rhoField = NSLattice.getField<DENSITY<T>>().getBlockField(b);
      auto& pressureField = NSLattice.getField<PRESSURE<T>>().getBlockField(b);
      auto& velocityField = NSLattice.getField<VELOCITY<T, 3>>().getBlockField(b);
      auto& forceField = NSLattice.getField<FORCE<T, 3>>().getBlockField(b);
      const int overlap = block.getOverlap();
      for (int k = overlap; k < block.getNz() - overlap; ++k) {
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            const std::size_t id = k * projection[2] + j * projection[1] + i;
            NSCELL cell(id, blockLat);
            T pressure = T{};
            Vector<T, 3> momentum{0, 0, 0};
            for (unsigned int q = 0; q < LatSet::q; ++q) {
              pressure += cell[q];
              momentum += latset::c<LatSet>(q) * cell[q];
            }
            const T rho = rhoField.get(id);
            const auto force = forceField.get(id);
            pressureField.get(id) = pressure;
            velocityField.get(id) =
                (momentum + T{0.5} * force / rho);
          }
        }
      }
    }
    NSLattice.getField<VELOCITY<T, 3>>().Communicate();
    NSLattice.getField<PRESSURE<T>>().Communicate();
    NSLattice.getField<DENSITY<T>>().Communicate();
    NSLattice.getField<OMEGA<T>>().Communicate();

    ++timer;
    ++outputTimer;
    if (timer() % OutputStep == 0) {
      outputTimer.Print_InnerLoopPerformance(Geo.getTotalCellNum(), OutputStep);
      Printer::Endl();
      MW.WriteBinary(timer());
    }
  }

  Printer::Print_BigBanner(std::string("Calculation Complete!"));
  timer.Print_MainLoopPerformance(Geo.getTotalCellNum());
  return 0;
}
