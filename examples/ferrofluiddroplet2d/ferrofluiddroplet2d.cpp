// ferrofluiddroplet2d.cpp
// 2D ferrofluid droplet deformation under uniform magnetic field
// Hu & Li (2018) Phys. Rev. E 98, 033301 — Case B
// Three coupled lattices: NS(D2Q9) + PF(D2Q9) + MF(D2Q5)
// φ=1: ferrofluid droplet (heavy),  φ=0: organic liquid (light)
// gravity=0 — static deformation, magnetic force vs surface tension

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

// droplet
T Droplet_Radius;
Vector<T, LatSet::d> Droplet_Center;

// phase field
T Interface_Width;
T Mobility;
T Tau_phi;
T Omega_phi;
T Kappa;
T Beta;

// two-phase — Hu & Li (2018) Case B lattice units
T rho_l;
T rho_h;
T eta_l;
T eta_h;
T sigma;
T gravity;
T Tau_ns;

// magnetic field — Hu & Li (2018)
T Mu_l;
T Mu_h;
T Chi_l;
T Chi_h;
T H0;
T Epsilon;
T Tau_mag_ref;

// simulation
int MaxStep;
int OutputStep;

std::string work_dir;

void readParam() {
  iniReader param_reader("ferrofluiddroplet2d.ini");
  work_dir = param_reader.getValue<std::string>("workdir", "workdir_");
  Thread_Num = param_reader.getValue<int>("parallel", "thread_num");

  Ni = param_reader.getValue<int>("Mesh", "Ni");
  Nj = param_reader.getValue<int>("Mesh", "Nj");
  Cell_Len = param_reader.getValue<T>("Mesh", "Cell_Len");
  BlockCellLen = param_reader.getValue<int>("Mesh", "BlockCellLen");

  Droplet_Radius = param_reader.getValue<T>("Droplet", "Radius");
  Droplet_Center[0] = param_reader.getValue<T>("Droplet", "CenterX");
  Droplet_Center[1] = param_reader.getValue<T>("Droplet", "CenterY");

  Interface_Width = param_reader.getValue<T>("Phase_Field", "Interface_Width");
  Mobility = param_reader.getValue<T>("Phase_Field", "Mobility");

  rho_l = param_reader.getValue<T>("Two_Phase", "rho_l");
  rho_h = param_reader.getValue<T>("Two_Phase", "rho_h");

  // Magnetic parameters
  Mu_l = param_reader.getValue<T>("Magnetic", "Mu_l");
  Mu_h = param_reader.getValue<T>("Magnetic", "Mu_h");
  Chi_l = param_reader.getValue<T>("Magnetic", "Chi_l");
  Chi_h = param_reader.getValue<T>("Magnetic", "Chi_h");
  H0 = param_reader.getValue<T>("Magnetic", "H0");

  MaxStep = param_reader.getValue<int>("Simulation_Settings", "TotalStep");
  OutputStep = param_reader.getValue<int>("Simulation_Settings", "OutputStep");

  // ---- Hu & Li (2018) Case B hardcoded lattice parameters ----
  // ρ_l=1.0, ρ_h=1.975, η_l=0.0025, η_h=0.05, σ=0.00191875
  // Reference: ρ_ref=800 kg/m³, L_ref=0.8 mm, U_ref=5 m/s
  eta_l = T(0.0025);
  eta_h = T(0.05);
  sigma = T(0.00191875);
  gravity = T(0);  // NO gravity — static deformation

  Beta = T(12.0) * sigma / Interface_Width;
  Kappa = T(3.0) * Interface_Width * sigma * T(0.5);
  Tau_phi = T(3.0) * Mobility + T(0.5);
  Omega_phi = T(1.0) / Tau_phi;
  Tau_ns = T(0.5) + eta_h / rho_h / LatSet::cs2;

  // Magnetic diffusion: τ(x) = 0.5 + ε·μ(x)/c_s²  (Eq.43)
  Tau_mag_ref = T{4.0};
  T mu_avg = (Mu_l + Mu_h) / T{2};
  Epsilon = (Tau_mag_ref - T{0.5}) * MFLatSet::cs2 / mu_avg;

  MPI_RANK(0) {
    std::cout << "-----Ferrofluid Droplet Deformation (Hu2018 Case B)-----\n";
    std::cout << "[Mesh]: " << Ni << "x" << Nj << "  BlockCellLen=" << BlockCellLen << "\n";
    std::cout << "[Droplet]: R=" << Droplet_Radius
              << "  Center=(" << Droplet_Center[0] << "," << Droplet_Center[1] << ")\n";
    std::cout << "[Fluid]: rho_l=" << rho_l << " rho_h=" << rho_h
              << " eta_l=" << eta_l << " eta_h=" << eta_h
              << " sigma=" << sigma << " gravity=" << gravity << "\n";
    std::cout << "[Phase]: W=" << Interface_Width << " M=" << Mobility
              << " tau_phi=" << Tau_phi << " beta=" << Beta << " kappa=" << Kappa << "\n";
    std::cout << "[Magnetic]: mu_l=" << Mu_l << " mu_h=" << Mu_h
              << " chi_l=" << Chi_l << " chi_h=" << Chi_h
              << " H0=" << H0 << " eps=" << Epsilon << " tau_ref=" << Tau_mag_ref << "\n";
    std::cout << "[Simulation]: MaxStep=" << MaxStep
              << "  OutputStep=" << OutputStep << "\n";
    T Bom_est = Chi_h * H0 * H0 * Droplet_Radius / sigma;
    std::cout << "[Estimate]: Bom ≈ " << Bom_est << "\n";
#ifdef _OPENMP
    std::cout << "[Parallel]: " << Thread_Num << " threads\n";
#endif
#ifdef MPI_ENABLED
    std::cout << "[Parallel]: " << mpi().getSize() << " MPI processes\n";
#endif
    std::cout << "---------------------------------------------------------\n";
  }
}

// Compute droplet aspect ratio b/a from φ=0.5 contour
// Returns {b, a, b/a} — b=semi-major (vertical), a=semi-minor (horizontal)
template <typename PFLAT>
std::array<T, 3> computeAspectRatio(PFLAT& PFLattice, const BlockGeometry2D<T>& Geo) {
  T x_min = T(1e10), x_max = T(-1e10);
  T y_min = T(1e10), y_max = T(-1e10);
  bool found = false;

  auto& phiField = PFLattice.template getField<PHI<T>>();
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    T voxelSize = block.getVoxelSize();
    int overlap = block.getOverlap();
    for (int j = overlap; j < block.getNy() - overlap; ++j) {
      for (int i = overlap; i < block.getNx() - overlap; ++i) {
        // Check neighbors for φ=0.5 crossing
        if (i > overlap && i < block.getNx() - overlap - 1 &&
            j > overlap && j < block.getNy() - overlap - 1) {
          T phi_l = blockPhi.get(j * proj[1] + (i-1));
          T phi_r = blockPhi.get(j * proj[1] + (i+1));
          T phi_b = blockPhi.get((j-1) * proj[1] + i);
          T phi_t = blockPhi.get((j+1) * proj[1] + i);
          bool crosses_x = (phi_l - T{0.5}) * (phi_r - T{0.5}) < T{0};
          bool crosses_y = (phi_b - T{0.5}) * (phi_t - T{0.5}) < T{0};
          if (crosses_x || crosses_y) {
            T x = block.getMin()[0] + static_cast<T>(i) * voxelSize;
            T y = block.getMin()[1] + static_cast<T>(j) * voxelSize;
            if (crosses_y) {
              if (x < x_min) x_min = x;
              if (x > x_max) x_max = x;
            }
            if (crosses_x) {
              if (y < y_min) y_min = y;
              if (y > y_max) y_max = y;
            }
            found = true;
          }
        }
      }
    }
  }

  // MPI reduction: each rank only searches its own blocks,
  // so aggregate min/max across all ranks to get global contour.
#ifdef MPI_ENABLED
  {
    int found_int = found ? 1 : 0;
    int found_global = 0;
    MPI_Allreduce(&found_int, &found_global, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    found = (found_global != 0);

    T buf[4] = {-x_min, x_max, -y_min, y_max};  // negate min → MPI_MAX = global min
    T rbuf[4];
    MPI_Allreduce(buf, rbuf, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    x_min = -rbuf[0];
    x_max =  rbuf[1];
    y_min = -rbuf[2];
    y_max =  rbuf[3];
  }
#endif

  if (!found) return {T{0}, T{0}, T{0}};

  T a = (x_max - x_min) / T{2};  // semi-minor (horizontal)
  T b = (y_max - y_min) / T{2};  // semi-major (vertical)
  T ba = (a > T{0}) ? b / a : T{0};
  return {b, a, ba};
}

int main(int argc, char* argv[]) {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);
  constexpr std::uint8_t BouncebackFlag = std::uint8_t(4);
  constexpr std::uint8_t PeriodicFlag = std::uint8_t(8);

  mpi().init(&argc, &argv);
  MPI_DEBUG_WAIT

  Printer::Print_BigBanner(std::string("Initializing Ferrofluid Droplet Deformation..."));

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
  GeoHelper.CreateBlocks(2, 2);
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
  vtmo::vtmWriter<T, 2> GeoWriter("GeoFlag_Droplet", Geo, 1);
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

  // ---- Init PHI: φ=1 inside droplet (ferrofluid), φ=0 outside (organic) ----
  T R_phys = Droplet_Radius * Cell_Len;
  T xc_phys = Droplet_Center[0] * Cell_Len;
  T yc_phys = Droplet_Center[1] * Cell_Len;
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
        // NOTE: φ=1 INSIDE droplet (ferrofluid), opposite of bubble case
        // tanh profile: φ → 1 when dist < R, φ → 0 when dist > R
        T phi = T(0.5) - T(0.5) * std::tanh(T(2.0) * (dist - R_phys) / W_phys);
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

  // Init OMEGA on MF lattice
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

  vtmo::vtmWriter<T, 2> MainWriter("ferrofluiddroplet2d", Geo);
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
  MF_Periodic.Apply();

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

  // Initial aspect ratio + output
  {
    auto ar = computeAspectRatio(PFLattice, Geo);
    if (mpi().getRank() == 0) {
      std::cout << "[Step " << MainLoopTimer() << "] b=" << ar[0]
                << " a=" << ar[1] << " b/a=" << ar[2] << "\n";
    }
  }
  MainWriter.WriteBinary(MainLoopTimer());

  Printer::Print_BigBanner(std::string("Start Calculation..."));

  // ================== Main Loop ==================
  while (MainLoopTimer() < MaxStep) {

    // ---- Phase A: Force accumulation ----
    RhoOmegaCoupling.ApplyInnerCellDynamics<RhoOmegaTaskSelector>(MainLoopTimer(), FlagFM);
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, 2>{T{0}, T{0}});

    // A1: F_s = λ·∇φ
    STCoupling.ApplyInnerCellDynamics<STForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A2: F_b = -ρ·g  (g=0 for Case B — static deformation, no gravity)
    GravCoupling.ApplyInnerCellDynamics<GravForceTaskSelector>(MainLoopTimer(), FlagFM);

    // A3: F_p = -(p/3)·Δρ·∇φ
    PreForceCoupling.ApplyInnerCellDynamics<PreForceTaskSelector>(MainLoopTimer(), FlagFM);

    // ---- Phase A_mag: Magnetic force ----
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

    // MF collision
    MFLattice.template ApplyInnerCellDynamics<MFCollisionSelector>(FlagFM);
    MF_Periodic.Apply();

    // MF bounceback + stream
    MF_BB.Apply(MainLoopTimer());
    MFLattice.Stream();
    MFLattice.NormalFullCommunicate();

    // MF Y ghost POPs = wall POPs
    {
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockLat = MFLattice.getBlockLat(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        if (minY < Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = 0 * proj[1] + i;
            std::size_t id_wall = overlap * proj[1] + i;
            MFCELL ghost_cell(id_ghost, blockLat);
            MFCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id_ghost = (ny - 1) * proj[1] + i;
            std::size_t id_wall = (ny - 1 - overlap) * proj[1] + i;
            MFCELL ghost_cell(id_ghost, blockLat);
            MFCELL wall_cell(id_wall, blockLat);
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              ghost_cell[k] = wall_cell[k];
            }
          }
        }
      }
    }

    // MF: ψ = Σ h_α
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
            T psi = T{0};
            for (unsigned int k = 0; k < MFLatSet::q; ++k) {
              psi += cell[k];
            }
            blockPsi.get(id) = psi;
          }
        }
      }
    }
    MFLattice.getField<ff::PSI<T>>().Communicate();

    // MF: ψ ghost for H₀ BC at top/bottom
    {
      auto& psiField = MFLattice.getField<ff::PSI<T>>();
      T H_global = T(Nj) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPsi = psiField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int overlap = block.getOverlap();
        T minY = block.getMin()[1];
        T maxY = block.getMax()[1];
        T dy = block.getVoxelSize();
        if (minY < Cell_Len * T(1.5)) {
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_ghost = j * proj[1] + i;
              std::size_t id_inner = overlap * proj[1] + i;
              int ghost_depth = overlap - j;
              blockPsi.get(id_ghost) = blockPsi.get(id_inner) + H0 * dy * T(ghost_depth);
            }
          }
        }
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              std::size_t id_ghost = j * proj[1] + i;
              std::size_t id_inner = (ny - 1 - overlap) * proj[1] + i;
              int ghost_depth = j - (ny - 1 - overlap);
              blockPsi.get(id_ghost) = blockPsi.get(id_inner) - H0 * dy * T(ghost_depth);
            }
          }
        }
      }
    }
    MFLattice.getField<ff::PSI<T>>().Communicate();

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
    PF_BB.Apply(MainLoopTimer());
    // PF Y ghost POP fix
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

    // φ=0 at walls (organic liquid wets walls)
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
        // Bottom wall: φ = 0
        if (minY < Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            blockPhi.get(overlap * proj[1] + i) = T{0};
          }
          for (int j = 0; j < overlap; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{0};
            }
          }
        }
        // Top wall: φ = 0
        if (maxY > H_global - Cell_Len * T(1.5)) {
          for (int i = overlap; i < nx - overlap; ++i) {
            blockPhi.get((ny - 1 - overlap) * proj[1] + i) = T{0};
          }
          for (int j = ny - overlap; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
              blockPhi.get(j * proj[1] + i) = T{0};
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
        auto& presField = NSLattice.getField<PRESSURE<T>>();
        auto& velField = NSLattice.getField<VELOCITY<T, 2>>();
        int overlap = block.getOverlap();
        for (int j = overlap; j < block.getNy() - overlap; ++j) {
          for (int i = overlap; i < block.getNx() - overlap; ++i) {
            std::size_t id = j * proj[1] + i;
            NSCELL cell(id, blockLat);
            T p_new = T{0};
            Vector<T, 2> u_new{T{0}, T{0}};
            for (unsigned int k = 0; k < LatSet::q; ++k) {
              p_new += cell[k];
              u_new[0] += latset::c<LatSet>(k)[0] * cell[k];
              u_new[1] += latset::c<LatSet>(k)[1] * cell[k];
            }
            presField.getBlockField(blockid).get(id) = p_new;
            velField.getBlockField(blockid).get(id) = u_new;
          }
        }
      }
    }

    // ---- Output ----
    ++MainLoopTimer;

    if (MainLoopTimer() % OutputStep == 0) {
      auto ar = computeAspectRatio(PFLattice, Geo);
      if (mpi().getRank() == 0) {
        std::cout << "[Step " << MainLoopTimer() << "] b=" << ar[0]
                  << " a=" << ar[1] << " b/a=" << ar[2] << "\n";
      }
      MainWriter.WriteBinary(MainLoopTimer());
    }
  }

  Printer::Print_BigBanner(std::string("Simulation Complete."));
  return 0;
}
