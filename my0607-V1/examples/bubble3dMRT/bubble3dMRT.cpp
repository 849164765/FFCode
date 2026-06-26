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

// Set phi=1 on Z-direction no-slip wall cells (top/bottom), including ghost
// layers.  X/Y boundaries are periodic, so their ghost cells are filled by
// FixedPeriodicBoundaryManager, NOT set to phi=1 here.
// Called once during init and again each time step after phi is reconstructed
// from the PF populations.
template <typename T, typename LatSet>
void applyWallPhi(BlockGeometry3D<T>& Geo,
                  BlockFieldManager<GenericField<GenericArray<T>, PHIBase>, T, LatSet::d>& phiField,
                  T Cell_Len) {
  T Lz_global = T(Geo.getNz()) * Cell_Len;
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    auto& blockPhi = phiField.getBlockField(blockid);
    int nx = block.getNx();
    int ny = block.getNy();
    int nz = block.getNz();
    int overlap = block.getOverlap();
    T minZ = block.getMin()[2];
    T maxZ = block.getMax()[2];
    // z = 0 face
    if (minZ < Cell_Len * T(1.5)) {
      for (int k = 0; k <= overlap; ++k) {
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id = k * proj[2] + j * proj[1] + i;
            blockPhi.get(id) = T{1};
          }
        }
      }
    }
    // z = Nz-1 face
    if (maxZ > Lz_global - Cell_Len * T(1.5)) {
      for (int k = nz - 1 - overlap; k <= nz - 1; ++k) {
        for (int j = 0; j < ny; ++j) {
          for (int i = 0; i < nx; ++i) {
            std::size_t id = k * proj[2] + j * proj[1] + i;
            blockPhi.get(id) = T{1};
          }
        }
      }
    }
  }
}

// Reset Z-direction domain-boundary ghost-layer POP to correct equilibrium
// values after stream.  X/Y ghost layers are filled by
// FixedPeriodicBoundaryManager (periodic), so only Z (top/bottom walls) need
// this reset.
// Stream (PULL scheme) leaves domain-boundary ghost cells reading out-of-domain
// memory, and NormalCommunicate does not fill them (no neighbour).
//
// PF POP  ← w[i] * 1          (phi_wall = 1, u = 0)
// NS POP  ← 0                 (IncompressibleSecondOrder p = 0, u = 0)
template <typename T, typename LatSet, typename PFLATMAN, typename NSLATMAN>
void applyWallPop(BlockGeometry3D<T>& Geo, PFLATMAN& PFLattice, NSLATMAN& NSLattice, T Cell_Len) {
  T Lz_global = T(Geo.getNz()) * Cell_Len;
  for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
    const auto& block = Geo.getBlock(blockid);
    const auto& proj = block.getProjection();
    int nx = block.getNx();
    int ny = block.getNy();
    int nz = block.getNz();
    int overlap = block.getOverlap();
    T minZ = block.getMin()[2];
    T maxZ = block.getMax()[2];

    auto& pfPopField = PFLattice.getBlockLat(blockid).template getField<POP<T, LatSet::q>>();
    auto& nsPopField = NSLattice.getBlockLat(blockid).template getField<POP<T, LatSet::q>>();
    auto& nsPresField = NSLattice.template getField<PRESSURE<T>>().getBlockField(blockid);
    auto& nsVelField = NSLattice.template getField<VELOCITY<T, LatSet::d>>().getBlockField(blockid);
    auto& nsDensField = NSLattice.template getField<DENSITY<T>>().getBlockField(blockid);

    auto setWallPop = [&](std::size_t id) {
      for (unsigned int i_pop = 0; i_pop < LatSet::q; ++i_pop) {
        pfPopField.getField(i_pop)[id] = latset::w<LatSet>(i_pop) * T{1};
        nsPopField.getField(i_pop)[id] = T{0};
      }
      nsPresField.get(id) = T{0};
      nsVelField.get(id) = Vector<T, LatSet::d>{T{0}};
      nsDensField.get(id) = rho_h;
    };

    // z = 0 face
    if (minZ < Cell_Len * T(1.5)) {
      for (int k = 0; k < overlap; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            setWallPop(k * proj[2] + j * proj[1] + i);
    }
    // z = Nz-1 face
    if (maxZ > Lz_global - Cell_Len * T(1.5)) {
      for (int k = nz - overlap; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
            setWallPop(k * proj[2] + j * proj[1] + i);
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
  AABB<T, LatSet::d> domain(Vector<T, LatSet::d>(T(0), T(0), T(0)),
                            Vector<T, LatSet::d>(T(Ni * Cell_Len), T(Nj * Cell_Len), T(Nz * Cell_Len)));
  // X-periodic ghost layers (left/right)
  AABB<T, LatSet::d> leftX(Vector<T, LatSet::d>(T(-Cell_Len), T(0), T(0)),
                           Vector<T, LatSet::d>(T(0), T(Nj * Cell_Len), T(Nz * Cell_Len)));
  AABB<T, LatSet::d> rightX(Vector<T, LatSet::d>(T(Ni * Cell_Len), T(0), T(0)),
                            Vector<T, LatSet::d>(T((Ni + 1) * Cell_Len), T(Nj * Cell_Len), T(Nz * Cell_Len)));
  // Y-periodic ghost layers (front/back)
  AABB<T, LatSet::d> leftY(Vector<T, LatSet::d>(T(0), T(-Cell_Len), T(0)),
                           Vector<T, LatSet::d>(T(Ni * Cell_Len), T(0), T(Nz * Cell_Len)));
  AABB<T, LatSet::d> rightY(Vector<T, LatSet::d>(T(0), T(Nj * Cell_Len), T(0)),
                            Vector<T, LatSet::d>(T(Ni * Cell_Len), T((Nj + 1) * Cell_Len), T(Nz * Cell_Len)));

  BlockGeometryHelper3D<T> GeoHelper(Ni, Nj, Nz, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks();
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());

  BlockGeometry3D<T> Geo(GeoHelper);

  // ------------------ define flag field ------------------
  BlockFieldManager<FLAG, T, LatSet::d> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain,
                 [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  // Z-direction remains no-slip (bounceback on top/bottom walls)
  // SetupBoundary marks domain-boundary cells as BouncebackFlag AND calls
  // AllNormalCommunicate, which syncs flags across inter-block ghost layers.
  // PeriodicFlag must be set AFTER SetupBoundary so AllNormalCommunicate
  // doesn't overwrite it.
  FlagFM.template SetupBoundary<LatSet>(domain, BouncebackFlag);

  // Fix: SetupBoundary marks ALL domain-boundary cells as BouncebackFlag,
  // including X/Y periodic boundaries.  Reset X/Y boundary INTERNAL cells
  // (excluding Z walls) back to BulkFlag so they participate in collision.
  // Without this, X/Y boundary cells are skipped by ApplyInnerCellDynamics,
  // and their POP retains constructor defaults (w[i]), producing p=1/6.
  {
    T z_inner_min = Cell_Len;
    T z_inner_max = T(Nz - 1) * Cell_Len;
    AABB<T, LatSet::d> xMinBulk(
      Vector<T, LatSet::d>(T(0), T(0), z_inner_min),
      Vector<T, LatSet::d>(Cell_Len, T(Nj) * Cell_Len, z_inner_max));
    AABB<T, LatSet::d> xMaxBulk(
      Vector<T, LatSet::d>(T(Ni - 1) * Cell_Len, T(0), z_inner_min),
      Vector<T, LatSet::d>(T(Ni) * Cell_Len, T(Nj) * Cell_Len, z_inner_max));
    AABB<T, LatSet::d> yMinBulk(
      Vector<T, LatSet::d>(T(0), T(0), z_inner_min),
      Vector<T, LatSet::d>(T(Ni) * Cell_Len, Cell_Len, z_inner_max));
    AABB<T, LatSet::d> yMaxBulk(
      Vector<T, LatSet::d>(T(0), T(Nj - 1) * Cell_Len, z_inner_min),
      Vector<T, LatSet::d>(T(Ni) * Cell_Len, T(Nj) * Cell_Len, z_inner_max));
    FlagFM.forEach(xMinBulk, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
    FlagFM.forEach(xMaxBulk, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
    FlagFM.forEach(yMinBulk, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
    FlagFM.forEach(yMaxBulk, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  }

  // Mark X/Y periodic ghost cells AFTER SetupBoundary (which calls
  // AllNormalCommunicate).  This prevents AllNormalCommunicate from
  // overwriting PeriodicFlag with the neighbour block's internal flag.
  // BBLikeFixedBlockBdManager only bounce-backs directions where the
  // neighbour is VoidFlag; since periodic ghost cells are PeriodicFlag
  // (not VoidFlag), X/Y directions are exempt.
  FlagFM.forEach(leftX, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(rightX, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(leftY, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });
  FlagFM.forEach(rightY, [&](FLAG& field, std::size_t id) { field.SetField(id, PeriodicFlag); });

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
    int overlap = block.getOverlap();
    for (int k = overlap; k < block.getNz() - overlap; ++k) {
      for (int j = overlap; j < block.getNy() - overlap; ++j) {
        for (int i = overlap; i < block.getNx() - overlap; ++i) {
          std::size_t id = k * proj[2] + j * proj[1] + i;
          T x = block.getMinCenter()[0] + static_cast<T>(i) * voxelSize;
          T y = block.getMinCenter()[1] + static_cast<T>(j) * voxelSize;
          T z = block.getMinCenter()[2] + static_cast<T>(k) * voxelSize;
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
    int overlap = block.getOverlap();
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
    int overlap = block.getOverlap();
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

  // X/Y periodic BCs (Z remains no-slip walls).
  // Lets pressure waves escape in X/Y, preventing the standing-wave feedback
  // loop that caused p_range explosion in the all-bounceback configuration.
  using LM_NS = BlockLatticeManager<T, LatSet, NSFIELDS>;
  using LM_PF = BlockLatticeManager<T, LatSet, PFFIELDPACK>;
  using FM = BlockFieldManager<FLAG, T, LatSet::d>;

  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_PeriodicX("NS_PeriodicX", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_PeriodicX.Setup(leftX, NbrDirection::XN, rightX, NbrDirection::XP);
  NS_PeriodicX.Setup(rightX, NbrDirection::XP, leftX, NbrDirection::XN);

  FixedPeriodicBoundaryManager<LM_NS, FM>
      NS_PeriodicY("NS_PeriodicY", NSLattice, FlagFM, PeriodicFlag, VoidFlag);
  NS_PeriodicY.Setup(leftY, NbrDirection::YN, rightY, NbrDirection::YP);
  NS_PeriodicY.Setup(rightY, NbrDirection::YP, leftY, NbrDirection::YN);

  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_PeriodicX("PF_PeriodicX", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_PeriodicX.Setup(leftX, NbrDirection::XN, rightX, NbrDirection::XP);
  PF_PeriodicX.Setup(rightX, NbrDirection::XP, leftX, NbrDirection::XN);

  FixedPeriodicBoundaryManager<LM_PF, FM>
      PF_PeriodicY("PF_PeriodicY", PFLattice, FlagFM, PeriodicFlag, VoidFlag);
  PF_PeriodicY.Setup(leftY, NbrDirection::YN, rightY, NbrDirection::YP);
  PF_PeriodicY.Setup(rightY, NbrDirection::YP, leftY, NbrDirection::YN);

  // Set up MPI communication for cross-process periodic pairs.
  // X periodic is same-block (X not split), but Y periodic is cross-block
  // (Y split into 2 blocks on different processes).  Without SetupMPI,
  // Y periodic ghost cells retain constructor defaults (w[i]), corrupting
  // the first stream and producing p=1/6 at Y-boundary cells.
  NS_PeriodicX.SetupMPI(GeoHelper);
  NS_PeriodicY.SetupMPI(GeoHelper);
  PF_PeriodicX.SetupMPI(GeoHelper);
  PF_PeriodicY.SetupMPI(GeoHelper);

  // Diagnostic: check ghost cell flags on each rank
  {
    int myRank = mpi().getRank();
    for (int bi = 0; bi < Geo.getBlockNum(); ++bi) {
      const auto& block = Geo.getBlock(bi);
      const auto& proj = block.getProjection();
      int nx = block.getNx();
      int ny = block.getNy();
      int nz = block.getNz();
      int overlap = block.getOverlap();
      auto& flagField = FlagFM.getBlockField(bi);
      // Check X-max ghost (i=nx-1) and Y-min ghost (j=0) at k=overlap
      std::size_t id_xmax = (overlap) * proj[2] + (overlap) * proj[1] + (nx - 1);
      std::size_t id_ymin = (overlap) * proj[2] + 0 * proj[1] + overlap;
      std::size_t id_xmin = (overlap) * proj[2] + (overlap) * proj[1] + 0;
      std::size_t id_ymax = (overlap) * proj[2] + (ny - 1) * proj[1] + overlap;
      std::cout << "[FLAG_DIAG] Rank=" << myRank << " block=" << bi
                << " nx=" << nx << " ny=" << ny << " nz=" << nz
                << " X-min-ghost(flag=" << int(flagField.get(0 * proj[2] + overlap * proj[1] + overlap)) << ")"
                << " X-max-ghost(flag=" << int(flagField.get(id_xmax)) << ")"
                << " Y-min-ghost(flag=" << int(flagField.get(id_ymin)) << ")"
                << " Y-max-ghost(flag=" << int(flagField.get(id_ymax)) << ")"
                << " BulkFlag=" << int(BulkFlag)
                << " PeriodicFlag=" << int(PeriodicFlag)
                << " BouncebackFlag=" << int(BouncebackFlag)
                << std::endl;
    }
  }

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

  // Fill X/Y periodic ghost cells with initial POP from opposite side.
  // Must be called before applyWallPop (which only resets Z ghost cells).
  PF_PeriodicX.Apply();
  PF_PeriodicY.Apply();
  NS_PeriodicX.Apply();
  NS_PeriodicY.Apply();

  // Reset Z-direction domain-boundary ghost-layer POP to equilibrium values.
  // NormalCommunicate only fills inter-block ghost cells (that have a
  // neighbour); Z domain-boundary ghost cells have no neighbour and retain
  // whatever stream left there.  Without this reset the corrupted POP
  // propagates inward each step, producing the boundary shadow.
  applyWallPop<T, LatSet, decltype(PFLattice), decltype(NSLattice)>(Geo, PFLattice, NSLattice, Cell_Len);

  // Communicate PHI so inter-block ghost cells are filled from neighbours
  // before the first gradient/normal/laplacian/chempot computation.
  // NormalCommunicate only exchanges POP, not PHI; without this call the
  // block-boundary ghost cells would retain phi=0 and corrupt the initial
  // gradients.
  PFLattice.getField<PHI<T>>().Communicate();

  // Set phi=1 on wall cells before the initial gradient computation.
  // T0 is written before the first main-loop E2, so without this call the
  // domain-boundary ghost cells would remain phi=0 and corrupt the initial
  // NORMAL field.  That error then feeds into the PF collision and produces
  // the boundary shadow that propagates inward.
  applyWallPhi<T, LatSet>(Geo, PFLattice.getField<PHI<T>>(), Cell_Len);

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
    NSLattice.getField<FORCE<T, LatSet::d>>().InitValue(Vector<T, LatSet::d>{T{0}});

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
    // Copy post-collision POP+PHI from opposite side to X/Y periodic ghost cells
    PF_PeriodicX.Apply();
    PF_PeriodicY.Apply();

    // ---- Phase C: NS collision (Fortran collisionDenDF2D) ----
    // C1: Viscous force F_v from non-equilibrium moments
    ViscoForceCoupling.ApplyInnerCellDynamics<ViscoForceTaskSelector>(MainLoopTimer(), FlagFM);

    // C2: MRTForce with total accumulated FORCE (F_s+F_b+F_p+F_v)
    NSLattice.template ApplyInnerCellDynamics<NSTaskSelector>(FlagFM);
    // Copy post-collision POP+DENSITY from opposite side to X/Y periodic ghost cells
    NS_PeriodicX.Apply();
    NS_PeriodicY.Apply();

    // ---- Phase D: Streaming (Fortran streamOrderDF + streamDenDF) ----
    // D1: PF bounceback + stream
    PF_BB.Apply(MainLoopTimer());
    PFLattice.Stream();
    PFLattice.NormalCommunicate();
    // D2: NS bounceback + stream
    NS_BB.Apply(MainLoopTimer());
    NSLattice.Stream();
    NSLattice.NormalCommunicate();

    // Reset domain-boundary ghost-layer POP after stream.
    // Stream (PULL) leaves domain-boundary ghost cells reading out-of-domain
    // memory; NormalCommunicate does not fill them (no neighbour).  Without
    // this reset the corrupted POP propagates inward each step, producing
    // the boundary shadow that diffuses toward the centre.
    applyWallPop<T, LatSet, decltype(PFLattice), decltype(NSLattice)>(Geo, PFLattice, NSLattice, Cell_Len);

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
              // NaN/Inf guard: IEEE 754 makes NaN<x and NaN>x both false,
              // so the [0,1] clamps below cannot catch NaN. isfinite first.
              if (!std::isfinite(phi_new)) phi_new = T{0};
              if (phi_new < T{0}) phi_new = T{0};
              if (phi_new > T{1}) phi_new = T{1};
              blockPhi.get(id) = phi_new;
            }
          }
        }
      }
    }
    PFLattice.getField<PHI<T>>().Communicate();


    // E2: phi=1 at all no-slip walls (including ghost layers, edges and corners)
    applyWallPhi<T, LatSet>(Geo, PFLattice.getField<PHI<T>>(), Cell_Len);
    // E2a: Re-communicate PHI after wall modification (ghost cells need updated phi)
    PFLattice.getField<PHI<T>>().Communicate();

    // E3: Gradients, normal, laplacian, chempot (Fortran computeOrderMacro)
    PFLattice.template ApplyInnerCellDynamics<FFNormalSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFLaplacianSelector>(FlagFM);
    PFLattice.template ApplyInnerCellDynamics<FFChemPotSelector>(FlagFM);
    PFLattice.getField<NORMAL<T, LatSet::d>>().Communicate();
    PFLattice.getField<GRAD<T, LatSet::d>>().Communicate();
    ff::CommunicateAllSelfFields<T>(PFLattice);

    // E3a: Wall grad_phi special handling: copy normal component from first interior cell
    // Only Z-direction walls (top/bottom).  X/Y are periodic, so grad at X/Y
    // boundaries is computed normally from periodic ghost-cell phi.
    {
      auto& gradField = PFLattice.getField<GRAD<T, LatSet::d>>();
      T Lz_global = T(Nz) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockGrad = gradField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int nz = block.getNz();
        int overlap = block.getOverlap();
        T minZ = block.getMin()[2];
        T maxZ = block.getMax()[2];
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

    // E4: Chempot extrapolation at Z walls (Fortran setMacroOrderBC chpoten)
    // Only Z-direction walls.  X/Y are periodic.
    {
      auto& chpotenField = PFLattice.getField<ff::CHEMICALPOTENTIAL<T>>();
      T Lz_global = T(Nz) * Cell_Len;
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockChpoten = chpotenField.getBlockField(blockid);
        int nx = block.getNx();
        int ny = block.getNy();
        int nz = block.getNz();
        int overlap = block.getOverlap();
        T minZ = block.getMin()[2];
        T maxZ = block.getMax()[2];
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

    // E5a: Zero out velocity at BouncebackFlag cells (wall cells should have u=0)
    // E5 computes u = u_raw + 0.5*F/rho for ALL internal cells including
    // BouncebackFlag, but BouncebackFlag cells are not collided (no MRTForce),
    // so their POP is just reflected values. The computed u = 0.5*F/rho = -0.5*g
    // is wrong (should be 0 for no-slip wall). This is a visualization artifact
    // that may appear as a boundary shadow in the Velocity field.
    {
      auto& velField = NSLattice.getField<VELOCITY<T, LatSet::d>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockVel = velField.getBlockField(blockid);
        auto& blockFlag = FlagFM.getBlockField(blockid);
        int overlap = block.getOverlap();
        for (int k = overlap; k < block.getNz() - overlap; ++k) {
          for (int j = overlap; j < block.getNy() - overlap; ++j) {
            for (int i = overlap; i < block.getNx() - overlap; ++i) {
              std::size_t id = k * proj[2] + j * proj[1] + i;
              if (util::isFlag(blockFlag.get(id), BouncebackFlag)) {
                blockVel.get(id) = Vector<T, LatSet::d>{T{0}};
              }
            }
          }
        }
      }
    }

    // ===== [DIAG] pressure profile along z-axis (x=48.5, y=48.5) =====
    // Checks pressure at bottom (z=0.5), bubble (z=48.5), mid (z=95.5), top (z=191.5)
    // Pressure wave travels at cs≈0.577 grid/step, so from bottom (z=0) to
    // bubble (z=48) takes ~83 steps.  Run for 200 steps to see gradient form
    // and verify p clamping stabilizes the feedback loop.
    if (MainLoopTimer() < 200 && MainLoopTimer() % 20 == 0) {
      T p_bottom = T{0}, p_bubble = T{0}, p_mid = T{0}, p_top = T{0};
      T pmin_l = T(1e30), pmax_l = T(-1e30);
      // Track pmax cell location for step-0 diagnosis
      T pmax_x = T{0}, pmax_y = T{0}, pmax_z = T{0};
      int pmax_block = -1, pmax_li = -1, pmax_lj = -1, pmax_lk = -1;
      auto& presField = NSLattice.getField<PRESSURE<T>>();
      for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
        const auto& block = Geo.getBlock(blockid);
        const auto& proj = block.getProjection();
        auto& blockPres = presField.getBlockField(blockid);
        int overlap = block.getOverlap();
        T gx = block.getMinCenter()[0];
        T gy = block.getMinCenter()[1];
        T gz = block.getMinCenter()[2];
        int nx = block.getNx(), ny = block.getNy(), nz = block.getNz();
        for (int k = overlap; k < nz - overlap; ++k) {
          for (int j = overlap; j < ny - overlap; ++j) {
            for (int i = overlap; i < nx - overlap; ++i) {
              T x = gx + T(i) * Cell_Len;
              T y = gy + T(j) * Cell_Len;
              T z = gz + T(k) * Cell_Len;
              std::size_t id = k * proj[2] + j * proj[1] + i;
              T p = blockPres.get(id);
              if (p < pmin_l) pmin_l = p;
              if (p > pmax_l) {
                pmax_l = p;
                pmax_x = x; pmax_y = y; pmax_z = z;
                pmax_block = blockid;
                pmax_li = i; pmax_lj = j; pmax_lk = k;
              }
              if (std::abs(x - T{48.5}) < T{0.5} && std::abs(y - T{48.5}) < T{0.5}) {
                if (std::abs(z - T{0.5}) < T{0.5}) p_bottom = p;
                else if (std::abs(z - T{48.5}) < T{0.5}) p_bubble = p;
                else if (std::abs(z - T{95.5}) < T{0.5}) p_mid = p;
                else if (std::abs(z - T{191.5}) < T{0.5}) p_top = p;
              }
            }
          }
        }
      }
#ifdef MPI_ENABLED
      mpi().reduceAndBcast<T>(pmin_l, MPI_MIN);
      mpi().reduceAndBcast<T>(pmax_l, MPI_MAX);
      mpi().reduceAndBcast<T>(p_bottom, MPI_SUM);
      mpi().reduceAndBcast<T>(p_bubble, MPI_SUM);
      mpi().reduceAndBcast<T>(p_mid, MPI_SUM);
      mpi().reduceAndBcast<T>(p_top, MPI_SUM);
#endif
      if (mpi().isMainProcessor()) {
        T dpdz_bubble = (p_top - p_bottom) / T{191};
        std::cout << "[DIAG] Step=" << MainLoopTimer()
                  << " p_range=[" << pmin_l << "," << pmax_l << "]"
                  << " p_bottom=" << p_bottom
                  << " p_bubble=" << p_bubble
                  << " p_mid=" << p_mid
                  << " p_top=" << p_top
                  << " dpdz_full=" << dpdz_bubble
                  << " (expected ~" << -rho_h * gravity << ")" << std::endl;
      }
      // ===== Step-0 pmax detailed diagnosis =====
      // Print pmax cell's POP values to identify the source of p=1/6 anomaly
      if (MainLoopTimer() == 0) {
        // Find pmax cell again on the rank that owns it
        T local_pmax = T(-1e30);
        int local_pmax_block = -1;
        int local_pmax_li = -1, local_pmax_lj = -1, local_pmax_lk = -1;
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& blockPres = presField.getBlockField(blockid);
          int overlap = block.getOverlap();
          int nx = block.getNx(), ny = block.getNy(), nz = block.getNz();
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              for (int i = overlap; i < nx - overlap; ++i) {
                std::size_t id = k * proj[2] + j * proj[1] + i;
                T p = blockPres.get(id);
                if (p > local_pmax) {
                  local_pmax = p;
                  local_pmax_block = blockid;
                  local_pmax_li = i; local_pmax_lj = j; local_pmax_lk = k;
                }
              }
            }
          }
        }
        // Gather all local pmax values to find the global winner
        struct PmaxInfo { T val; int rank; int block; int li; int lj; int lk; };
        PmaxInfo my_info{local_pmax, mpi().getRank(), local_pmax_block,
                         local_pmax_li, local_pmax_lj, local_pmax_lk};
#ifdef MPI_ENABLED
        PmaxInfo all_info[8];
        MPI_Allgather(&my_info, sizeof(PmaxInfo), MPI_BYTE,
                      all_info, sizeof(PmaxInfo), MPI_BYTE, MPI_COMM_WORLD);
        int winner = 0;
        for (int r = 1; r < mpi().getSize(); ++r) {
          if (all_info[r].val > all_info[winner].val) winner = r;
        }
#else
        PmaxInfo all_info[1] = {my_info};
        int winner = 0;
#endif
        // The winning rank prints its pmax cell details
        if (mpi().getRank() == all_info[winner].rank) {
          const auto& block = Geo.getBlock(all_info[winner].block);
          const auto& proj = block.getProjection();
          auto& blockLat = NSLattice.getBlockLat(all_info[winner].block);
          auto& blockPres2 = presField.getBlockField(all_info[winner].block);
          auto& blockFlag = FlagFM.getBlockField(all_info[winner].block);
          int overlap = block.getOverlap();
          int li = all_info[winner].li, lj = all_info[winner].lj, lk = all_info[winner].lk;
          std::size_t id = lk * proj[2] + lj * proj[1] + li;
          T gx = block.getMinCenter()[0], gy = block.getMinCenter()[1], gz = block.getMinCenter()[2];
          T x = gx + T(li) * Cell_Len, y = gy + T(lj) * Cell_Len, z = gz + T(lk) * Cell_Len;

          std::cout << "[PMAX_DIAG] Step=0 pmax=" << blockPres2.get(id)
                    << " at global=(" << x << "," << y << "," << z << ")"
                    << " block=" << all_info[winner].block
                    << " local=(" << li << "," << lj << "," << lk << ")"
                    << " overlap=" << overlap
                    << " nx=" << block.getNx() << " ny=" << block.getNy() << " nz=" << block.getNz()
                    << " flag=" << int(blockFlag.get(id))
                    << std::endl;

          // Print all 19 POP values
          NSCELL cell(id, blockLat);
          T p_sum = T{0};
          std::cout << "[PMAX_DIAG] POP values:" << std::endl;
          for (unsigned int qi = 0; qi < LatSet::q; ++qi) {
            const auto& ci = latset::c<LatSet>(qi);
            T wi = latset::w<LatSet>(qi);
            T fi = cell[qi];
            p_sum += fi;
            std::cout << "  i=" << qi
                      << " c=(" << int(ci[0]) << "," << int(ci[1]) << "," << int(ci[2]) << ")"
                      << " w=" << wi
                      << " f=" << fi
                      << " (f-w=" << (fi - wi) << ")"
                      << std::endl;
          }
          std::cout << "[PMAX_DIAG] sum_f=" << p_sum << " pressure_field=" << blockPres2.get(id) << std::endl;

          // Also check if this cell is near a block boundary
          int dist_to_min_x = li - overlap;
          int dist_to_min_y = lj - overlap;
          int dist_to_min_z = lk - overlap;
          int dist_to_max_x = block.getNx() - 1 - overlap - li;
          int dist_to_max_y = block.getNy() - 1 - overlap - lj;
          int dist_to_max_z = block.getNz() - 1 - overlap - lk;
          std::cout << "[PMAX_DIAG] dist_to_boundary: x_min=" << dist_to_min_x
                    << " x_max=" << dist_to_max_x
                    << " y_min=" << dist_to_min_y
                    << " y_max=" << dist_to_max_y
                    << " z_min=" << dist_to_min_z
                    << " z_max=" << dist_to_max_z
                    << " (0=first internal, 1=second internal)" << std::endl;

          // Check block's global position
          std::cout << "[PMAX_DIAG] block min=(" << block.getMin()[0] << "," << block.getMin()[1] << "," << block.getMin()[2] << ")"
                    << " max=(" << block.getMax()[0] << "," << block.getMax()[1] << "," << block.getMax()[2] << ")"
                    << " minCenter=(" << block.getMinCenter()[0] << "," << block.getMinCenter()[1] << "," << block.getMinCenter()[2] << ")"
                    << std::endl;
        }
      }
    }

    // ===== phi 切片诊断（每 20 步打印界面锐度统计）=====
    if (MainLoopTimer() % 20 == 0 && MainLoopTimer() > 0) {
        // 打印 y=48 切片的 phi 值（x-z 平面）
        auto& phiField = PFLattice.getField<PHI<T>>();
        T phi_min_l = T(2), phi_max_l = T(-1);
        T phi_sum_l = T{0};
        int phi_count_l = 0;
        // 统计 phi < 0.1 和 phi > 0.9 的格点数（检查界面锐度）
        int phi_low_count_l = 0;  // phi < 0.1（气泡内部）
        int phi_high_count_l = 0;  // phi > 0.9（重流体）
        int phi_mid_count_l = 0;   // 0.1 <= phi <= 0.9（界面）

        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
            const auto& block = Geo.getBlock(blockid);
            const auto& proj = block.getProjection();
            auto& blockPhi = phiField.getBlockField(blockid);
            int overlap = block.getOverlap();
            int nx = block.getNx(), ny = block.getNy(), nz = block.getNz();
            for (int k = overlap; k < nz - overlap; ++k) {
                for (int j = overlap; j < ny - overlap; ++j) {
                    for (int i = overlap; i < nx - overlap; ++i) {
                        std::size_t id = k * proj[2] + j * proj[1] + i;
                        T phi = blockPhi.get(id);
                        if (phi < phi_min_l) phi_min_l = phi;
                        if (phi > phi_max_l) phi_max_l = phi;
                        phi_sum_l += phi;
                        phi_count_l++;
                        if (phi < T{0.1}) phi_low_count_l++;
                        else if (phi > T{0.9}) phi_high_count_l++;
                        else phi_mid_count_l++;
                    }
                }
            }
        }
#ifdef MPI_ENABLED
        mpi().reduceAndBcast<T>(phi_min_l, MPI_MIN);
        mpi().reduceAndBcast<T>(phi_max_l, MPI_MAX);
        mpi().reduceAndBcast<T>(phi_sum_l, MPI_SUM);
        mpi().reduceAndBcast<int>(phi_count_l, MPI_SUM);
        mpi().reduceAndBcast<int>(phi_low_count_l, MPI_SUM);
        mpi().reduceAndBcast<int>(phi_high_count_l, MPI_SUM);
        mpi().reduceAndBcast<int>(phi_mid_count_l, MPI_SUM);
#endif
        if (mpi().isMainProcessor()) {
            T phi_avg = phi_sum_l / T(phi_count_l);
            std::cout << "[SHAPE] Step=" << MainLoopTimer()
                      << " phi_avg=" << phi_avg
                      << " low(<0.1)=" << phi_low_count_l
                      << " mid(0.1-0.9)=" << phi_mid_count_l
                      << " high(>0.9)=" << phi_high_count_l
                      << " ratio=" << T(phi_mid_count_l) / T(phi_count_l)
                      << std::endl;
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
      // === TEMP VERIFICATION: phi & density range check ===
      {
        auto& phiField = PFLattice.getField<PHI<T>>();
        auto& densField = NSLattice.getField<DENSITY<T>>();
        T phi_min_l = T(2), phi_max_l = T(-1);
        T dens_min_l = T(1e30), dens_max_l = T(-1e30);
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& blockPhi = phiField.getBlockField(blockid);
          auto& blockDens = densField.getBlockField(blockid);
          int overlap = block.getOverlap();
          for (int k = overlap; k < block.getNz() - overlap; ++k) {
            for (int j = overlap; j < block.getNy() - overlap; ++j) {
              for (int i = overlap; i < block.getNx() - overlap; ++i) {
                std::size_t id = k * proj[2] + j * proj[1] + i;
                T phi = blockPhi.get(id);
                T dens = blockDens.get(id);
                if (phi < phi_min_l) phi_min_l = phi;
                if (phi > phi_max_l) phi_max_l = phi;
                if (dens < dens_min_l) dens_min_l = dens;
                if (dens > dens_max_l) dens_max_l = dens;
              }
            }
          }
        }
#ifdef MPI_ENABLED
        mpi().reduceAndBcast<T>(phi_min_l, MPI_MIN);
        mpi().reduceAndBcast<T>(phi_max_l, MPI_MAX);
        mpi().reduceAndBcast<T>(dens_min_l, MPI_MIN);
        mpi().reduceAndBcast<T>(dens_max_l, MPI_MAX);
#endif
        if (mpi().isMainProcessor()) {
          std::cout << "[Step " << MainLoopTimer() << "] phi range: [" << phi_min_l << ", " << phi_max_l
                    << "]  dens range: [" << dens_min_l << ", " << dens_max_l << "]" << std::endl;
        }
      }
      // === END TEMP VERIFICATION ===

      // 气泡形状监控：计算质心、等效半径、体积（使用全局/物理坐标）
      // 与气泡初始化一致：x = block.getMinCenter()[0] + i*voxelSize
      // 这样跨块 MPI_SUM 才有意义，且中心可直接与 Bubble_Center=(48,48,48) 比较
      {
        T sum_phi = T{0};
        T sum_phi_x = T{0}, sum_phi_y = T{0}, sum_phi_z = T{0};
        int bubble_count = 0;

        auto& phiField = PFLattice.getField<PHI<T>>();
        for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
          const auto& block = Geo.getBlock(blockid);
          const auto& proj = block.getProjection();
          auto& blockPhi = phiField.getBlockField(blockid);
          int overlap = block.getOverlap();
          int nx = block.getNx(), ny = block.getNy(), nz = block.getNz();
          T voxelSize = block.getVoxelSize();
          // 块的全局起始物理坐标（含 overlap 偏移），与初始化公式完全一致
          T gx0 = block.getMinCenter()[0];
          T gy0 = block.getMinCenter()[1];
          T gz0 = block.getMinCenter()[2];
          for (int k = overlap; k < nz - overlap; ++k) {
            for (int j = overlap; j < ny - overlap; ++j) {
              for (int i = overlap; i < nx - overlap; ++i) {
                std::size_t id = k * proj[2] + j * proj[1] + i;
                T phi = blockPhi.get(id);
                // phi < 0.5 是气泡内部（轻流体）
                T bubble_indicator = T{1} - phi;  // 气泡指示器：1=气泡内，0=气泡外
                if (bubble_indicator > T{0.5}) {
                  bubble_count++;
                  sum_phi += bubble_indicator;
                  // 全局物理坐标 = MinCenter + 局部索引 * voxelSize
                  sum_phi_x += bubble_indicator * (gx0 + T(i) * voxelSize);
                  sum_phi_y += bubble_indicator * (gy0 + T(j) * voxelSize);
                  sum_phi_z += bubble_indicator * (gz0 + T(k) * voxelSize);
                }
              }
            }
          }
        }

        // MPI reduce
#ifdef MPI_ENABLED
        mpi().reduceAndBcast<T>(sum_phi, MPI_SUM);
        mpi().reduceAndBcast<T>(sum_phi_x, MPI_SUM);
        mpi().reduceAndBcast<T>(sum_phi_y, MPI_SUM);
        mpi().reduceAndBcast<T>(sum_phi_z, MPI_SUM);
        mpi().reduceAndBcast<int>(bubble_count, MPI_SUM);
#endif

        if (mpi().isMainProcessor() && bubble_count > 0) {
          T cx = sum_phi_x / sum_phi;
          T cy = sum_phi_y / sum_phi;
          T cz = sum_phi_z / sum_phi;
          T volume = sum_phi;  // 近似体积
          T r_eff = std::pow(T{3} * volume / (T{4} * M_PI), T{1}/T{3});  // 等效半径
          std::cout << "[BUBBLE] Step=" << MainLoopTimer()
                    << " center=(" << cx << "," << cy << "," << cz << ")"
                    << " volume=" << volume
                    << " r_eff=" << r_eff
                    << " count=" << bubble_count
                    << std::endl;
        }
      }

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
