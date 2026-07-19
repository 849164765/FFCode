// test_mrt_diffusion.cpp — 1000-step test for MRTDiffusion collision (D2Q5)
// Tests collision in isolation (no streaming): psi conservation + relaxation
#include "freelb.h"
#include "freelb.hh"

using T = double;
using MFLatSet = D2Q5<T>;
using namespace mfield;

constexpr int N = 8;
constexpr T Cell_Len = 1.0;
constexpr int BlockCellLen = 8;
constexpr T omega_test = T{1.0};
constexpr T psi_init = T{1.0};
constexpr T tol = 1e-12;

int main() {
  constexpr std::uint8_t VoidFlag = std::uint8_t(1);
  constexpr std::uint8_t BulkFlag = std::uint8_t(2);

  // ------------------ geometry ------------------
  AABB<T, 2> domain(Vector<T, 2>{T(0), T(0)},
                    Vector<T, 2>{T(N * Cell_Len), T(N * Cell_Len)});
  BlockGeometryHelper2D<T> GeoHelper(N, N, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, 1);
  GeoHelper.AdaptiveOptimization(1);
  GeoHelper.LoadBalancing(1);
  BlockGeometry2D<T> Geo(GeoHelper);

  // ------------------ flag ------------------
  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(domain, [&](FLAG& f, std::size_t id) { f.SetField(id, BulkFlag); });

  // ------------------ MF lattice ------------------
  BaseConverter<T> MFBaseConv(MFLatSet::cs2);
  MFBaseConv.SimplifiedConverterFromRT(N, T(0.01), T{1.0} / omega_test);

  using MFF = TypePack<PSI<T>, OMEGA_PSI<T>, POP<T, MFLatSet::q>>;
  ValuePack MFInit(T{psi_init}, omega_test, T{});
  using MFCELL = Cell<T, MFLatSet, MFF>;
  BlockLatticeManager<T, MFLatSet, MFF> MFLattice(Geo, MFInit, MFBaseConv);

  // ------------------ init psi=1.0, g_k = random perturbation from eq ------------------
  MFLattice.getField<PSI<T>>().InitValue(psi_init);
  MFLattice.getField<OMEGA_PSI<T>>().InitValue(omega_test);

  // Store initial non-eq measure: sum(|g_k - w_k*psi|)
  {
    auto& psiField = MFLattice.getField<PSI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      auto& blockLat = MFLattice.getBlockLat(blockid);
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& blockPsi = psiField.getBlockField(blockid);
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          std::size_t id = j * proj[1] + i;
          MFCELL cell(id, blockLat);
          T psi_cur = psi_init;
          // Perturb: g_k = w_k*psi + random deviation (mass-conserving)
          T accum = T{0};
          for (unsigned int k = 0; k < MFLatSet::q; ++k) {
            // 10% perturbation on each direction except k=0 (which absorbs mass)
            T w = latset::w<MFLatSet>(k);
            T pert = (k > 0) ? T{0.1} * w * psi_cur * (k%2 ? T{1} : T{-1}) : T{0};
            cell[k] = w * psi_cur + pert;
            accum += pert;
          }
          // Correct k=0 to conserve mass
          cell[0] -= accum;
          blockPsi.get(id) = psi_cur;
        }
      }
    }
  }

  // Run 1000 collision steps (no streaming — collision-only test)
  Timer timer;
  for (int step = 0; step < 1000; ++step) {
    // Direct cell iteration — bypass TaskSelector for clarity
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      auto& blockLat = MFLattice.getBlockLat(blockid);
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          std::size_t id = j * proj[1] + i;
          MFCELL cell(id, blockLat);
          collision::MRTDiffusion<MFCELL, OMEGA_PSI<T>>::apply(cell);
        }
      }
    }
    // After collision, update PSI = Σg_i
    auto& psiField = MFLattice.getField<PSI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      auto& blockLat = MFLattice.getBlockLat(blockid);
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& blockPsi = psiField.getBlockField(blockid);
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          std::size_t id = j * proj[1] + i;
          MFCELL cell(id, blockLat);
          T psi = T{0};
          for (unsigned int k = 0; k < MFLatSet::q; ++k) psi += cell[k];
          blockPsi.get(id) = psi;
        }
      }
    }
    ++timer;
  }

  // ------------------ verify ------------------
  T max_psi_err = T{0};
  T sum_psi = T{0};
  int count = 0;

  // Also check that populations have relaxed to equilibrium
  T max_g_err = T{0};

  {
    auto& psiField = MFLattice.getField<PSI<T>>();
    for (int blockid = 0; blockid < Geo.getBlockNum(); ++blockid) {
      auto& blockLat = MFLattice.getBlockLat(blockid);
      const auto& block = Geo.getBlock(blockid);
      const auto& proj = block.getProjection();
      auto& blockPsi = psiField.getBlockField(blockid);
      for (int j = 0; j < block.getNy(); ++j) {
        for (int i = 0; i < block.getNx(); ++i) {
          std::size_t id = j * proj[1] + i;
          MFCELL cell(id, blockLat);
          T psi_val = blockPsi.get(id);
          T err = std::abs(psi_val - psi_init);
          if (err > max_psi_err) max_psi_err = err;
          sum_psi += psi_val;
          ++count;

          // g_k should be close to w_k * psi
          for (unsigned int k = 0; k < MFLatSet::q; ++k) {
            T ge = std::abs(cell[k] - latset::w<MFLatSet>(k) * psi_val);
            if (ge > max_g_err) max_g_err = ge;
          }
        }
      }
    }
  }

  T avg_psi = sum_psi / T(count);
  bool psi_ok = (max_psi_err < tol);
  bool eq_ok  = (max_g_err < T{1e-6});  // after 1000 collisions, should be near equilibrium

  std::cout << "Step 1000 (collision only):" << std::endl;
  std::cout << "  max|psi-1.0| = " << max_psi_err << (psi_ok ? " OK" : " FAIL") << std::endl;
  std::cout << "  avg_psi      = " << avg_psi << std::endl;
  std::cout << "  max|g-w*psi| = " << max_g_err << (eq_ok ? " OK" : " FAIL") << std::endl;
  std::cout << "  cells        = " << count << std::endl;
  std::cout << (psi_ok && eq_ok ? "PASS" : "FAIL") << std::endl;

  return (psi_ok && eq_ok) ? 0 : 1;
}
