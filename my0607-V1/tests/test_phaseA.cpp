// test_phaseA.cpp — Phase A gradient & Laplacian tests
// Validates Eqs. 34 & 35 against analytical solutions
//
// Test 1: φ(x,y) = x²  →  ∇φ = (2x, 0),  ∇²φ = 2
// Test 2: φ(x,y) = 3x - 2y  →  ∇φ = (3, -2),  ∇²φ = 0

#include "freelb.h"
#include "freelb.hh"
#include "lbm/phase_gradient.ur.h"

using T = double;
using LatSet = D2Q9<T>;

static int runTest(const char* name,
                   void (*initFunc)(BlockLatticeManager<T, LatSet, TypePack<PHI<T>, GRAD<T, 2>, ff::LAPLACIAN<T>>>&,
                                   const BlockGeometry2D<T>&),
                   T expected_lap,
                   void (*expectedGrad)(T x, T y, Vector<T, 2>& out)) {
  constexpr std::uint8_t VoidFlag = 1;
  constexpr std::uint8_t BulkFlag = 2;
  constexpr int MeshN = 8;
  const T Cell_Len = T{1};
  constexpr int BlockCellLen = 8;

  AABB<T, 2> domain(Vector<T, 2>(T{0}, T{0}),
                    Vector<T, 2>(T(MeshN * Cell_Len), T(MeshN * Cell_Len)));
  BlockGeometryHelper2D<T> GeoHelper(MeshN, MeshN, domain, Cell_Len, BlockCellLen);
  GeoHelper.CreateBlocks(1, mpi().getSize());
  GeoHelper.AdaptiveOptimization(mpi().getSize());
  GeoHelper.LoadBalancing(mpi().getSize());
  BlockGeometry2D<T> Geo(GeoHelper);

  if (Geo.getBlockNum() == 0) {
    std::cerr << "[" << name << "] FAIL: no blocks" << std::endl;
    return 1;
  }

  BlockFieldManager<FLAG, T, 2> FlagFM(Geo, VoidFlag);
  FlagFM.forEach(
    domain, [&](FLAG& field, std::size_t id) { field.SetField(id, BulkFlag); });
  FlagFM.template SetupBoundary<LatSet>(domain, BulkFlag);

  BaseConverter<T> dummyConv(LatSet::cs2);
  dummyConv.SimplifiedConverterFromRT(MeshN, T{1}, T{1.0});

  using PFIELDS = TypePack<PHI<T>, GRAD<T, 2>, ff::LAPLACIAN<T>>;
  ValuePack InitValues(T{}, Vector<T, 2>{T{0}, T{0}}, T{});
  BlockLatticeManager<T, LatSet, PFIELDS> Lat(Geo, InitValues, dummyConv);

  initFunc(Lat, Geo);

  using MyCell = Cell<T, LatSet, PFIELDS>;
  Lat.template ApplyInnerCellDynamics<phase_gradient::ComputeGradientPhi<MyCell>>();
  Lat.template ApplyInnerCellDynamics<phase_gradient::ComputeLaplacianPhi<MyCell>>();

  int pass = 0, fail = 0;
  T eps = T{1e-10};
  int ol = Geo.getBlock(0).getOverlap();
  const auto& blk = Geo.getBlock(0);
  const auto& proj = blk.getProjection();
  T vs = blk.getVoxelSize();
  T mx = blk.getMin()[0];
  T my = blk.getMin()[1];
  auto& blockLat = Lat.getBlockLat(0);

  for (int j = ol; j < blk.getNy() - ol; ++j) {
    for (int i = ol; i < blk.getNx() - ol; ++i) {
      std::size_t id = j * proj[1] + i;
      T x = mx + static_cast<T>(i) * vs;
      T y = my + static_cast<T>(j) * vs;

      MyCell cell(id, blockLat);
      const auto& cgrad = cell.template get<GRAD<T, 2>>();
      T clap = cell.template get<ff::LAPLACIAN<T>>();

      Vector<T, 2> egrad;
      expectedGrad(x, y, egrad);

      bool gok = true;
      for (unsigned int d = 0; d < 2; ++d) {
        if (std::abs(cgrad[d] - egrad[d]) > eps) gok = false;
      }
      bool lok = std::abs(clap - expected_lap) < eps;

      if (gok && lok) {
        pass++;
      } else {
        fail++;
        std::cerr << "[" << name << "] FAIL at i=" << i << " j=" << j
                  << " x=" << x << " y=" << y
                  << " grad=(" << cgrad[0] << "," << cgrad[1]
                  << ") exp=(" << egrad[0] << "," << egrad[1]
                  << ") lap=" << clap << " exp=" << expected_lap << std::endl;
      }
    }
  }

  if (mpi().isMainProcessor()) {
    std::cout << "[" << name << "] " << pass << " passed, " << fail << " failed" << std::endl;
  }
  return fail;
}


// ---- Test 1: φ(x,y) = x² ----
namespace test1 {
void init(BlockLatticeManager<T, LatSet, TypePack<PHI<T>, GRAD<T, 2>, ff::LAPLACIAN<T>>>& Lat,
          const BlockGeometry2D<T>& Geo) {
  auto& phiField = Lat.template getField<PHI<T>>();
  for (int bi = 0; bi < Geo.getBlockNum(); ++bi) {
    auto& bPhi = phiField.getBlockField(bi);
    const auto& b = Geo.getBlock(bi);
    const auto& p = b.getProjection();
    T vs = b.getVoxelSize();
    T mx = b.getMin()[0];
    for (int j = 0; j < b.getNy(); ++j)
      for (int i = 0; i < b.getNx(); ++i) {
        std::size_t id = j * p[1] + i;
        T x = mx + static_cast<T>(i) * vs;
        bPhi.get(id) = x * x;
      }
  }
}
void expectedGrad(T x, T /*y*/, Vector<T, 2>& out) {
  out[0] = T{2} * x;
  out[1] = T{0};
}
}  // namespace test1


// ---- Test 2: φ(x,y) = 3x - 2y ----
namespace test2 {
void init(BlockLatticeManager<T, LatSet, TypePack<PHI<T>, GRAD<T, 2>, ff::LAPLACIAN<T>>>& Lat,
          const BlockGeometry2D<T>& Geo) {
  auto& phiField = Lat.template getField<PHI<T>>();
  for (int bi = 0; bi < Geo.getBlockNum(); ++bi) {
    auto& bPhi = phiField.getBlockField(bi);
    const auto& b = Geo.getBlock(bi);
    const auto& p = b.getProjection();
    T vs = b.getVoxelSize();
    T mx = b.getMin()[0];
    T my = b.getMin()[1];
    for (int j = 0; j < b.getNy(); ++j)
      for (int i = 0; i < b.getNx(); ++i) {
        std::size_t id = j * p[1] + i;
        T x = mx + static_cast<T>(i) * vs;
        T y = my + static_cast<T>(j) * vs;
        bPhi.get(id) = T{3} * x - T{2} * y;
      }
  }
}
void expectedGrad(T /*x*/, T /*y*/, Vector<T, 2>& out) {
  out[0] = T{3};
  out[1] = T{-2};
}
}  // namespace test2


int main(int argc, char* argv[]) {
  mpi().init(&argc, &argv);

  int total_fail = 0;
  total_fail += runTest("test1:x^2",
                        test1::init, T{2}, test1::expectedGrad);
  total_fail += runTest("test2:3x-2y",
                        test2::init, T{0}, test2::expectedGrad);

  if (mpi().isMainProcessor()) {
    if (total_fail == 0)
      std::cout << "\n=== Phase A: ALL TESTS PASSED ===" << std::endl;
    else
      std::cout << "\n=== Phase A: " << total_fail << " TEST(S) FAILED ===" << std::endl;
  }
  return total_fail;
}
