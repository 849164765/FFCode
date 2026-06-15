// tests/test_pf.cpp
// Unit tests for pf/pf2d.h — Phase field physics functors (D2Q9)

#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>
#include <vector>
#include <cstdlib>

#include "utils/alias.h"
#include "lbm/lattice_set.h"
#include "pf/pf2d.h"

using T = double;
using LatType = D2Q9<T>;
static constexpr T TOL = 1e-14;
static int failures = 0;

#define ASSERT_NEAR(val, expected, tol)                                      \
  do {                                                                       \
    auto _v = (val); auto _e = (expected);                                   \
    auto _d = (_v > _e) ? (_v - _e) : (_e - _v);                            \
    if (_d > (tol)) {                                                        \
      std::cerr << "FAIL: " << #val << "=" << _v << " exp=" << _e            \
                << " diff=" << _d << " @" << __FILE__ << ":" << __LINE__ << "\n"; \
      ++failures;                                                            \
    }                                                                        \
  } while (0)

// ===========================================================================
// MockCell with explicit per-direction neighbor phi values
//
// Rather than building a full 2D grid, the caller populates phiNbr[k]
// with the phi value the gradient/laplacian should see in direction k.
// getNeighbor(k) returns a thin proxy whose get<PHI<T>>() yields phiNbr[k].
// ===========================================================================
struct PfMockCell {
  using FloatType  = T;
  using LatticeSet = LatType;
  static constexpr unsigned int D = 2;
  static constexpr unsigned int Q = 9;

  T  phi_{};
  T  ifaceWidth_{};
  Vector<T, D> u_{};
  Vector<T, D> grad_{};
  Vector<T, D> normal_{};
  T  phiNbr_[9]{};

  T getOmega() const { return 1.5; }

  template <typename FieldType, unsigned int = 0>
  auto& get() {
    if constexpr (std::is_same_v<FieldType, PHI<T>>)          return phi_;
    else if constexpr (std::is_same_v<FieldType, INTERFACEWIDTH<T>>) return ifaceWidth_;
    else if constexpr (std::is_same_v<FieldType, VELOCITY<T, D>>)    return u_;
    else if constexpr (std::is_same_v<FieldType, GRAD<T, D>>)        return grad_;
    else if constexpr (std::is_same_v<FieldType, NORMAL<T, D>>)      return normal_;
    else { static_assert(sizeof(FieldType) == 0, "unsupported"); }
  }

  // Neighbor proxy
  struct NbrProxy {
    T phi;
    template <typename FT, unsigned int = 0> T& get() { return phi; }
    template <typename FT, unsigned int = 0> const T& get() const { return phi; }
  };
  NbrProxy getNeighbor(unsigned int k) const { return NbrProxy{phiNbr_[k]}; }
};

static void printHeader(const std::string& s) {
  std::cout << "\n=== " << s << " ===" << std::endl;
}

// ===========================================================================
// 1. PFGradient2D — exact linear profile φ = a·x + b·y + c
//    grad = (a, b)
// ===========================================================================
static int testGradientLinear() {
  printHeader("PFGradient2D linear phi=2x+3y+1");
  PfMockCell cell;
  cell.phi_ = T(1);  // phi at origin
  // Set neighbor values: φ(x+e_i) = 2*c_i_x + 3*c_i_y + 1
  // c from codebase: (0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)
  for (unsigned int k = 0; k < 9; ++k) {
    const auto& ck = latset::c<LatType>(k);
    cell.phiNbr_[k] = T(2)*T(ck[0]) + T(3)*T(ck[1]) + T(1);
  }
  pf::PFGradient2D<PfMockCell>::apply(cell);
  // Expected: grad = (2, 3)
  ASSERT_NEAR(cell.grad_[0], T(2), TOL);
  ASSERT_NEAR(cell.grad_[1], T(3), TOL);
  std::cout << "  grad=(" << cell.grad_[0] << "," << cell.grad_[1] << ")" << std::endl;
  return 0;
}

// ===========================================================================
// 2. PFGradient2D uniform field → zero gradient
// ===========================================================================
static int testGradientUniform() {
  printHeader("PFGradient2D uniform");
  PfMockCell cell;
  cell.phi_ = T(0.8);
  for (unsigned int k = 0; k < 9; ++k) cell.phiNbr_[k] = T(0.8);
  pf::PFGradient2D<PfMockCell>::apply(cell);
  ASSERT_NEAR(cell.grad_[0], T(0), TOL);
  ASSERT_NEAR(cell.grad_[1], T(0), TOL);
  std::cout << "  passed." << std::endl;
  return 0;
}

// ===========================================================================
// 3. PFLaplacian2D φ = x² + y² → ∇² = 4
//    φ(x+c) = (c_x)² + (c_y)²,  φ(0,0) = 0
// ===========================================================================
static int testLaplacianQuadratic() {
  printHeader("PFLaplacian2D phi=x^2+y^2");
  PfMockCell cell;
  cell.phi_ = T(0);
  for (unsigned int k = 0; k < 9; ++k) {
    const auto& ck = latset::c<LatType>(k);
    cell.phiNbr_[k] = T(ck[0]*ck[0] + ck[1]*ck[1]);
  }
  T lap = pf::PFLaplacian2D<PfMockCell>::compute(cell);
  // Analytical: ∇²(x²+y²) = 4
  ASSERT_NEAR(lap, T(4), T(0.01));
  std::cout << "  lap=" << lap << " (exp 4)" << std::endl;
  return 0;
}

// ===========================================================================
// 4. PFLaplacian2D φ = x² → ∇² = 2
// ===========================================================================
static int testLaplacianX2() {
  printHeader("PFLaplacian2D phi=x^2");
  PfMockCell cell;
  cell.phi_ = T(0);
  for (unsigned int k = 0; k < 9; ++k) {
    const auto& ck = latset::c<LatType>(k);
    cell.phiNbr_[k] = T(ck[0]*ck[0]);
  }
  T lap = pf::PFLaplacian2D<PfMockCell>::compute(cell);
  ASSERT_NEAR(lap, T(2), T(0.01));
  std::cout << "  lap=" << lap << " (exp 2)" << std::endl;
  return 0;
}

// ===========================================================================
// 5. PFChemPotential2D φ = 0.5 (constant, equilibrium)
//    μ = 4β*0.5*(0.5-1)*(0.5-0.5) - κ*0 = 0
// ===========================================================================
static int testChemPotentialEquilibrium() {
  printHeader("PFChemPotential phi=0.5 equilibrium");
  pf::Beta<T> = T(12) * T(0.072) / T(4);
  pf::Kappa<T> = T(3) * T(4) * T(0.072) * T(0.5);

  PfMockCell cell;
  cell.phi_ = T(0.5);
  for (unsigned int k = 0; k < 9; ++k) cell.phiNbr_[k] = T(0.5);
  pf::PFChemPotential2D<PfMockCell>::apply(cell);
  ASSERT_NEAR(cell.ifaceWidth_, T(0), TOL);
  std::cout << "  mu=" << cell.ifaceWidth_ << " (exp 0)" << std::endl;
  return 0;
}

// ===========================================================================
// 6. PFNormal2D grad=(3,4) → n=(0.6, 0.8)
// ===========================================================================
static int testNormal() {
  printHeader("PFNormal2D grad=(3,4)");
  PfMockCell cell;
  cell.grad_ = {T(3), T(4)};
  pf::PFNormal2D<PfMockCell>::apply(cell);
  ASSERT_NEAR(cell.normal_[0], T(0.6), TOL);
  ASSERT_NEAR(cell.normal_[1], T(0.8), TOL);
  std::cout << "  n=(" << cell.normal_[0] << "," << cell.normal_[1] << ")" << std::endl;
  return 0;
}

// ===========================================================================
// 7. PFNormal2D zero gradient → zero vector
// ===========================================================================
static int testNormalZero() {
  printHeader("PFNormal2D zero grad");
  PfMockCell cell;
  cell.grad_ = {T(0), T(0)};
  pf::PFNormal2D<PfMockCell>::apply(cell);
  ASSERT_NEAR(cell.normal_[0], T(0), TOL);
  ASSERT_NEAR(cell.normal_[1], T(0), TOL);
  std::cout << "  passed." << std::endl;
  return 0;
}

// ===========================================================================
int main() {
  std::cout << std::setprecision(16);
  std::cout << "========================================================\n";
  std::cout << "  pf/pf2d.h Unit Tests\n";
  std::cout << "========================================================" << std::endl;
  testGradientLinear();
  testGradientUniform();
  testLaplacianQuadratic();
  testLaplacianX2();
  testChemPotentialEquilibrium();
  testNormal();
  testNormalZero();
  std::cout << "\n========================================================" << std::endl;
  if (failures == 0)
    std::cout << "  ALL TESTS PASSED" << std::endl;
  else
    std::cout << "  " << failures << " TEST(S) FAILED" << std::endl;
  std::cout << "========================================================" << std::endl;
  return failures;
}
