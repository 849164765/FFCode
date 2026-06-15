// tests/test_collision.cpp
// Unit tests for lbm/collision.h (D2Q9 Lattice Boltzmann collision module)
//
// Test-Driven Development: RED phase
// These tests will fail to compile until collision.h is created.
//
// Compile:
//   g++ -std=c++17 -I$PROJECT_ROOT/src/ -DThread_Num=4 \
//       tests/test_collision.cpp -o tests/test_collision
//
// Run:
//   tests/test_collision

#include <cmath>
#include <iostream>
#include <iomanip>
#include <array>
#include <cstdlib>

// ---------------------------------------------------------------------------
// Project headers (lattice constants, field types, Vector)
// ---------------------------------------------------------------------------
// alias.h must precede lattice_set.h because tmp.h (included via util.h -> lattice_set.h)
// references RHO<T>, PHI<T>, etc. without forward-declaring them.
#include "utils/alias.h"
#include "lbm/lattice_set.h"

// ---------------------------------------------------------------------------
// Module under test (will NOT compile until collision.h exists -- RED phase)
// ---------------------------------------------------------------------------
#include "lbm/collision.h"
#include "lbm/force.h"

// ---------------------------------------------------------------------------
// Test tolerance
// ---------------------------------------------------------------------------
static constexpr double TOL = 1e-14;

// ---------------------------------------------------------------------------
// Assertion macro
// ---------------------------------------------------------------------------
#define ASSERT_NEAR(val, expected, tol)                                        \
  do {                                                                         \
    auto _v = (val);                                                           \
    auto _e = (expected);                                                      \
    auto _d = (_v > _e) ? (_v - _e) : (_e - _v);                               \
    if (_d > (tol)) {                                                          \
      std::cerr << "FAIL: " << #val << " = " << std::setprecision(16) << _v   \
                << " expected " << _e << " +/- " << (tol)                      \
                << " (diff=" << _d << ")"                                      \
                << " at " << __FILE__ << ":" << __LINE__ << std::endl;         \
      ++failures;                                                              \
    }                                                                          \
  } while (0)

#define ASSERT_TRUE(cond)                                                      \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL: " << #cond << " is false"                            \
                << " at " << __FILE__ << ":" << __LINE__ << std::endl;         \
      ++failures;                                                              \
    }                                                                          \
  } while (0)

// ---------------------------------------------------------------------------
// Helper for static_assert in discarded if constexpr branches
// ---------------------------------------------------------------------------
template <typename T>
struct dependent_false : std::false_type {};

// ===========================================================================
// Mock Cell
// ===========================================================================
// Implements the Cell interface required by collision.h without needing the
// real BlockLattice / GenericField infrastructure.
template <typename T, typename LatSet_>
struct MockCell {
  using FloatType   = T;
  using LatticeSet  = LatSet_;
  static constexpr unsigned int D = LatSet_::d;
  static constexpr unsigned int Q = LatSet_::q;

  // Stored fields
  T                         rho_{};
  T                         phi_{};
  Vector<T, D>              u_{};
  Vector<T, D>              force_{};       // for FORCE source field
  Vector<T, D>              constForce_{};  // for CONSTFORCE source field
  T                         omega_{};
  std::array<T, Q>          f_{};

  // Population access
  T& operator[](int k)             { return f_[k]; }
  const T& operator[](int k) const { return f_[k]; }

  // Relaxation frequency
  T getOmega() const { return omega_; }
  T get_Omega() const { return omega_; }
  T getfOmega() const { return omega_; }

  // -----------------------------------------------------------------------
  // Field access
  // Detects the template FieldType via if constexpr + std::is_same_v.
  // Supports: PHI<T>, RHO<T>, VELOCITY<T,D>, FORCE<T,D>, CONSTFORCE<T,D>
  // -----------------------------------------------------------------------
  template <typename FieldType, unsigned int /*i*/ = 0>
  auto& get() {
    if constexpr (std::is_same_v<FieldType, PHI<T>>) {
      return phi_;
    } else if constexpr (std::is_same_v<FieldType, RHO<T>>) {
      return rho_;
    } else if constexpr (std::is_same_v<FieldType, VELOCITY<T, D>>) {
      return u_;
    } else if constexpr (std::is_same_v<FieldType, FORCE<T, D>>) {
      return force_;
    } else if constexpr (std::is_same_v<FieldType, CONSTFORCE<T, D>>) {
      return constForce_;
    } else {
      static_assert(dependent_false<FieldType>::value,
                    "MockCell: unsupported FieldType in get()");
    }
  }

  // Const overload (not strictly needed by collision.h but good practice)
  template <typename FieldType, unsigned int /*i*/ = 0>
  const auto& get() const {
    if constexpr (std::is_same_v<FieldType, PHI<T>>) {
      return phi_;
    } else if constexpr (std::is_same_v<FieldType, RHO<T>>) {
      return rho_;
    } else if constexpr (std::is_same_v<FieldType, VELOCITY<T, D>>) {
      return u_;
    } else if constexpr (std::is_same_v<FieldType, FORCE<T, D>>) {
      return force_;
    } else if constexpr (std::is_same_v<FieldType, CONSTFORCE<T, D>>) {
      return constForce_;
    } else {
      static_assert(dependent_false<FieldType>::value,
                    "MockCell: unsupported FieldType in get() const");
    }
  }
};

// ===========================================================================
// Test helpers
// ===========================================================================

// Convenience alias for the D2Q9 double-precision lattice
using LatType = D2Q9<double>;
using Cell = MockCell<double, LatType>;
using TestD2Q9 = LatType;  // for use in testD2Q9Constants

// ---------------------------------------------------------------------------
// Print a test section header
// ---------------------------------------------------------------------------
static int failures = 0;

static void printHeader(const std::string& name) {
  std::cout << "\n=== " << name << " ===" << std::endl;
}

// ===========================================================================
// Tests
// ===========================================================================

// ---------------------------------------------------------------------------
// 1. D2Q9 lattice constant sanity checks
// ---------------------------------------------------------------------------
static int testD2Q9Constants() {
  printHeader("D2Q9 Lattice Constants");

  // Discrete velocities c[k] (codebase ordering)
  // k=0: (0,0)
  ASSERT_NEAR(latset::c<LatType>(0)[0], 0, TOL);
  ASSERT_NEAR(latset::c<LatType>(0)[1], 0, TOL);
  // k=1: (1,0)
  ASSERT_NEAR(latset::c<LatType>(1)[0], 1, TOL);
  ASSERT_NEAR(latset::c<LatType>(1)[1], 0, TOL);
  // k=2: (-1,0)
  ASSERT_NEAR(latset::c<LatType>(2)[0], -1, TOL);
  ASSERT_NEAR(latset::c<LatType>(2)[1], 0, TOL);
  // k=3: (0,1)
  ASSERT_NEAR(latset::c<LatType>(3)[0], 0, TOL);
  ASSERT_NEAR(latset::c<LatType>(3)[1], 1, TOL);
  // k=4: (0,-1)
  ASSERT_NEAR(latset::c<LatType>(4)[0], 0, TOL);
  ASSERT_NEAR(latset::c<LatType>(4)[1], -1, TOL);
  // k=5: (1,1)
  ASSERT_NEAR(latset::c<LatType>(5)[0], 1, TOL);
  ASSERT_NEAR(latset::c<LatType>(5)[1], 1, TOL);
  // k=6: (-1,-1)
  ASSERT_NEAR(latset::c<LatType>(6)[0], -1, TOL);
  ASSERT_NEAR(latset::c<LatType>(6)[1], -1, TOL);
  // k=7: (1,-1)
  ASSERT_NEAR(latset::c<LatType>(7)[0], 1, TOL);
  ASSERT_NEAR(latset::c<LatType>(7)[1], -1, TOL);
  // k=8: (-1,1)
  ASSERT_NEAR(latset::c<LatType>(8)[0], -1, TOL);
  ASSERT_NEAR(latset::c<LatType>(8)[1], 1, TOL);

  // Weights  w[k]
  ASSERT_NEAR(latset::w<LatType>(0), 4.0/9.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(1), 1.0/9.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(2), 1.0/9.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(3), 1.0/9.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(4), 1.0/9.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(5), 1.0/36.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(6), 1.0/36.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(7), 1.0/36.0, TOL);
  ASSERT_NEAR(latset::w<LatType>(8), 1.0/36.0, TOL);

  // Lattice constants
  ASSERT_NEAR(LatType::cs2, 1.0/3.0, TOL);
  ASSERT_NEAR(LatType::InvCs2, 3.0, TOL);
  ASSERT_NEAR(LatType::d, 2, TOL);
  ASSERT_NEAR(LatType::q, 9, TOL);

  // Sum of weights should be 1
  double wsum = 0.0;
  for (unsigned int k = 0; k < 9; ++k) wsum += latset::w<LatType>(k);
  ASSERT_NEAR(wsum, 1.0, TOL);

  std::cout << "  All D2Q9 constant checks passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 2. FirstOrder equilibrium: rest state (u = 0)
//    feq(k) = w_k * phi * (1 + c_k.u / cs^2)
//    When u=0: feq(k) = w_k * phi
// ---------------------------------------------------------------------------
static int testFirstOrderRest() {
  printHeader("FirstOrder equilibrium at rest (u=0)");

  Cell cell;
  cell.phi_ = 2.0;                    // phi = 2
  cell.u_   = {0.0, 0.0};            // u = (0, 0)
  // feq(k) = w_k * 2.0

  for (unsigned int k = 0; k < 9; ++k) {
    double feq = equilibrium::FirstOrder<Cell>::get(cell, k);
    double expected = latset::w<LatType>(k) * 2.0;
    ASSERT_NEAR(feq, expected, TOL);
  }

  // Sum of feq should equal phi = 2.0
  double sum = 0.0;
  for (unsigned int k = 0; k < 9; ++k)
    sum += equilibrium::FirstOrder<Cell>::get(cell, k);
  ASSERT_NEAR(sum, 2.0, TOL);

  std::cout << "  FirstOrder rest state passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 3. FirstOrder equilibrium: non-zero velocity
//    feq(k) = w_k * phi * (1 + c_k.u / cs^2)
//    With phi=1.0, u=(0.1, 0.05)
// ---------------------------------------------------------------------------
static int testFirstOrderMoving() {
  printHeader("FirstOrder equilibrium with u=(0.1, 0.05)");

  Cell cell;
  cell.phi_ = 1.0;
  cell.u_   = {0.1, 0.05};

  // c_k dot u  and  3 * (c_k dot u)  [since 1/cs^2 = 3]
  // Actual ordering: c[0]=(0,0), c[1]=(1,0), c[2]=(-1,0), c[3]=(0,1), c[4]=(0,-1),
  //                   c[5]=(1,1), c[6]=(-1,-1), c[7]=(1,-1), c[8]=(-1,1)
  // k=0: cu=0     -> 3*cu=0
  // k=1: cu=0.1   -> 3*cu=0.3
  // k=2: cu=-0.1  -> 3*cu=-0.3
  // k=3: cu=0.05  -> 3*cu=0.15
  // k=4: cu=-0.05 -> 3*cu=-0.15
  // k=5: cu=0.15  -> 3*cu=0.45
  // k=6: cu=-0.15 -> 3*cu=-0.45
  // k=7: cu=0.05  -> 3*cu=0.15
  // k=8: cu=-0.05 -> 3*cu=-0.15
  //
  // feq = w_k * (1 + 3*cu)

  // k=0: w=4/9, 1+0=1         -> feq = 4/9
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 0),
              4.0/9.0, TOL);
  // k=1: w=1/9, 1+0.3=1.3    -> feq = 1.3/9
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 1),
              1.3/9.0, TOL);
  // k=2: w=1/9, 1-0.3=0.7    -> feq = 0.7/9
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 2),
              0.7/9.0, TOL);
  // k=3: w=1/9, 1+0.15=1.15  -> feq = 1.15/9
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 3),
              1.15/9.0, TOL);
  // k=4: w=1/9, 1-0.15=0.85  -> feq = 0.85/9
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 4),
              0.85/9.0, TOL);
  // k=5: w=1/36, 1+0.45=1.45 -> feq = 1.45/36
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 5),
              1.45/36.0, TOL);
  // k=6: w=1/36, 1-0.45=0.55 -> feq = 0.55/36
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 6),
              0.55/36.0, TOL);
  // k=7: w=1/36, 1+0.15=1.15 -> feq = 1.15/36
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 7),
              1.15/36.0, TOL);
  // k=8: w=1/36, 1-0.15=0.85 -> feq = 0.85/36
  ASSERT_NEAR(equilibrium::FirstOrder<Cell>::get(cell, 8),
              0.85/36.0, TOL);

  std::cout << "  FirstOrder moving state passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 4. SecondOrder equilibrium: rest state (u = 0)
//    geq(k) = w_k * rho * (1 + 3(c.u) + 4.5(c.u)^2 - 1.5|u|^2)
//    When u=0: geq(k) = w_k * rho
// ---------------------------------------------------------------------------
static int testSecondOrderRest() {
  printHeader("SecondOrder equilibrium at rest (u=0)");

  Cell cell;
  cell.rho_ = 1.5;
  cell.u_   = {0.0, 0.0};

  for (unsigned int k = 0; k < 9; ++k) {
    double geq = equilibrium::SecondOrder<Cell>::get(cell, k);
    ASSERT_NEAR(geq, latset::w<LatType>(k) * 1.5, TOL);
  }

  // Sum of geq should equal rho = 1.5
  double sum = 0.0;
  for (unsigned int k = 0; k < 9; ++k)
    sum += equilibrium::SecondOrder<Cell>::get(cell, k);
  ASSERT_NEAR(sum, 1.5, TOL);

  std::cout << "  SecondOrder rest state passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 5. SecondOrder equilibrium: non-zero velocity
//    geq(k) = w_k * rho * (1 + 3(c.u) + 4.5(c.u)^2 - 1.5|u|^2)
//    With rho=1.0, u=(0.1, 0.05)
// ---------------------------------------------------------------------------
static int testSecondOrderMoving() {
  printHeader("SecondOrder equilibrium with u=(0.1, 0.05)");

  Cell cell;
  cell.rho_ = 1.0;
  cell.u_   = {0.1, 0.05};

  // Pre-computed values:
  // |u|^2 = 0.01 + 0.0025 = 0.0125
  // 1.5 * |u|^2 = 0.01875
  //
  // For each direction k (actual ordering):
  // geq = w_k * (1 + 3*cu + 4.5*cu^2 - 0.01875)
  //
  // k=0: cu=0,       cu^2=0       -> 1 + 0 + 0 - 0.01875 = 0.98125
  // k=1: cu=0.1,     cu^2=0.01    -> 1 + 0.3 + 0.045 - 0.01875 = 1.32625
  // k=2: cu=-0.1,    cu^2=0.01    -> 1 - 0.3 + 0.045 - 0.01875 = 0.72625
  // k=3: cu=0.05,    cu^2=0.0025  -> 1 + 0.15 + 0.01125 - 0.01875 = 1.1425
  // k=4: cu=-0.05,   cu^2=0.0025  -> 1 - 0.15 + 0.01125 - 0.01875 = 0.8425
  // k=5: cu=0.15,    cu^2=0.0225  -> 1 + 0.45 + 0.10125 - 0.01875 = 1.5325
  // k=6: cu=-0.15,   cu^2=0.0225  -> 1 - 0.45 + 0.10125 - 0.01875 = 0.6325
  // k=7: cu=0.05,    cu^2=0.0025  -> 1 + 0.15 + 0.01125 - 0.01875 = 1.1425
  // k=8: cu=-0.05,   cu^2=0.0025  -> 1 - 0.15 + 0.01125 - 0.01875 = 0.8425

  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 0),
              (4.0/9.0) * 0.98125, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 1),
              (1.0/9.0) * 1.32625, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 2),
              (1.0/9.0) * 0.72625, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 3),
              (1.0/9.0) * 1.1425, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 4),
              (1.0/9.0) * 0.8425, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 5),
              (1.0/36.0) * 1.5325, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 6),
              (1.0/36.0) * 0.6325, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 7),
              (1.0/36.0) * 1.1425, TOL);
  ASSERT_NEAR(equilibrium::SecondOrder<Cell>::get(cell, 8),
              (1.0/36.0) * 0.8425, TOL);

  std::cout << "  SecondOrder moving state passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 6. SecondOrder: macroscopic density conservation
//    sum(geq(k)) over k should equal rho
// ---------------------------------------------------------------------------
static int testSecondOrderDensityConservation() {
  printHeader("SecondOrder density conservation (sum geq = rho)");

  Cell cell;
  cell.rho_ = 1.0;

  // Test with several velocity vectors
  Vector<double, 2> testVelocities[] = {
    {0.0,   0.0},
    {0.1,   0.0},
    {0.0,   0.05},
    {0.1,   0.05},
    {0.2,  -0.15},
    {-0.05, 0.03},
  };

  for (const auto& u : testVelocities) {
    cell.u_ = u;
    double sum = 0.0;
    for (unsigned int k = 0; k < 9; ++k)
      sum += equilibrium::SecondOrder<Cell>::get(cell, k);
    ASSERT_NEAR(sum, 1.0, 1e-14);
  }

  // Also test with rho = 2.5
  cell.rho_ = 2.5;
  cell.u_   = {0.1, 0.05};
  double sum = 0.0;
  for (unsigned int k = 0; k < 9; ++k)
    sum += equilibrium::SecondOrder<Cell>::get(cell, k);
  ASSERT_NEAR(sum, 2.5, 1e-14);

  std::cout << "  Density conservation passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 7. BGK collision without source (hasSource = false)
//    f'_k = f_k - omega * (f_k - feq_k)
// ---------------------------------------------------------------------------
static int testBGKNoSource() {
  printHeader("BGK collision (hasSource=false)");

  Cell cell;
  cell.rho_   = 1.0;
  cell.u_     = {0.1, 0.05};
  cell.omega_ = 1.5;   // tau = 1/omega = 2/3

  // Set pre-collision populations to arbitrary values
  double initPop[9] = {
    0.5, 0.12, 0.10, 0.06, 0.08, 0.035, 0.020, 0.015, 0.030
  };
  for (unsigned int k = 0; k < 9; ++k)
    cell.f_[k] = initPop[k];

  // Compute expected post-collision values
  double feq[9];
  for (unsigned int k = 0; k < 9; ++k)
    feq[k] = equilibrium::SecondOrder<Cell>::get(cell, k);

  double expected[9];
  for (unsigned int k = 0; k < 9; ++k)
    expected[k] = initPop[k] - 1.5 * (initPop[k] - feq[k]);

  // Apply BGK collision (no source)
  collision::BGKSource<equilibrium::SecondOrder<Cell>,
                       void, false>::apply(cell);

  // Verify post-collision populations
  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], expected[k], TOL);
  }

  std::cout << "  BGK no-source collision passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 8. BGK collision with source (hasSource = true)
//    SourceField = CONSTFORCE (which has isField = false)
//    f'_k = f_k - omega*(f_k - feq_k) + (1-0.5*omega)*w_k*(c_k . F)
// ---------------------------------------------------------------------------
static int testBGKWithConstForceSource() {
  printHeader("BGK collision with CONSTFORCE source (hasSource=true)");

  Cell cell;
  cell.rho_        = 1.0;
  cell.u_          = {0.1, 0.05};
  cell.omega_      = 1.0;          // tau = 1, makes arithmetic simpler
  cell.constForce_ = {0.01, 0.02}; // F = (0.01, 0.02)

  // Set pre-collision populations
  double initPop[9] = {
    0.5, 0.12, 0.10, 0.06, 0.08, 0.035, 0.020, 0.015, 0.030
  };
  for (unsigned int k = 0; k < 9; ++k)
    cell.f_[k] = initPop[k];

  // Compute expected values
  double feq[9];
  for (unsigned int k = 0; k < 9; ++k)
    feq[k] = equilibrium::SecondOrder<Cell>::get(cell, k);

  // Source term: (1 - 0.5*omega) * w_k * (c_k . F)
  // With omega=1.0: (1 - 0.5) = 0.5
  double srcFactor = 0.5;
  double expected[9];
  for (unsigned int k = 0; k < 9; ++k) {
    const auto& ck = latset::c<LatType>(k);
    double wk      = latset::w<LatType>(k);
    double cdotF   = ck * cell.constForce_;
    double srcTerm = srcFactor * wk * cdotF;

    // f' = f - omega*(f - feq) + srcTerm
    // With omega=1: f' = feq + srcTerm
    expected[k] = feq[k] + srcTerm;
  }

  // Apply BGK collision with source
  collision::BGKSource<equilibrium::SecondOrder<Cell>,
                       CONSTFORCE<double, 2>, true>::apply(cell);

  // Verify
  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], expected[k], TOL);
  }

  std::cout << "  BGK with CONSTFORCE source passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 9. BGK collision with FORCE source (hasSource=true, isField=true)
// ---------------------------------------------------------------------------
static int testBGKWithForceSource() {
  printHeader("BGK collision with FORCE source (hasSource=true)");

  Cell cell;
  cell.rho_   = 1.0;
  cell.u_     = {0.1, 0.05};
  cell.omega_ = 1.5;               // tau = 2/3
  cell.force_ = {-0.02, 0.03};    // F = (-0.02, 0.03)

  // Set pre-collision populations to initial guess
  double initPop[9] = {
    0.45, 0.14, 0.13, 0.07, 0.09, 0.04, 0.025, 0.018, 0.035
  };
  for (unsigned int k = 0; k < 9; ++k)
    cell.f_[k] = initPop[k];

  // Compute expected values
  double feq[9];
  for (unsigned int k = 0; k < 9; ++k)
    feq[k] = equilibrium::SecondOrder<Cell>::get(cell, k);

  // Source term: (1 - 0.5*omega) * w_k * (c_k . F)
  // With omega=1.5: (1 - 0.75) = 0.25
  double srcFactor = 1.0 - 0.5 * 1.5;   // = 0.25
  double expected[9];
  for (unsigned int k = 0; k < 9; ++k) {
    const auto& ck = latset::c<LatType>(k);
    double wk      = latset::w<LatType>(k);
    double cdotF   = ck * cell.force_;
    double srcTerm = srcFactor * wk * cdotF;

    expected[k] = initPop[k] - 1.5 * (initPop[k] - feq[k]) + srcTerm;
  }

  // Apply
  collision::BGKSource<equilibrium::SecondOrder<Cell>,
                       FORCE<double, 2>, true>::apply(cell);

  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], expected[k], TOL);
  }

  std::cout << "  BGK with FORCE source passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 10. BGK FirstOrder collision (no source, for phase-field applications)
// ---------------------------------------------------------------------------
static int testBGKFirstOrderNoSource() {
  printHeader("BGK collision with FirstOrder equilibrium (hasSource=false)");

  Cell cell;
  cell.phi_   = 0.8;
  cell.u_     = {0.05, 0.025};
  cell.omega_ = 1.8;

  double initPop[9] = {
    0.40, 0.10, 0.08, 0.05, 0.07, 0.025, 0.018, 0.012, 0.022
  };
  for (unsigned int k = 0; k < 9; ++k)
    cell.f_[k] = initPop[k];

  double feq[9];
  for (unsigned int k = 0; k < 9; ++k)
    feq[k] = equilibrium::FirstOrder<Cell>::get(cell, k);

  double expected[9];
  for (unsigned int k = 0; k < 9; ++k)
    expected[k] = initPop[k] - 1.8 * (initPop[k] - feq[k]);

  collision::BGKSource<equilibrium::FirstOrder<Cell>,
                       void, false>::apply(cell);

  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], expected[k], TOL);
  }

  std::cout << "  BGK FirstOrder no-source passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 11. Edge cases: zero omega (no relaxation), omega=2 (over-relaxation),
//     negative phi
// ---------------------------------------------------------------------------
static int testEdgeCases() {
  printHeader("Edge cases");

  // Case 1: omega = 0 -> populations should be unchanged
  {
    Cell cell;
    cell.rho_   = 1.0;
    cell.u_     = {0.0, 0.0};
    cell.omega_ = 0.0;

    for (unsigned int k = 0; k < 9; ++k)
      cell.f_[k] = static_cast<double>(k + 1) * 0.1;

    double before[9];
    for (unsigned int k = 0; k < 9; ++k)
      before[k] = cell[k];

    collision::BGKSource<equilibrium::SecondOrder<Cell>,
                         void, false>::apply(cell);

    for (unsigned int k = 0; k < 9; ++k) {
      ASSERT_NEAR(cell[k], before[k], TOL);
    }
    std::cout << "  omega=0 (no relaxation) passed." << std::endl;
  }

  // Case 2: omega = 2 -> f'_k = 2*feq_k - f_k (full over-relaxation)
  {
    Cell cell;
    cell.rho_   = 1.0;
    cell.u_     = {0.1, 0.05};
    cell.omega_ = 2.0;

    double initPop[9] = {
      0.5, 0.12, 0.10, 0.06, 0.08, 0.035, 0.020, 0.015, 0.030
    };
    for (unsigned int k = 0; k < 9; ++k)
      cell.f_[k] = initPop[k];

    double feq[9];
    for (unsigned int k = 0; k < 9; ++k)
      feq[k] = equilibrium::SecondOrder<Cell>::get(cell, k);

    collision::BGKSource<equilibrium::SecondOrder<Cell>,
                         void, false>::apply(cell);

    for (unsigned int k = 0; k < 9; ++k) {
      double expected = 2.0 * feq[k] - initPop[k];
      ASSERT_NEAR(cell[k], expected, TOL);
    }
    std::cout << "  omega=2 (over-relaxation) passed." << std::endl;
  }

  // Case 3: negative phi value for FirstOrder
  {
    Cell cell;
    cell.phi_   = -0.5;
    cell.u_     = {0.0, 0.0};

    for (unsigned int k = 0; k < 9; ++k) {
      double feq = equilibrium::FirstOrder<Cell>::get(cell, k);
      ASSERT_NEAR(feq, latset::w<LatType>(k) * (-0.5), TOL);
    }
    std::cout << "  negative phi passed." << std::endl;
  }

  return 0;
}

// ---------------------------------------------------------------------------
// 12. Guo force: zero force -> populations unchanged
// ---------------------------------------------------------------------------
static int testGuoForceZero() {
  printHeader("Guo force with F=(0,0)");
  Cell cell;
  cell.u_ = {0.1, 0.05};
  cell.omega_ = 1.5;

  double before[9];
  for (unsigned int k = 0; k < 9; ++k) {
    cell.f_[k] = static_cast<double>(k + 1) * 0.1;
    before[k] = cell.f_[k];
  }
  cell.constForce_ = {0.0, 0.0};

  force::Guo<CONSTFORCE<double, 2>>::apply(cell);

  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], before[k], TOL);
  }
  std::cout << "  zero force passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 13. Guo force: known gravity force F=(0,-g)
// ---------------------------------------------------------------------------
static int testGuoForceGravity() {
  printHeader("Guo force with F=(0, -0.02) [gravity]");
  Cell cell;
  cell.u_ = {0.0, 0.0};          // rest state
  cell.omega_ = 1.0;             // omega=1 -> factor = 0.5
  cell.constForce_ = {0.0, -0.02};

  for (unsigned int k = 0; k < 9; ++k)
    cell.f_[k] = 0.0;            // start from zero

  force::Guo<CONSTFORCE<double, 2>>::apply(cell);

  // factor = 1 - 0.5 = 0.5, u = (0,0)
  // G_k = 0.5 * w_k * (c_k·F / cs² + 0)
  //     = 0.5 * w_k * (c_k·F) * 3
  //     = 1.5 * w_k * (c_k·F)
  //
  // c_k·F = c_k_y * (-0.02):
  // k=0: 0, k=1: 0, k=2: 0, k=3: -0.02, k=4: 0.02,
  // k=5: -0.02, k=6: 0.02, k=7: 0.02, k=8: -0.02
  //
  // Expected:
  // k=0: 0
  // k=1: 0
  // k=2: 0
  // k=3: 1.5 * (1/9) * (-0.02) = -0.02/6 = -0.00333333...
  // k=4: 1.5 * (1/9) * 0.02   = 0.02/6   = 0.00333333...
  // k=5: 1.5 * (1/36) * (-0.02) = -0.03/36
  // ...

  double expected[9] = {0};
  for (unsigned int k = 0; k < 9; ++k) {
    double ck_F = latset::c<LatType>(k) * cell.constForce_;
    expected[k] = 1.5 * latset::w<LatType>(k) * ck_F;
  }

  for (unsigned int k = 0; k < 9; ++k) {
    ASSERT_NEAR(cell[k], expected[k], TOL);
  }
  std::cout << "  gravity force passed." << std::endl;
  return 0;
}

// ---------------------------------------------------------------------------
// 14. Guo force: density conservation (sum G_k = 0)
// ---------------------------------------------------------------------------
static int testGuoForceDensityConservation() {
  printHeader("Guo force density conservation (sum G_k = 0)");
  Cell cell;
  cell.omega_ = 1.5;
  cell.constForce_ = {0.01, -0.03};

  Vector<double, 2> testVels[] = {
    {0.0, 0.0}, {0.1, 0.05}, {-0.05, 0.03}, {0.2, -0.15}, {0.001, -0.002}
  };

  for (const auto& u : testVels) {
    cell.u_ = u;
    for (unsigned int k = 0; k < 9; ++k)
      cell.f_[k] = static_cast<double>(k) * 0.01;  // arbitrary init

    double sumBefore = 0;
    for (unsigned int k = 0; k < 9; ++k)
      sumBefore += cell[k];

    force::Guo<CONSTFORCE<double, 2>>::apply(cell);

    double sumAfter = 0;
    for (unsigned int k = 0; k < 9; ++k)
      sumAfter += cell[k];

    ASSERT_NEAR(sumAfter, sumBefore, 1e-14);
  }
  std::cout << "  density conservation passed." << std::endl;
  return 0;
}

// ===========================================================================
// main
// ===========================================================================
int main() {
  std::cout << std::setprecision(16);
  std::cout << "========================================================\n";
  std::cout << "  collision.h Unit Tests (D2Q9)\n";
  std::cout << "========================================================"
            << std::endl;

  failures = 0;

  testD2Q9Constants();
  testFirstOrderRest();
  testFirstOrderMoving();
  testSecondOrderRest();
  testSecondOrderMoving();
  testSecondOrderDensityConservation();
  testBGKNoSource();
  testBGKWithConstForceSource();
  testBGKWithForceSource();
  testBGKFirstOrderNoSource();
  testEdgeCases();
  testGuoForceZero();
  testGuoForceGravity();
  testGuoForceDensityConservation();

  std::cout << "\n========================================================"
            << std::endl;
  if (failures == 0) {
    std::cout << "  ALL TESTS PASSED" << std::endl;
  } else {
    std::cout << "  " << failures << " TEST(S) FAILED" << std::endl;
  }
  std::cout << "========================================================"
            << std::endl;

  return failures;
}
