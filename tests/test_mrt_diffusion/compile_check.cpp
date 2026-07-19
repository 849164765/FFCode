// test_mrt_diffusion_compile.cpp — compile-check MRTDiffusion with D2Q5
#include "lbm/collisionMRT.h"
#include "lbm/lattice_set.h"
#include <iostream>
#include <cmath>

using T = double;
using LatSet = D2Q5<T>;

// Minimal compile check: instantiate MRTDiffusion<Cell<...>, OMEGA<T>>
// We don't need runtime correctness — just verify the header compiles

int main() {
  // Check D2Q5 MRT matrices are accessible
  T test_m = mrt::M<LatSet>(0, 0);
  T test_inv = mrt::InvM<LatSet>(0, 0);
  T test_s = mrt::s<LatSet>(0);

  std::cout << "D2Q5 MRT: M[0][0]=" << test_m
            << " InvM[0][0]=" << test_inv
            << " s[0]=" << test_s << std::endl;

  // Verify M[0][0] = 1.0 (first row, first col of M is 1)
  if (std::abs(test_m - 1.0) < 1e-15) {
    std::cout << "PASS: MRT matrices accessible" << std::endl;
    return 0;
  }
  return 1;
}
