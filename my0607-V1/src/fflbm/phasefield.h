// phasefield.h
// Phase-field pre-computation operators
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   Eq.10: ComputeNormal       (n = grad_phi / |grad_phi|)
//   Eq.13: ChemicalPotential   (lambda = 4*beta*phi*(phi-1)*(phi-0.5) - kappa*laplacian)

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// ComputeNormal (Eq.10): n = grad_phi / |grad_phi|
// Reads:  GRAD
// Writes: PHASENORMAL
// ================================================================
template <typename CELLTYPE>
struct ComputeNormal {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const auto& grad = cell.template get<GRAD<T, LatSet::d>>();
    T norm = grad.getnorm();
    if (norm > T{0}) {
      cell.template get<PHASENORMAL<T, LatSet::d>>() = grad / norm;
    } else {
      cell.template get<PHASENORMAL<T, LatSet::d>>().clear();
    }
  }
};

// ================================================================
// ChemicalPotential (Eq.13): lambda = 4*beta*phi*(phi-1)*(phi-0.5) - kappa*laplacian
// where kappa = beta * W^2 / 8  (from Eq.14: beta=12*sigma/W, kappa=3*W*sigma/2)
// Reads:  PHI, PHASELAPLACIAN, GBETA, INTERFACEWIDTH
// Writes: PHASECHEMPOTENTIAL
// ================================================================
template <typename CELLTYPE>
struct ChemicalPotential {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T phi = cell.template get<PHI<T>>();
    T lap = cell.template get<PHASELAPLACIAN<T>>();
    T beta = cell.template get<GBETA<T>>();
    T W = cell.template get<INTERFACEWIDTH<T>>();

    // kappa = beta * W^2 / 8  (Eq.14)
    T kappa = beta * W * W / T{8};

    // lambda = 4*beta*phi*(phi-1)*(phi-0.5) - kappa*lap
    T lambda = T{4} * beta * phi * (phi - T{1}) * (phi - T{0.5}) - kappa * lap;

    cell.template get<PHASECHEMPOTENTIAL<T>>() = lambda;
  }
};
