// gradient.h
// Isotropic gradient and Laplacian operators (Eq.52-53)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//
// These operators sum over α = 1..q-1 (excluding zero-velocity direction)
// and require neighbor cell access via cell.getNeighbor(k).

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// IsotropicGradient (Eq.52)
// ∇φ = (1/cs^2) * Σ_{α≠0} w_α * c_α * φ(x + c_α*dt)
// Reads:  PHI (self + neighbors)
// Writes: GRAD
// ================================================================
template <typename CELLTYPE>
struct IsotropicGradient {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    Vector<T, LatSet::d> grad{};
    for (int k = 1; k < LatSet::q; ++k) {
      T phi_neighbor =
          cell.getNeighbor(k).template get<PHI<T>>();
      grad += latset::w<LatSet>(k) * latset::c<LatSet>(k) * phi_neighbor;
    }
    cell.template get<GRAD<T, LatSet::d>>() = grad * LatSet::InvCs2;
  }
};

// ================================================================
// IsotropicLaplacian (Eq.53)
// ∇²φ = (2/cs^2) * Σ_{α≠0} w_α * [φ(x + c_α*dt) - φ(x)]
// Reads:  PHI (self + neighbors)
// Writes: PHASELAPLACIAN
// ================================================================
template <typename CELLTYPE>
struct IsotropicLaplacian {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T phi_self = cell.template get<PHI<T>>();
    T lap = T{0};
    for (int k = 1; k < LatSet::q; ++k) {
      T phi_neighbor =
          cell.getNeighbor(k).template get<PHI<T>>();
      lap += latset::w<LatSet>(k) * (phi_neighbor - phi_self);
    }
    cell.template get<PHASELAPLACIAN<T>>() = lap * T{2} * LatSet::InvCs2;
  }
};
