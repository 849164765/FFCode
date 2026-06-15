// magnetic_field.h
// Magnetic field post-processing operators (D2Q5)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   Eq.4:  MagGradient  (H = -grad(psi))
//   Eq.8:  MagHSq       (|H|^2 / 2)
//   Eq.8:  MagForce     (F_m = chi * grad(|H|^2/2))

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// MagGradient (Eq.4): H = -grad(psi)
// H = -(1/cs^2) * Σ_{α≠0} w_α * c_α * psi(x + c_α*dt)
// Reads:  PSI (neighbors)
// Writes: MAGH
// ================================================================
template <typename CELLTYPE>
struct MagGradient {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    Vector<T, LatSet::d> H{};
    for (int k = 1; k < LatSet::q; ++k) {
      T psi_neighbor =
          cell.getNeighbor(k).template get<PSI<T>>();
      H += latset::w<LatSet>(k) * latset::c<LatSet>(k) * psi_neighbor;
    }
    // H = -grad(psi)
    cell.template get<MAGH<T, LatSet::d>>() = H * (-LatSet::InvCs2);
  }
};

// ================================================================
// MagHSq (Eq.8): |H|^2 / 2
// Reads:  MAGH
// Writes: MAGHSQ
// ================================================================
template <typename CELLTYPE>
struct MagHSq {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const auto& H = cell.template get<MAGH<T, LatSet::d>>();
    cell.template get<MAGHSQ<T>>() = (H * H) * T{0.5};
  }
};

// ================================================================
// MagForce (Eq.8): F_m = chi * grad(|H|^2 / 2)
// grad(HSQ) via isotropic gradient of MAGHSQ
// chi = mu - 1  (mu0 = 1 in lattice units)
// Reads:  MAGHSQ (neighbors), MAGPERMEABILITY
// Writes: FORCE (on magnetic cell; T14 coupling copies to NS cell)
// ================================================================
template <typename CELLTYPE>
struct MagForce {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    // Gradient of |H|^2/2  (isotropic, excluding zero-velocity direction)
    Vector<T, LatSet::d> gradHSQ{};
    for (int k = 1; k < LatSet::q; ++k) {
      T hsq_neighbor =
          cell.getNeighbor(k).template get<MAGHSQ<T>>();
      gradHSQ += latset::w<LatSet>(k) * latset::c<LatSet>(k) * hsq_neighbor;
    }
    gradHSQ = gradHSQ * LatSet::InvCs2;

    // Susceptibility: chi = mu - 1  (mu0 = 1 in lattice units)
    T mu = cell.template get<MAGPERMEABILITY<T>>();
    T chi = mu - T{1};

    // F_m = chi * grad(|H|^2/2)
    cell.template get<FORCE<T, LatSet::d>>() = chi * gradHSQ;
  }
};
