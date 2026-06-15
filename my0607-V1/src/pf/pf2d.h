// pf2d.h
// Phase field physics functors for Allen-Cahn two-phase flow (2D, D2Q9)
//
// Each functor implements `static void apply(CELL&)` for TaskSelector dispatch.
//
// Allen-Cahn equation:
//   ∂φ/∂t + ∇·(φu) = M * ∇²μ
//   μ = 4β φ(φ-1)(φ-0.5) - κ ∇²φ

#pragma once

#include "pf/pf2d.hh"
#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include <cmath>

namespace pf {

// =========================================================================
// PFGradient2D<CELL>
//   Isotropic D2Q9 finite-difference gradient:
//     grad(φ) = InvCs2 * Σ_i w_i * c_i * φ(x + e_i)
//   Reads:  PHI<T> (from cell and all neighbors)
//   Writes: GRAD<T, D>
// =========================================================================
template <typename CELL>
struct PFGradient2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int D = LatSet::d;
  static constexpr unsigned int Q = LatSet::q;

  static void apply(CELL& cell) {
    Vector<T, D> grad{};
    for (unsigned int k = 0; k < Q; ++k) {
      auto nbr = cell.getNeighbor(k);
      T phi_nbr = nbr.template get<PHI<T>>();
      grad += (latset::w<LatSet>(k) * phi_nbr) * latset::c<LatSet>(k);
    }
    grad = grad * LatSet::InvCs2;
    cell.template get<GRAD<T, D>>() = grad;
  }
};

// =========================================================================
// PFLaplacian2D<CELL>
//   D2Q9 Laplacian:
//     ∇²φ = 2 * InvCs2 * Σ_i w_i * (φ(x+e_i) - φ(x))
//   Returns T (not stored in any field — used by PFChemPotential2D).
// =========================================================================
template <typename CELL>
struct PFLaplacian2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int Q = LatSet::q;

  static T compute(CELL& cell) {
    T lap = T(0);
    T phi_c = cell.template get<PHI<T>>();
    for (unsigned int k = 0; k < Q; ++k) {
      auto nbr = cell.getNeighbor(k);
      T phi_n = nbr.template get<PHI<T>>();
      lap += latset::w<LatSet>(k) * (phi_n - phi_c);
    }
    lap = lap * (T(2) * LatSet::InvCs2);
    return lap;
  }
};

// =========================================================================
// PFChemPotential2D<CELL>
//   Allen-Cahn chemical potential:
//     μ = 4β φ(φ-1)(φ-0.5) - κ ∇²φ
//   Reads:  PHI<T>
//   Writes: INTERFACEWIDTH<T> (as temporary storage for μ)
//   Requires: pf::Beta<T> and pf::Kappa<T> set globally
// =========================================================================
template <typename CELL>
struct PFChemPotential2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static void apply(CELL& cell) {
    T phi = cell.template get<PHI<T>>();
    T lap = PFLaplacian2D<CELL>::compute(cell);
    T beta = pf::Beta<T>;
    T kappa = pf::Kappa<T>;
    T mu = T(4) * beta * phi * (phi - T(1)) * (phi - T(0.5)) - kappa * lap;
    cell.template get<INTERFACEWIDTH<T>>() = mu;
  }
};

// =========================================================================
// PFNormal2D<CELL>
//   Interface normal vector from gradient:
//     n = grad(φ) / |grad(φ)|
//   Reads:  GRAD<T, D>
//   Writes: NORMAL<T, D>
// =========================================================================
template <typename CELL>
struct PFNormal2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int D = LatSet::d;

  static void apply(CELL& cell) {
    const auto& grad = cell.template get<GRAD<T, D>>();
    T mag = std::sqrt(grad * grad);
    const T eps = T(1e-12);
    cell.template get<NORMAL<T, D>>() = (mag > eps)
                                        ? grad * (T(1) / mag)
                                        : Vector<T, D>{};
  }
};

}  // namespace pf
