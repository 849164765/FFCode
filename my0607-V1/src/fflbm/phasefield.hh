// phasefield.hh — order parameter pre-processing (paper-aligned)
//
// PhiUpdate<CELL>            : φ = Σg_i, clamp [0,1]
// PhiGradient<CELL>          : ∇φ = (1/cs²)·Σ w_i·c_i·φ_neighbor
// PhiLaplacian<CELL>         : ∇²φ = (2/cs²)·Σ w_i·(φ_neighbor - φ_self)
// ChemPotential<CELL>        : μ = 4β·φ(φ-1)(φ-0.5) - κ·∇²φ
// ChemPotentialGradient<CELL>: ∇μ = (1/cs²)·Σ w_i·c_i·μ_neighbor

#pragma once

#include "ff/ff2d.h"

namespace fflbm {

template <typename CELL>
struct PhiUpdate {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    T phi = T{0};
    for (unsigned int k = 0; k < LatSet::q; ++k) phi += cell[k];
    if (phi < T{0}) phi = T{0};
    if (phi > T{1}) phi = T{1};
    cell.template get<GenericRho>() = phi;
  }
};

template <typename CELL>
struct PhiGradient {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    Vector<T, LatSet::d> grad;
    grad[0] = T{0};
    grad[1] = T{0};
    for (unsigned int i = 1; i < LatSet::q; ++i) {
      T phi_i = cell.getNeighbor(i).template get<typename CELL::GenericRho>();
      T wi = latset::w<LatSet>(i);
      const auto& ci = latset::c<LatSet>(i);
      grad[0] += wi * ci[0] * phi_i;
      grad[1] += wi * ci[1] * phi_i;
    }
    grad[0] /= LatSet::cs2;
    grad[1] /= LatSet::cs2;
    cell.template get<GRAD<T, LatSet::d>>() = grad;

    T grad_mag = grad.getnorm();
    T epsilon{T(0.005)};
    Vector<T, LatSet::d> n;
    if (grad_mag < epsilon) {
      n[0] = T{0};
      n[1] = T{0};
    } else {
      n[0] = grad[0] / grad_mag;
      n[1] = grad[1] / grad_mag;
    }
    cell.template get<NORMAL<T, LatSet::d>>() = n;
  }
};

template <typename CELL>
struct PhiLaplacian {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    T phi_self = cell.template get<typename CELL::GenericRho>();
    T laplacian = T{0};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T phi_k = cell.getNeighbor(k).template get<typename CELL::GenericRho>();
      T wk = latset::w<LatSet>(k);
      laplacian += wk * (phi_k - phi_self);
    }
    laplacian *= T{2} / LatSet::cs2;
    cell.template get<ff::LAPLACIAN<T>>() = laplacian;
  }
};

template <typename CELL>
struct ChemPotential {
  using T = typename CELL::FloatType;

  __any__ static void apply(CELL& cell) {
    T phi = cell.template get<typename CELL::GenericRho>();
    T laplacian = cell.template get<ff::LAPLACIAN<T>>();
    T beta = cell.template get<ff::BETA<T>>();
    T kappa = cell.template get<ff::KAPPA<T>>();
    T double_well = phi * (phi - T{1}) * (phi - T{0.5});
    cell.template get<ff::CHEMICALPOTENTIAL<T>>() =
        T{4} * beta * double_well - kappa * laplacian;
  }
};

// ∇μ = (1/cs²)·Σ_{i=1}^{q-1} w_i·c_i·μ(x+c_i)
// Writes to NORMAL field — consumed by BGKSourceCollision as source direction
template <typename CELL>
struct ChemPotentialGradient {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    Vector<T, LatSet::d> grad_mu;
    grad_mu[0] = T{0};
    grad_mu[1] = T{0};
    for (unsigned int i = 1; i < LatSet::q; ++i) {
      T mu_i = cell.getNeighbor(i).template get<ff::CHEMICALPOTENTIAL<T>>();
      T wi = latset::w<LatSet>(i);
      const auto& ci = latset::c<LatSet>(i);
      grad_mu[0] += wi * ci[0] * mu_i;
      grad_mu[1] += wi * ci[1] * mu_i;
    }
    grad_mu[0] /= LatSet::cs2;
    grad_mu[1] /= LatSet::cs2;
    cell.template get<NORMAL<T, LatSet::d>>() = grad_mu;
  }
};

}  // namespace fflbm
