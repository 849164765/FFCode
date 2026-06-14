// phase_gradient.ur.h
// Phase A: Gradient and Laplacian computation for phase-field LBM
// Eqs. 34 & 35 from Fakhari et al., Phys. Rev. E 96, 053301 (2017)
//
// ∇φ  = (1/cs²) * Σ_α e_α w_α φ(x+e_α)         — Eq.34, isotropic central diff
// ∇²φ = (2/cs²) * Σ_α w_α [φ(x+e_α) - φ(x)]    — Eq.35, isotropic Laplacian

#pragma once

#include "data_struct/cell.h"
#include "lbm/lattice_set.h"
#include "ff/ff2d.h"  // for LAPLACIAN field type

namespace phase_gradient {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;
#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// ===================================================================
// ComputeGradientPhi — Eq.34 isotropic central difference gradient
// Reads:  PHI<T> (phase field)
// Writes: GRAD<T, D> (gradient vector per cell)
// ===================================================================

template <typename CELLTYPE>
struct ComputeGradientPhi {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // Full-cell apply: reads neighbor phi, writes cell gradient
  __any__ static inline void apply(CELL& cell) {
    // Accumulate gradient on the stack; write at the end for immutability
    Vector<T, LatSet::d> grad;
    for (unsigned int d = 0; d < LatSet::d; ++d) grad[d] = T{0};

    // Skip α=0 because e_0=(0,0) contributes zero to gradient
    for (unsigned int alpha = 1; alpha < LatSet::q; ++alpha) {
      const T phi_neighbor =
        cell.getNeighbor(alpha).template get<PHI<T>>();
      const T w_alpha = latset::w<LatSet>(alpha);
      const auto& c_alpha = latset::c<LatSet>(alpha);
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        grad[d] += w_alpha * c_alpha[d] * phi_neighbor;
      }
    }
    // Eq.34 prefactor:  c/(cs²·δx) = 1/cs²  (lattice units: δx=c=1)
    const T inv_cs2 = T{1} / LatSet::cs2;
    for (unsigned int d = 0; d < LatSet::d; ++d) grad[d] *= inv_cs2;

    cell.template get<GRAD<T, LatSet::d>>() = grad;
  }
};


// ===================================================================
// ComputeLaplacianPhi — Eq.35 isotropic Laplacian
// Reads:  PHI<T> (phase field)
// Writes: LAPLACIAN<T> (scalar Laplacian per cell)
// ===================================================================

template <typename CELLTYPE>
struct ComputeLaplacianPhi {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T phi_self = cell.template get<PHI<T>>();
    T laplacian = T{0};

    // Sum over all directions (including α=0 for consistency with Eq.35)
    // w_0 * [φ(x) - φ(x)] = 0, so α=0 contributes zero in practice
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const T phi_neighbor =
        cell.getNeighbor(alpha).template get<PHI<T>>();
      const T w_alpha = latset::w<LatSet>(alpha);
      laplacian += w_alpha * (phi_neighbor - phi_self);
    }
    // Eq.35 prefactor:  2c²/(cs²·δx²) = 2/cs²  (lattice units)
    const T prefactor = T{2} / LatSet::cs2;
    laplacian *= prefactor;

    cell.template get<ff::LAPLACIAN<T>>() = laplacian;
  }
};

}  // namespace phase_gradient


// ===================================================================
// ComputeNormalFromGradient — compute unit normal from gradient
// Reads:  GRAD<T, D>
// Writes: NORMAL<T, D> = ∇φ/|∇φ|  (zero vector if |∇φ| < ε)
// Required by Eq.7 (F^φ_α depends on e_α · n)
// ===================================================================

namespace phase_gradient {

template <typename CELLTYPE>
struct ComputeNormalFromGradient {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const Vector<T, LatSet::d>& grad = cell.template get<GRAD<T, LatSet::d>>();
    Vector<T, LatSet::d> normal;
    T grad_mag = grad.getnorm();
    // Threshold to suppress noise in bulk (same as FF2D: 0.005)
    T epsilon = T{0.005};
    if (grad_mag < epsilon) {
      for (unsigned int d = 0; d < LatSet::d; ++d) normal[d] = T{0};
    } else {
      for (unsigned int d = 0; d < LatSet::d; ++d) normal[d] = grad[d] / grad_mag;
    }
    cell.template get<NORMAL<T, LatSet::d>>() = normal;
  }
};

}  // namespace phase_gradient
