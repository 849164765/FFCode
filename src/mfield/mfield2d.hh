#pragma once

#include "mfield/mfield2d.h"

namespace mfield {

// ---- MFUpdateCoeffs2D ----
// Update per-cell mu, chi, omega_psi from phi
// mu  = mu_l  + phi * (mu_h  - mu_l)
// chi = chi_l + phi * (chi_h - chi_l)
// tau_psi = 0.5 + mu / cs²,  omega_psi = 1 / tau_psi (clamped)
template <typename PFCELL, typename MFCELL>
__any__ void MFUpdateCoeffs2D<PFCELL, MFCELL>::apply(PFCELL& pf_cell,
                                                     MFCELL& mf_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;

  T phi   = pf_cell.template get<PHI<T>>();
  T mu_l  = mf_cell.template get<MU_L<T>>();
  T mu_h  = mf_cell.template get<MU_H<T>>();
  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();

  T mu  = mu_l  + phi * (mu_h  - mu_l);
  T chi = chi_l + phi * (chi_h - chi_l);

  // tau = 0.5 + mu / cs²  (diffusion: D = mu,  tau = 0.5 + D/cs²)
  T tau_psi = T{0.5} + mu / LatSet::cs2;
  T omega_psi = T{1} / tau_psi;
  if (omega_psi > T{1.95}) omega_psi = T{1.95};
  if (omega_psi < T{0.01}) omega_psi = T{0.01};

  mf_cell.template get<MU_PERCELL<T>>()   = mu;
  mf_cell.template get<CHI_PERCELL<T>>()  = chi;
  mf_cell.template get<OMEGA_PSI<T>>()    = omega_psi;
}

// ---- MFComputeH2D ----
// H = -∇ψ using D2Q5 4-direction isotropic gradient
// ∇ψ = Σ w_k * c_k * ψ_k / cs²
// D2Q5: c_k ∈ {(0,0),(1,0),(-1,0),(0,1),(0,-1)}
template <typename CELL>
__any__ void MFComputeH2D<CELL>::apply(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  T grad_x = T{0};
  T grad_y = T{0};

  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T psi_k = cell.getNeighbor(k).template get<PSI<T>>();
    T wk    = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_x += wk * ck[0] * psi_k;
    grad_y += wk * ck[1] * psi_k;
  }
  grad_x /= LatSet::cs2;
  grad_y /= LatSet::cs2;

  T Hx = -grad_x;
  T Hy = -grad_y;
  T Hmag = std::sqrt(Hx * Hx + Hy * Hy);

  cell.template get<HX<T>>()   = Hx;
  cell.template get<HY<T>>()   = Hy;
  cell.template get<HMAG<T>>() = Hmag;
}

// ---- MFMagneticForce2D ----
// Kelvin magnetic body force (Guo et al., Phys. Fluids 37, 022148 (2025), Eq. (8)):
//   F_m = (μ₀χ/2) ∇(|H|²)   (with μ₀ ≡ 1)
//
// ∇(|H|²) is computed from the MF lattice (D2Q5 isotropic gradient of HMAG²).
// This is PREFERRED over the solenoidal decomposition F = -(|H|²/2)∇χ because:
//   - The solenoidal form uses ∇φ from the PF lattice (D2Q9), which has a 7.6%
//     gradient error for W=4 tanh profiles. This error is applied to BOTH the
//     magnetic driving force AND the surface tension restoring force, doubling
//     the effective error and pushing H0=8.2 below the instability threshold.
//   - The Kelvin form computes ∇(|H|²) from the MF field (smooth Laplace
//     solution), so the gradient error only affects surface tension, halving
//     the effective error and restoring growth at H0=8.2.
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  (void)pf_cell;  // Kelvin force does not use PF fields

  T chi = mf_cell.template get<CHI_PERCELL<T>>();

  // Compute ∇(|H|²) using D2Q5 isotropic gradient from MF field
  T grad_H2_x = T{0};
  T grad_H2_y = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T H2_k = Hmag_k * Hmag_k;
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_H2_x += wk * ck[0] * H2_k;
    grad_H2_y += wk * ck[1] * H2_k;
  }
  grad_H2_x /= LatSet::cs2;
  grad_H2_y /= LatSet::cs2;

  // F = (χ/2) ∇(|H|²)
  T coeff = T{0.5} * chi;
  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += coeff * grad_H2_x;
  ns_force[1] += coeff * grad_H2_y;
}

}  // namespace mfield
