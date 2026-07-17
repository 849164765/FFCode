#pragma once

#include "mf/mf2d.h"
#include "lbm/collisionMRT.h"

namespace mf {

// ---- MFCollision: D2Q5 MRT collision for ∇·(μ∇ψ)=0 ----
// Uses per-cell omega_m(x) = 1/(0.5 + μ(x)/cs²) for variable permeability.
// ω_m is read from OMEGA<T> field (per-cell) if available, fallback to block-level.
// Solves ∂ψ/∂t = ∇·(μ∇ψ) → steady state gives ∇·(μ∇ψ) = 0
// Equilibrium: m_eq = [ψ, 0, 0, 0, 0]
// Relaxation: s = [0, ω_m(x), ω_m(x), ω_m(x), ω_m(x)]
template <typename CELL>
__any__ void MFCollision<CELL>::apply(CELL& cell) {
  T psi = T{0};
  for (unsigned int k = 0; k < LatSet::q; ++k) psi += cell[k];
  cell.template get<typename CELL::GenericRho>() = psi;

  T omega_m{};
  if constexpr (cell.template hasField<OMEGA<T>>()) {
    omega_m = cell.template get<OMEGA<T>>();
  } else {
    omega_m = cell.getOmega();
  }

  T h[5];
  for (unsigned int k = 0; k < 5; ++k) h[k] = cell[k];

  T m[5];
  m[0] = h[0] + h[1] + h[2] + h[3] + h[4];
  m[1] = h[1] - h[2];
  m[2] = h[3] - h[4];
  m[3] = T{-4} * h[0] + h[1] + h[2] + h[3] + h[4];
  m[4] = h[1] + h[2] - h[3] - h[4];

  T m_eq[5] = {psi, T{0}, T{0}, T{0}, T{0}};

  T rtvec[5] = {T{0}, omega_m, omega_m, omega_m, omega_m};

  for (unsigned int i = 0; i < 5; ++i) {
    T coll{};
    for (unsigned int j = 0; j < 5; ++j) {
      coll += mrt::InvM<LatSet>(i, j) * rtvec[j] * (m[j] - m_eq[j]);
    }
    cell[i] = h[i] - coll;
  }
}

// ---- MFMacro: ψ = Σh_α ----
template <typename CELL>
__any__ void MFMacro<CELL>::apply(CELL& cell) {
  T psi = T{0};
  for (unsigned int k = 0; k < LatSet::q; ++k) psi += cell[k];
  cell.template get<GenericRho>() = psi;
}

// ---- MFInitPsi: ψ = -H₀·y  (homogeneous solution for H_y = H₀) ----
// This is NOT an ApplyInnerCellDynamics functor — it's called from the
// main loop to set initial ψ and equilibrium populations on the MF lattice.
// After init, pseudo-time iterations converge the population field to the
// steady-state solution with Neumann BCs.

// ---- MFHField: H = -∇ψ via 4-point central difference ----
// Hx = -∂ψ/∂x ≈ -(ψ(i+1,j) - ψ(i-1,j)) / 2
// Hy = -∂ψ/∂y ≈ -(ψ(i,j+1) - ψ(i,j-1)) / 2
// |H|² written to H2<T>, H vector written to GRAD<T,2> (reuse field)
template <typename MFCELL>
__any__ void MFHField<MFCELL>::apply(MFCELL& cell) {
  T psi_E = cell.getNeighbor(1).template get<typename MFCELL::GenericRho>();
  T psi_W = cell.getNeighbor(2).template get<typename MFCELL::GenericRho>();
  T psi_N = cell.getNeighbor(3).template get<typename MFCELL::GenericRho>();
  T psi_S = cell.getNeighbor(4).template get<typename MFCELL::GenericRho>();

  T dpsi_dx = (psi_E - psi_W) / T{2};
  T dpsi_dy = (psi_N - psi_S) / T{2};

  T Hx = -dpsi_dx;
  T Hy = -dpsi_dy;

  T H2_val = Hx * Hx + Hy * Hy;
  cell.template get<H2<T>>() = H2_val;
  cell.template get<GRAD<T, LatSet::d>>()[0] = Hx;
  cell.template get<GRAD<T, LatSet::d>>()[1] = Hy;
}

// ---- MFGradH2: ∇(|H|²) via 4-point central difference ----
// ∂(H²)/∂x ≈ (H²(i+1,j) - H²(i-1,j)) / 2
// ∂(H²)/∂y ≈ (H²(i,j+1) - H²(i,j-1)) / 2
template <typename MFCELL>
__any__ void MFGradH2<MFCELL>::apply(MFCELL& cell) {
  T H2_E = cell.getNeighbor(1).template get<H2<T>>();
  T H2_W = cell.getNeighbor(2).template get<H2<T>>();
  T H2_N = cell.getNeighbor(3).template get<H2<T>>();
  T H2_S = cell.getNeighbor(4).template get<H2<T>>();

  auto& gradH2 = cell.template get<GRAD_H2<T>>();
  gradH2[0] = (H2_E - H2_W) / T{2};
  gradH2[1] = (H2_N - H2_S) / T{2};
}

// ---- MFChiUpdate: PF → MF coupling ----
// χ_cell = χ_ferro · φ
// PFCELL provides phi (via GenericRho/PHI), MFCELL receives per-cell χ
template <typename PFCELL, typename MFCELL>
__any__ void MFChiUpdate<PFCELL, MFCELL>::apply(PFCELL& pf_cell, MFCELL& mf_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T chi_ferro = mf_cell.template get<CHI_FERRO<T>>();
  mf_cell.template get<CHI_CELL<T>>() = chi_ferro * phi;
}

// ---- MFKelvinForce: MF → NS coupling ----
// F_m = (μ₀·χ_cell/2)·∇(|H|²)
// MFCELL provides per-cell χ and ∇(|H|²), NSCELL receives force
template <typename MFCELL, typename NSCELL>
__any__ void MFKelvinForce<MFCELL, NSCELL>::apply(MFCELL& mf_cell, NSCELL& ns_cell) {
  T chi_cell = mf_cell.template get<CHI_CELL<T>>();
  T mu0 = mf_cell.template get<MU0<T>>();
  const auto& gradH2 = mf_cell.template get<GRAD_H2<T>>();

  T coeff = mu0 * chi_cell / T{2};

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += coeff * gradH2[0];
  ns_force[1] += coeff * gradH2[1];
}

}  // namespace mf
