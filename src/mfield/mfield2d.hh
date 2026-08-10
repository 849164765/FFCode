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

  T phi   = pf_cell.template get<PHI<T>>();
  T mu_l  = mf_cell.template get<MU_L<T>>();
  T mu_h  = mf_cell.template get<MU_H<T>>();
  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();

  T mu  = mu_l  + phi * (mu_h  - mu_l);
  T chi = chi_l + phi * (chi_h - chi_l);

  // tau = 0.5 + PsiSolver_K * mu  (diffusion D = cs²·(tau-0.5) = cs²·K·mu ∝ mu)
  // Any K leaves the continuum ∇·(μ∇ψ)=0 fixed point unchanged; K=1/cs² is the
  // classical choice. K=0.5 (omega(μ=9)=0.2, omega(μ=1)=1.0) boosts the
  // relaxation so the discrete fixed point is smooth (grid-scale |H|
  // oscillation 0.040→0.008·H0, numerical spike 2.58→1.91; the residual ~1.9
  // is the physical pole concentration 2μ_f/(μ_f+μ_b)) and the sub-iterated
  // solve converges in ~100 iterations — see .scratch/psi_boost_experiment.py.
  T PsiSolver_K = mf_cell.template get<PSI_K<T>>();
  T tau_psi = T{0.5} + PsiSolver_K * mu;
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
// Magnetic body force (Guo et al., Phys. Fluids 37, 022148 (2025), Eq. (8)):
//   F_m = (μ₀χ/2) ∇(|H|²)   (with μ₀ ≡ 1)
//
// This form assumes linear magnetization (χ constant within each phase,
// interpolated by φ across the interface). The ∇χ (magnetostrictive) term
// is NOT included, following the reference.
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5 for neighbor access

  T chi{};
  if constexpr (mf_cell.template hasField<CHI_PERCELL<T>>()) {
    chi = mf_cell.template get<CHI_PERCELL<T>>();
  } else {
    T chi_l = mf_cell.template get<CHI_L<T>>();
    T chi_h = mf_cell.template get<CHI_H<T>>();
    T phi   = pf_cell.template get<typename PFCELL::GenericRho>();
    chi = chi_l + phi * (chi_h - chi_l);
  }

  // Compute ∇(|H|²) directly. It is analytically equivalent to
  // 2|H|∇|H|, but is less sensitive to product-rule errors on a diffuse
  // interface when evaluated with the same centered stencil.
  Vector<T, 2> grad_H2{0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T H2_k = Hmag_k * Hmag_k;
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_H2[0] += wk * ck[0] * H2_k;
    grad_H2[1] += wk * ck[1] * H2_k;
  }
  grad_H2[0] /= LatSet::cs2;
  grad_H2[1] /= LatSet::cs2;

  T coeff = T{0.5} * chi;
  T Fmag_x = coeff * grad_H2[0];
  T Fmag_y = coeff * grad_H2[1];

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

}  // namespace mfield
