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
// Complete Korteweg-Helmholtz magnetic body force:
//   F_mag = χ(φ)·|H|·∇|H| - ½|H|²·∇χ  (with μ₀ ≡ 1)
//   ≡ divergence of Maxwell stress tensor: ∇·(μ H⊗H - ½ μ |H|² I)
//
// Reference: Guo et al., Phys. Fluids 37, 022148 (2025), Eqs. (64-66)
//
// The ∇χ term (magnetostrictive pressure gradient) is essential for correct
// bubble elongation along the field direction when μ_h > μ_l:
//   - At the poles: ∇χ points from bubble → fluid, and -½ H² ∇χ points
//     fluid → bubble (pushing the interface OUTWARD → stretching)
//   - At the equator: ∇|H| points toward the bubble, and χ|H|∇|H|
//     points fluid → bubble (pushing the interface INWARD → compression)
//   - Net effect: elongation along field + compression perpendicular
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5 for neighbor access

  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T Hmag  = mf_cell.template get<HMAG<T>>();
  T phi   = pf_cell.template get<typename PFCELL::GenericRho>();
  T chi   = chi_l + phi * (chi_h - chi_l);

  // ∇|H| via D2Q5 4-direction isotropic gradient
  Vector<T, 2> grad_Hmag{0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_Hmag[0] += wk * ck[0] * Hmag_k;
    grad_Hmag[1] += wk * ck[1] * Hmag_k;
  }
  grad_Hmag[0] /= LatSet::cs2;
  grad_Hmag[1] /= LatSet::cs2;

  // ∇χ via D2Q5 4-direction isotropic gradient (from per-cell CHI field)
  Vector<T, 2> grad_chi{0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T chi_k = mf_cell.getNeighbor(k).template get<CHI_PERCELL<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_chi[0] += wk * ck[0] * chi_k;
    grad_chi[1] += wk * ck[1] * chi_k;
  }
  grad_chi[0] /= LatSet::cs2;
  grad_chi[1] /= LatSet::cs2;

  // F_mag = χ·|H|·∇|H| - ½·|H|²·∇χ  (complete KH force divergence)
  T Hsq = Hmag * Hmag;
  T Fmag_x = chi * Hmag * grad_Hmag[0] - T{0.5} * Hsq * grad_chi[0];
  T Fmag_y = chi * Hmag * grad_Hmag[1] - T{0.5} * Hsq * grad_chi[1];

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

}  // namespace mfield
