#pragma once

#include "mfield/mfield3d.h"

namespace mfield {

// ---- MFUpdateCoeffs3D ----
// Update per-cell mu, chi, omega_psi from phi
// mu  = mu_l  + phi * (mu_h  - mu_l)
// chi = chi_l + phi * (chi_h - chi_l)
// tau_psi = 0.5 + mu / cs²,  omega_psi = 1 / tau_psi (clamped)
template <typename PFCELL, typename MFCELL>
__any__ void MFUpdateCoeffs3D<PFCELL, MFCELL>::apply(PFCELL& pf_cell,
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

  // tau = 0.5 + PsiSolver_K * mu  (diffusion D = Q2·(τ-0.5) = Q2·K·μ ∝ μ,
  // Q2 = Σ w c² = 1/4 for D3Q7).  Any K leaves the continuum fixed point
  // ∇·(μ∇ψ)=0 unchanged; the sub-iterated solve converges to a smooth
  // discrete fixed point that suppresses the grid-scale |H| oscillation.
  T PsiSolver_K = mf_cell.template get<PSI_K<T>>();
  T tau_psi = T{0.5} + PsiSolver_K * mu;
  T omega_psi = T{1} / tau_psi;
  if (omega_psi > T{1.95}) omega_psi = T{1.95};
  if (omega_psi < T{0.01}) omega_psi = T{0.01};

  mf_cell.template get<MU_PERCELL<T>>()   = mu;
  mf_cell.template get<CHI_PERCELL<T>>()  = chi;
  mf_cell.template get<OMEGA_PSI<T>>()    = omega_psi;
}

// ---- MFComputeH3D ----
// H = -∇ψ using D3Q7 6-direction isotropic gradient (central difference)
// ∇ψ = Σ_k w_k * c_k * ψ_k / Q2,  Q2 = Σ_j w_j * c_jα²
// D3Q7: w=1/8, Q2=1/4  →  ∇ψ_α = (ψ_α⁺ - ψ_α⁻)/2
// (D3Q7 的 Σ w c_α² = 1/4 ≠ cs² = 1/3，不能直接除以 cs²)
template <typename CELL>
__any__ void MFComputeH3D<CELL>::apply(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  T Q2 = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    Q2 += wk * ck[0] * ck[0];
  }

  T grad_x = T{0};
  T grad_y = T{0};
  T grad_z = T{0};

  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T psi_k = cell.getNeighbor(k).template get<PSI<T>>();
    T wk    = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_x += wk * ck[0] * psi_k;
    grad_y += wk * ck[1] * psi_k;
    grad_z += wk * ck[2] * psi_k;
  }
  grad_x /= Q2;
  grad_y /= Q2;
  grad_z /= Q2;

  T Hx = -grad_x;
  T Hy = -grad_y;
  T Hz = -grad_z;
  T Hmag = std::sqrt(Hx * Hx + Hy * Hy + Hz * Hz);

  cell.template get<HX<T>>()   = Hx;
  cell.template get<HY<T>>()   = Hy;
  cell.template get<HZ<T>>()   = Hz;
  cell.template get<HMAG<T>>() = Hmag;
}

// ---- MFMagneticForce3D ----
// Magnetic body force (Guo et al., Phys. Fluids 37, 022148 (2025), Eq. (8)):
//   F_m = (μ₀χ/2) ∇(|H|²) = μ₀χ|H|∇|H|   (with μ₀ ≡ 1)
// ∇|H| uses the same D3Q7 Q2-normalized isotropic gradient as MFComputeH3D.
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce3D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D3Q7 for neighbor access

  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T Hmag  = mf_cell.template get<HMAG<T>>();
  T phi   = pf_cell.template get<typename PFCELL::GenericRho>();
  T chi   = chi_l + phi * (chi_h - chi_l);

  T Q2 = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    Q2 += wk * ck[0] * ck[0];
  }

  // ∇|H| via D3Q7 6-direction isotropic gradient
  Vector<T, 3> grad_Hmag{0, 0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_Hmag[0] += wk * ck[0] * Hmag_k;
    grad_Hmag[1] += wk * ck[1] * Hmag_k;
    grad_Hmag[2] += wk * ck[2] * Hmag_k;
  }
  grad_Hmag[0] /= Q2;
  grad_Hmag[1] /= Q2;
  grad_Hmag[2] /= Q2;

  // F_mag = χ·|H|·∇|H|  (Guo 2025, Eq. (8), μ₀ ≡ 1)
  T Fmag_x = chi * Hmag * grad_Hmag[0];
  T Fmag_y = chi * Hmag * grad_Hmag[1];
  T Fmag_z = chi * Hmag * grad_Hmag[2];

  auto& ns_force = ns_cell.template get<FORCE<T, 3>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
  ns_force[2] += Fmag_z;
}

}  // namespace mfield
