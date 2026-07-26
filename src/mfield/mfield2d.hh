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
// Magnetic body force from the divergence of the Maxwell stress tensor:
//   τ_m = μHH - (μ/2)|H|²I
//   F_m = ∇·τ_m = (H·∇μ)H - (1/2)|H|²∇μ
//
// With μ = 1+χ (μ₀ ≡ 1 in LBM units), ∇μ = ∇χ:
//   F_m = (H·∇χ)H - (1/2)|H|²∇χ
//
// This is the COMPLETE expression for a linear magnetic material with
// spatially varying permeability, derived from the Maxwell stress tensor
// without approximation. It captures both the anisotropic interfacial
// force (H·∇χ)H and the isotropic "magnetic pressure" term (1/2)|H|²∇χ.
//
// IMPORTANT: This is the ONLY formula that produces real magnetic deformation.
// The (H·∇χ)H term is NOT a gradient and cannot be absorbed by pressure.
// For a uniform field H = (0, H₀):
//   - At the top pole (n ∥ H):  F_m =  (1/2)H₀²∇χ  (outward, elongates)
//   - At the side equator (n ⟂ H): F_m = -(1/2)H₀²∇χ  (inward, compresses)
//   → Net effect: bubble elongates in the field direction.
//
// Reference: Rosensweig, "Ferrohydrodynamics" (1985)
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5 for neighbor access

  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T Hx    = mf_cell.template get<HX<T>>();
  T Hy    = mf_cell.template get<HY<T>>();
  T Hmag  = mf_cell.template get<HMAG<T>>();

  // ∇χ = (chi_h - chi_l) * ∇φ via D2Q5 4-direction isotropic gradient
  T dchi_x = T{0}, dchi_y = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T phi_k = mf_cell.getNeighbor(k).template get<PHI<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    dchi_x += wk * ck[0] * phi_k;
    dchi_y += wk * ck[1] * phi_k;
  }
  dchi_x = (chi_h - chi_l) * dchi_x / LatSet::cs2;
  dchi_y = (chi_h - chi_l) * dchi_y / LatSet::cs2;

  // F_m = (H·∇χ)H - (1/2)|H|²∇χ  (Maxwell stress tensor divergence)
  T HdotDchi = Hx * dchi_x + Hy * dchi_y;
  T half_Hsq = T{0.5} * Hmag * Hmag;
  T Fmag_x = Hx * HdotDchi - half_Hsq * dchi_x;
  T Fmag_y = Hy * HdotDchi - half_Hsq * dchi_y;

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

// ---- MFMagneticForceInterfacial2D ----
// Interfacial force only: F_m = -(1/2)|H|²∇χ
// This is the interfacial part of the Maxwell stress tensor divergence,
// representing the "magnetic pressure" term. It captures the normal
// magnetic force at the interface due to the jump in |H|².
//
// NOTE: This formula is mathematically equivalent to the Paper formula
// (χ/2)∇|H|² in incompressible flow. The difference is a pure gradient
// (1/2)∇(χ|H|²) that gets absorbed by the pressure field. Neither
// produces real magnetic deformation.
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForceInterfacial2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;

  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T Hx    = mf_cell.template get<HX<T>>();
  T Hy    = mf_cell.template get<HY<T>>();
  T Hmag  = mf_cell.template get<HMAG<T>>();

  // ∇χ = (chi_h - chi_l) * ∇φ via D2Q5 4-direction isotropic gradient
  T dchi_x = T{0}, dchi_y = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T phi_k = mf_cell.getNeighbor(k).template get<PHI<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    dchi_x += wk * ck[0] * phi_k;
    dchi_y += wk * ck[1] * phi_k;
  }
  dchi_x = (chi_h - chi_l) * dchi_x / LatSet::cs2;
  dchi_y = (chi_h - chi_l) * dchi_y / LatSet::cs2;

  // F_m = -(1/2)|H|²∇χ  (interfacial force only)
  T half_Hsq = T{0.5} * Hmag * Hmag;
  T Fmag_x = -half_Hsq * dchi_x;
  T Fmag_y = -half_Hsq * dchi_y;

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

// ---- MFMagneticForcePaper2D ----
// Paper Eq.(8): F_m = (χ/2)∇|H|²  (Kelvin body force)
// This is the formula used in Guo et al. (2025), Phys. Fluids 37, 022148.
// It is valid when χ is approximately constant (sharp interface limit).
// The force is applied as a volume force proportional to χ.
//
// NOTE: In incompressible flow, this formula is mathematically equivalent
// to the Interfacial formula -(1/2)|H|²∇χ. The difference is the pure
// gradient (1/2)∇(χ|H|²) which is absorbed by pressure. Neither produces
// real magnetic deformation.
//
// Reference:
//   Guo et al., "Phase-field lattice Boltzmann model with adaptive mesh
//   refinement for ferrofluid interfacial dynamics", Phys. Fluids 37,
//   022148 (2025), Eq.(8).
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForcePaper2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5 for neighbor access

  T chi    = mf_cell.template get<CHI_PERCELL<T>>();
  T Hx     = mf_cell.template get<HX<T>>();
  T Hy     = mf_cell.template get<HY<T>>();
  T Hmag   = mf_cell.template get<HMAG<T>>();

  // ∇|H|² via D2Q5 4-direction isotropic gradient
  T dHsq_x = T{0}, dHsq_y = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T Hsq_k  = Hmag_k * Hmag_k;
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    dHsq_x += wk * ck[0] * Hsq_k;
    dHsq_y += wk * ck[1] * Hsq_k;
  }
  dHsq_x = dHsq_x / LatSet::cs2;
  dHsq_y = dHsq_y / LatSet::cs2;

  // F_m = (χ/2) * ∇|H|²
  T half_chi = T{0.5} * chi;
  T Fmag_x = half_chi * dHsq_x;
  T Fmag_y = half_chi * dHsq_y;

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

}  // namespace mfield
