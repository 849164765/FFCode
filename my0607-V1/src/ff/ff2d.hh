#pragma once

#include "ff/ff2d.h"

namespace ff {

template <typename CELL>
__any__ void FF2D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad;
  grad[0] = T{0};
  grad[1] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T phi_i = cell.getNeighbor(i).template get<GenericRho>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad[0] += wi * ci[0] * phi_i;
    grad[1] += wi * ci[1] * phi_i;
  }
  grad[0] /= LatSet::cs2;
  grad[1] /= LatSet::cs2;

  cell.template get<GRAD<T, LatSet::d>>() = grad;

  T grad_mag = grad.getnorm();
  // ε threshold for noise suppression: bulk regions have |∇φ| << 0.005
  // At the interface center: |∇φ| ≈ 1/W ≈ 0.33 (W=3) >> 0.005
  // This filters out noise-driven source terms while keeping interface dynamics
  T delta = T(0.005); 
  Vector<T, LatSet::d> n;
  if (grad_mag < delta) {
    n[0] = T{0};
    n[1] = T{0};
  } else {
    n[0] = grad[0] / grad_mag;
    n[1] = grad[1] / grad_mag;
  }
  cell.template get<NORMAL<T, LatSet::d>>() = n;
}

// ---- FFLaplacian2D ----
// ∇²φ = (2/(cs²*Δt²)) * Σ_{α=0}^{8} w_α * [φ(x+e_α) - φ(x)]
// For lattice units (Δt=1, Δx=1): factor = 2/cs²
template <typename CELL>
__any__ void FFLaplacian2D<CELL>::apply(CELL& cell) {
  T phi_self = cell.template get<GenericRho>();
  T laplacian = T{0};

  for (unsigned int k = 0; k < LatSet::q; ++k) {
    T phi_k = cell.getNeighbor(k).template get<GenericRho>();
    T wk = latset::w<LatSet>(k);
    // sum w_k * (phi_neighbor - phi_self)
    laplacian += wk * (phi_k - phi_self);
  }
  // factor 2/cs² for lattice units
  laplacian *= T{2} / LatSet::cs2;

  cell.template get<LAPLACIAN<T>>() = laplacian;
}

// ---- FFChemPotential2D ----
// λ = 4β * φ*(φ-1)*(φ-0.5) - κ * ∇²φ
template <typename CELL>
__any__ void FFChemPotential2D<CELL>::apply(CELL& cell) {
  T phi = cell.template get<GenericRho>();
  T laplacian = cell.template get<LAPLACIAN<T>>();
  T beta = cell.template get<BETA<T>>();
  T kappa = cell.template get<KAPPA<T>>();

  // φ*(φ-1)*(φ-0.5)
  T double_well = phi * (phi - T{1}) * (phi - T{0.5});
  T chem_potential = T{4} * beta * double_well - kappa * laplacian;

  cell.template get<CHEMICALPOTENTIAL<T>>() = chem_potential;
}

// ---- FFChemPotentialGradient2D ----
// ∇λ via isotropically centered FD, written to NORMAL (overwriting n=∇φ/|∇φ|)
// BGKSource then uses ∇λ (via NORMAL) as the Cahn-Hilliard source
// source_i = fomega * w_i * (c_i·∇λ) * 4φ(1-φ)/W
// The 4φ(1-φ) factor localizes the source to the interface; ∇λ provides the
// correct driving force from both phase separation and surface tension.
template <typename CELL>
__any__ void FFChemPotentialGradient2D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad_lambda;
  grad_lambda[0] = T{0};
  grad_lambda[1] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T lambda_i = cell.getNeighbor(i).template get<CHEMICALPOTENTIAL<T>>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad_lambda[0] += wi * ci[0] * lambda_i;
    grad_lambda[1] += wi * ci[1] * lambda_i;
  }
  grad_lambda[0] /= LatSet::cs2;
  grad_lambda[1] /= LatSet::cs2;

  // Store in NORMAL (replaces interface normal, used by BGKSource as source field)
  cell.template get<NORMAL<T, LatSet::d>>() = grad_lambda;
}

// ---- FFSurfaceTension2D ----
// F_s = λ * ∇φ  → added to ns_cell FORCE (uniform scaling, no per-cell density division)
template <typename PFCELL, typename NSCELL>
__any__ void FFSurfaceTension2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T chem_potential = pf_cell.template get<CHEMICALPOTENTIAL<T>>();
  const Vector<T, LatSet::d>& grad = pf_cell.template get<GRAD<T, LatSet::d>>();

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += chem_potential * grad[0];
  ns_force[1] += chem_potential * grad[1];
}

// ---- FFSurfaceTensionChemPot2D ----
// F_s = -φ * ∇λ  → added to ns_cell FORCE (uniform scaling, no per-cell density division)
template <typename PFCELL, typename NSCELL>
__any__ void FFSurfaceTensionChemPot2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();

  // Compute ∇λ from neighbor chemical potentials
  Vector<T, LatSet::d> grad_lambda;
  grad_lambda[0] = T{0};
  grad_lambda[1] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T lambda_i = pf_cell.getNeighbor(i).template get<CHEMICALPOTENTIAL<T>>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad_lambda[0] += wi * ci[0] * lambda_i;
    grad_lambda[1] += wi * ci[1] * lambda_i;
  }
  grad_lambda[0] /= LatSet::cs2;
  grad_lambda[1] /= LatSet::cs2;

  // F_s = -φ * ∇λ
  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] -= phi * grad_lambda[0];
  ns_force[1] -= phi * grad_lambda[1];
}

// ---- FFGravityForce2D ----
// F_b[1] = -|gravity| * ρ(φ)  where ρ(φ) = ρ_l + φ*(ρ_h-ρ_l)
template <typename PFCELL, typename NSCELL>
__any__ void FFGravityForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T g = pf_cell.template get<GRAVITY<T>>();

  // g is negative (downward).
  // Boussinesq recalibration: divide by ρ_l to restore physical acceleration ratio.
  // Effective buoyancy acceleration = (ρ_phys - ρ_l) / ρ_l * g
  // The NS solver handles the sign convention.
  T rho = rho_l + phi * (rho_h - rho_l);
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += (rho - rho_l) * g / rho_l;
}

// ---- FFRhoOmegaUpdate2D ----
// Density: smooth polynomial  g(φ) = φ²(3-2φ),  ρ(φ) = ρ_l + g(φ)·(ρ_h-ρ_l)
// Viscosity: harmonic average of relaxation times
//   Boussinesq fix: ν = η/ρ₀ = η (ρ₀ = 1.0 is the LBM inertial density)
//   NOT η/ρ_phys, because the LBM momentum equation uses ρ₀=1 everywhere.
//   τ_L = 0.5 + ν_L/cs², τ_G = 0.5 + ν_G/cs²
//   1/(τ(φ)-0.5) = (1-φ)/(τ_L-0.5) + φ/(τ_G-0.5)
//   ω(φ) = 1/τ(φ), clamped to [0.01, 1.95]
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaUpdate2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  // ---- Density: smooth polynomial g(φ) = φ²(3-2φ) ----
  T g_phi = phi * phi * (T{3} - T{2} * phi);
  T rho = rho_l + g_phi * (rho_h - rho_l);
  ns_cell.template get<typename NSCELL::GenericRho>() = rho;

  // ---- Viscosity: harmonic mean, ν = η (Boussinesq: ρ₀=1) ----
  T nu_L = eta_l;  // was: eta_l / rho_l
  T nu_G = eta_h;  // was: eta_h / rho_h
  T inv_cs2 = T{1} / LatSet::cs2;

  T tau_L_m05 = nu_L * inv_cs2;  // τ_L - 0.5
  T tau_G_m05 = nu_G * inv_cs2;  // τ_G - 0.5

  T inv_tau_m05 = (T{1} - phi) / tau_L_m05 + phi / tau_G_m05;
  T tau = T{0.5} + T{1} / inv_tau_m05;

  T omega = T{1} / tau;
  if (omega > T{1.95}) omega = T{1.95};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<OMEGA<T>>() = omega;
}

}  // namespace ff
