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
// F_s = λ * ∇φ  → added to ns_cell FORCE
template <typename PFCELL, typename NSCELL>
__any__ void FFSurfaceTension2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T chem_potential = pf_cell.template get<CHEMICALPOTENTIAL<T>>();
  const Vector<T, LatSet::d>& grad = pf_cell.template get<GRAD<T, LatSet::d>>();

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += chem_potential * grad[0];
  ns_force[1] += chem_potential * grad[1];
}

// ---- FFGravityForce2D ----
// F_b[1] = -|gravity| * ρ(φ)  where ρ(φ) = ρ_l + φ*(ρ_h-ρ_l)
template <typename PFCELL, typename NSCELL>
__any__ void FFGravityForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T g = pf_cell.template get<GRAVITY<T>>();

  // g is negative (downward), so |g| = -g
  T rho = rho_l + phi * (rho_h - rho_l);
  // F_b in y-direction. Gravity points downward (negative y).
  // Since g is negative: buoyancy = rho*g (negative * negative = positive upward force if rho_l < rho_h)
  // Wait: F_b = ρ * G_y where G_y is the gravitational acceleration.
  // For bubble rising, the bubble is lighter (ρ_l << ρ_h), so light fluid rises.
  // F_b[1] = rho * g where g is negative (pointing down).
  // The NS solver handles the sign convention.
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += (rho - rho_h) * g;
}

// ---- FFRhoOmegaUpdate2D ----
// WARNING: Cell::getOmega() returns the block-level scalar (from converter),
// NOT the per-cell OMEGA<T> field. So per-cell omega updates have NO effect on
// the collision. For variable viscosity, modify the collision operator to call
// getOmegaf() instead of getOmega(). See src/data_struct/cell.h:168.
//
// ρ = ρ_l + φ*(ρ_h - ρ_l)
// ν = η/ρ, τ = 0.5 + ν/cs², omega = 1/τ
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaUpdate2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  // Interpolate rho
  T rho = rho_l + phi * (rho_h - rho_l);
  ns_cell.template get<typename NSCELL::GenericRho>() = rho;

  // Interpolate eta (dynamic viscosity)
  T eta = eta_l + phi * (eta_h - eta_l);

  // Compute kinematic viscosity and omega
  T nu = eta / rho;
  T tau = T{0.5} + nu / LatSet::cs2;
  // Clamp omega to avoid instability
  T omega = T{1} / tau;
  if (omega > T{1.95}) omega = T{1.95};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<OMEGA<T>>() = omega;
}

// ---- FFPressForce2D ----
// F_p = -(p/ρ)·∇ρ  where ∇ρ = (ρ_h-ρ_l)·∇φ
template <typename PFCELL, typename NSCELL>
__any__ void FFPressForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  const auto& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T p = ns_cell.template get<PRESSURE<T>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  T inv_rho = T{1} / rho;
  T coeff = -p * inv_rho * (rho_h - rho_l);

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += coeff * grad_phi[0];
  ns_force[1] += coeff * grad_phi[1];
}

// ---- FFVisForce2D ----
// F_v = ν(∇u + u∇)·∇ρ    ∇ρ = (ρ_h-ρ_l)·∇φ
// ∇u computed via D2Q9 isotropic FD on ns_cell neighbor VELOCITY
template <typename PFCELL, typename NSCELL>
__any__ void FFVisForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();
  const auto& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  T eta = eta_l + phi * (eta_h - eta_l);
  T nu = eta / rho;

  // ∇u via isotropic FD using neighbor VELOCITY (D2Q9 stencil)
  T grad_u[2][2] = {};  // grad_u[a][b] = ∂u_a/∂x_b
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const auto& u_neighbor = ns_cell.getNeighbor(k).template get<VELOCITY<T, LatSet::d>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    for (unsigned int a = 0; a < 2; ++a) {
      for (unsigned int b = 0; b < 2; ++b) {
        grad_u[a][b] += wk * ck[b] * u_neighbor[a];
      }
    }
  }
  for (unsigned int a = 0; a < 2; ++a)
    for (unsigned int b = 0; b < 2; ++b)
      grad_u[a][b] /= LatSet::cs2;

  // ∇ρ = (ρ_h-ρ_l)·∇φ
  T grad_rho[2];
  grad_rho[0] = (rho_h - rho_l) * grad_phi[0];
  grad_rho[1] = (rho_h - rho_l) * grad_phi[1];

  // F_v^a = nu * Σ_b (∂u_a/∂x_b + ∂u_b/∂x_a) * ∂ρ/∂x_b
  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  for (unsigned int a = 0; a < 2; ++a) {
    T fv = T{};
    for (unsigned int b = 0; b < 2; ++b) {
      fv += nu * (grad_u[a][b] + grad_u[b][a]) * grad_rho[b];
    }
    ns_force[a] += fv;
  }
}

}  // namespace ff
