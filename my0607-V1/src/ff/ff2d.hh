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


template <typename PFCELL, typename NSCELL>
__any__ void FFPressForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T press = ns_cell.template get<PRESSURE<T>>();
  T rho = ns_cell.template get<typename NSCELL::GenericRho>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  const Vector<T, LatSet::d>& grad = pf_cell.template get<GRAD<T, LatSet::d>>();
  T temp = -press * (rho_h - rho_l)  / rho; // F = -∇p/ρ, grad is ∇p

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += temp * grad[0];
  ns_force[1] += temp * grad[1];
}

template <typename PFCELL, typename NSCELL>
__any__ void FFVisForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T dux_dx = T{};
  T dux_dy = T{};
  T duy_dx = T{};
  T duy_dy = T{};
  T dphi_dx = T{};
  T dphi_dy = T{};
  for (unsigned int i = 1; i < LatSet::q; ++i) {
    const auto& ci = latset::c<LatSet>(i);
    T wi = latset::w<LatSet>(i);
    T ux_i = ns_cell.getNeighbor(i).template get<VELOCITY<T, LatSet::d>>()[0];
    T uy_i = ns_cell.getNeighbor(i).template get<VELOCITY<T, LatSet::d>>()[1];
    T phi_i = pf_cell.getNeighbor(i).template get<typename PFCELL::GenericRho>();
    dux_dx += wi * ci[0] * ux_i;
    dux_dy += wi * ci[1] * ux_i;
    duy_dx += wi * ci[0] * uy_i;
    duy_dy += wi * ci[1] * uy_i;
    dphi_dx += wi * ci[0] * phi_i;
    dphi_dy += wi * ci[1] * phi_i;
  }
  dux_dx /= LatSet::cs2;
  dux_dy /= LatSet::cs2;
  duy_dx /= LatSet::cs2;
  duy_dy /= LatSet::cs2;
  dphi_dx /= LatSet::cs2;
  dphi_dy /= LatSet::cs2;
  T omega = ns_cell.template get<OMEGA<T>>();                         
  T nu = LatSet::cs2 * (T(1) / omega - T(0.5));
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T rho_l = pf_cell.template get<RHO_L<T>>();

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  // Eq.(27): F_ν = ν (∇u + u∇) · ∇ρ
  // F_ν,x = ν(ρ_h-ρ_l) · [2·∂_xu_x·∂_xφ + (∂_yu_x+∂_xu_y)·∂_yφ]
  ns_force[0] += nu * (rho_h - rho_l) * (T{2} * dux_dx * dphi_dx + (dux_dy + duy_dx) * dphi_dy);
  ns_force[1] += nu * (rho_h - rho_l) * (T{2} * duy_dy * dphi_dy + (duy_dx + dux_dy) * dphi_dx);

}

// ---- FFRhoOmegaUpdate2D ----
// For velocity-based LBM: interpolates ρ(φ) and ν(φ) from PF φ field.
// Writes: NS RHO = ρ_l + φ*(ρ_h-ρ_l) — used in Eq.(35) u = Σge + F/(2ρ)
//         NS OMEGA = 1/(0.5 + ν/cs²) — per-cell variable viscosity in BGK collision.
// ν = η/ρ, τ = 0.5 + ν/cs², omega = 1/τ
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaUpdate2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  // Interpolate rho and write to NS field (for Eq.35 denominator)
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

// ---- FFVelocityPressureUpdate2D ----
// Velocity-based NS macro update from post-stream populations.
// Eq.(35): u = Σ e_α g_α + F/(2ρ)    (ρ = ρ(φ) from RHO field)
// Eq.(36): p = ρ c_s^2 Σ g_α
// Reads: POP (post-stream), FORCE, RHO (= ρ(φ))
// Writes: VELOCITY, PRESSURE.  Does NOT touch RHO.
template <typename CELL>
__any__ void FFVelocityPressureUpdate2D<CELL>::apply(CELL& cell) {
  T g_sum = T{};
  Vector<T, LatSet::d> momentum{};
  for (unsigned int k = 0; k < LatSet::q; ++k) {
    T gk = cell[k];
    g_sum += gk;
    const auto& ck = latset::c<LatSet>(k);
    for (unsigned int d = 0; d < LatSet::d; ++d) {
      momentum[d] += gk * static_cast<T>(ck[d]);
    }
  }
  const auto& force = cell.template get<FORCE<T, LatSet::d>>();
  T rho = cell.template get<typename CELL::GenericRho>();  // ρ(φ)

  // Eq.(35): u = Σ e_α g_α + F/(2ρ) — momentum IS already velocity,
  // only the force correction gets divided by ρ.
  auto& vel = cell.template get<VELOCITY<T, LatSet::d>>();
  for (unsigned int d = 0; d < LatSet::d; ++d) {
    vel[d] = momentum[d] + force[d] * T{0.5} / rho;
  }

  // Eq.(36): p = ρ c_s^2 Σ g_α
  cell.template get<PRESSURE<T>>() = rho * LatSet::cs2 * g_sum;
}

}  // namespace ff
