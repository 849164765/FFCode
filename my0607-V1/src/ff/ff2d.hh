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

// ---- FFRhoUCompute2D ----
// Post-stream macroscopic computation: ρ = Σ f_i, u* = Σ c_i f_i / ρ
// Writes results to RHO and VELOCITY fields. No force correction.
template <typename CELL>
__any__ void FFRhoUCompute2D<CELL>::apply(CELL& cell) {
  T rho_value = T{};
  Vector<T, LatSet::d> u_value{};

  for (unsigned int i = 0; i < LatSet::q; ++i) {
    rho_value += cell[i];
    u_value += latset::c<LatSet>(i) * cell[i];
  }
  u_value /= rho_value;

  cell.template get<RHO<T>>() = rho_value;
  cell.template get<VELOCITY<T, LatSet::d>>() = u_value;
}

// ---- FFFinalVelocityCompute2D ----
// Force-corrected velocity: u = u* + F/(2·ρ(φ))
// Reads RHO (ρ(φ)) and FORCE; adds force correction to VELOCITY.
template <typename CELL>
__any__ void FFFinalVelocityCompute2D<CELL>::apply(CELL& cell) {
  const T rho = cell.template get<RHO<T>>();
  const auto& force = cell.template get<FORCE<T, LatSet::d>>();

  auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
  u += force * (T{0.5} / rho);
}

// ---- FFDensityForce2D ----
// Computes F_p and F_v per Eqs.(27)-(28) and adds them to ns_cell FORCE.
//
// F_p = -cs²·∇ρ                    ∇ρ = (ρ_h - ρ_l)·∇φ
// F_v = ν·(∇u + u∇)·∇ρ            ν = η(φ) / ρ(φ)
//
// Requires pre-computed PF fields: PHI, GRAD, RHO_L, RHO_H, ETA_L, ETA_H
// Requires pre-computed NS fields: RHO, VELOCITY (u*)
template <typename PFCELL, typename NSCELL>
__any__ void FFDensityForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  // ---- PF side: φ, ∇φ, ρ(φ), ν(φ) ----
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  T rho_phi = rho_l + phi * (rho_h - rho_l);
  T eta_phi = eta_l + phi * (eta_h - eta_l);
  T nu_phi = eta_phi / rho_phi;

  const auto& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T delta_rho = rho_h - rho_l;
  Vector<T, LatSet::d> grad_rho;
  grad_rho[0] = delta_rho * grad_phi[0];
  grad_rho[1] = delta_rho * grad_phi[1];

  // ---- NS side: velocity gradient ∇u via isotropic FD stencil ----
  T dux_dx = T{}, dux_dy = T{};
  T duy_dx = T{}, duy_dy = T{};

  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const auto& u_neighbor =
        ns_cell.getNeighbor(k).template get<VELOCITY<T, LatSet::d>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);

    dux_dx += wk * ck[0] * u_neighbor[0];
    dux_dy += wk * ck[1] * u_neighbor[0];
    duy_dx += wk * ck[0] * u_neighbor[1];
    duy_dy += wk * ck[1] * u_neighbor[1];
  }
  T inv_cs2 = T{1} / LatSet::cs2;
  dux_dx *= inv_cs2;
  dux_dy *= inv_cs2;
  duy_dx *= inv_cs2;
  duy_dy *= inv_cs2;

  // ---- F_p = -cs² · ∇ρ   Eq.(28), p/ρ = cs² in standard LBM EOS ----
  Vector<T, LatSet::d> fp;
  fp[0] = -LatSet::cs2 * grad_rho[0];
  fp[1] = -LatSet::cs2 * grad_rho[1];

  // ---- F_v = ν · (∇u + u∇) · ∇ρ ----
  // (∇u + u∇)_{ab} = ∂_a u_b + ∂_b u_a
  // F_v^a = ν · Σ_b (∂_a u_b + ∂_b u_a) · ∂_b ρ
  T vxx = T{2} * dux_dx;        // ∂_x u_x + ∂_x u_x
  T vxy = dux_dy + duy_dx;      // ∂_x u_y + ∂_y u_x
  T vyx = duy_dx + dux_dy;      // ∂_y u_x + ∂_x u_y (= vxy symmetric)
  T vyy = T{2} * duy_dy;        // ∂_y u_y + ∂_y u_y

  Vector<T, LatSet::d> fv;
  fv[0] = nu_phi * (vxx * grad_rho[0] + vxy * grad_rho[1]);
  fv[1] = nu_phi * (vyx * grad_rho[0] + vyy * grad_rho[1]);

  // ---- Accumulate to NS FORCE ----
  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += fp[0] + fv[0];
  ns_force[1] += fp[1] + fv[1];
}

// ---- FFRhoOmegaInterp2D ----
// Interpolate ρ(φ)=ρ_l+φ(ρ_h-ρ_l) → NS RHO, ω(φ)=1/(0.5+ν(φ)/cs²) → NS OMEGA
// Called after PF φ is updated, before NS collision.
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaInterp2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  T eta = eta_l + phi * (eta_h - eta_l);
  T nu = eta / rho;
  T tau = T{0.5} + nu / LatSet::cs2;
  T omega = T{1} / tau;
  if (omega > T{1.95}) omega = T{1.95};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<RHO<T>>() = rho;
  ns_cell.template get<OMEGA<T>>() = omega;
}

// ---- FFVelocityCompute2D ----
// Compute u* = Σ c_i·g_i / Σ g_i from post-stream populations.
// Reads populations; writes VELOCITY only (does NOT touch RHO).
template <typename CELL>
__any__ void FFVelocityCompute2D<CELL>::apply(CELL& cell) {
  T rho_val = T{};
  Vector<T, LatSet::d> u_val{};

  for (unsigned int i = 0; i < LatSet::q; ++i) {
    rho_val += cell[i];
    u_val += latset::c<LatSet>(i) * cell[i];
  }
  u_val /= rho_val;

  cell.template get<VELOCITY<T, LatSet::d>>() = u_val;
}

// ---- FFVelocityPressureCompute2D ----
// Velocity + Pressure per Guo et al. Eqs.(35)-(36):
//   u = Σ_α g_α e_α + (dt/(2ρ)) F_total   (Eq.35, dt=1)
//   p = ρ · cs² · Σ_α g_α                  (Eq.36)
// Reads NS distributions (cell[i]), RHO (= ρ(φ)), FORCE.
// Writes VELOCITY and PRESSURE.
template <typename CELL>
__any__ void FFVelocityPressureCompute2D<CELL>::apply(CELL& cell) {
  T rho_sum = T{};
  Vector<T, LatSet::d> momentum{};

  for (unsigned int i = 0; i < LatSet::q; ++i) {
    rho_sum += cell[i];
    momentum += latset::c<LatSet>(i) * cell[i];
  }

  const T rho_phi = cell.template get<RHO<T>>();
  const auto& force = cell.template get<FORCE<T, LatSet::d>>();

  // Eq.(35): u = (Σ g_i·e_i) + F/(2ρ)
  // In paper: Σg ≈ 1, so Σg·e = velocity directly
  // In our code: Σg = ρ, Σg·e = ρu* (momentum), must divide by Σg to get velocity
  Vector<T, LatSet::d> vel;
  if (rho_sum > T{1e-12}) {
    vel = momentum  + force * (T{0.5} / rho_phi);
  } else {
    vel = Vector<T, LatSet::d>{T{0}, T{0}};
  }
  cell.template get<VELOCITY<T, LatSet::d>>() = vel;

  // Eq.(36): p = ρ(φ) · cs² · Σ g_i
  T pres = rho_phi * LatSet::cs2 * rho_sum;
  cell.template get<PRESSURE<T>>() = pres;
}

}  // namespace ff
