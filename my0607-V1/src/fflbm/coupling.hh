// coupling.hh — PF→NS force coupling
//
// SurfaceTension<PFCELL, NSCELL> : F_s = λ·∇φ
// GravityBuoyancy<PFCELL, NSCELL> : F_b[1] = (ρ(φ)-ρ_h)·g

#pragma once

#include "ff/ff2d.h"

namespace fflbm {

template <typename PFCELL, typename NSCELL>
struct SurfaceTension {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;

  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    const T chem_potential = pf_cell.template get<ff::CHEMICALPOTENTIAL<T>>();
    const auto& grad = pf_cell.template get<GRAD<T, LatSet::d>>();
    auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
    for (unsigned int d = 0; d < LatSet::d; ++d) {
      ns_force[d] += chem_potential * grad[d];
    }
  }
};

template <typename PFCELL, typename NSCELL>
struct GravityBuoyancy {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;

  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    const T phi = pf_cell.template get<typename PFCELL::GenericRho>();
    const T rho_l = pf_cell.template get<ff::RHO_L<T>>();
    const T rho_h = pf_cell.template get<ff::RHO_H<T>>();
    const T g = pf_cell.template get<ff::GRAVITY<T>>();
    const T rho = rho_l + phi * (rho_h - rho_l);
    auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
    ns_force[1] += (rho - rho_h) * g;
  }
};

// Spatially varying viscosity update:
//   μ(φ) = eta_l + φ·(eta_h - eta_l)         dynamic viscosity
//   ρ(φ) = rho_l + φ·(rho_h - rho_l)         density
//   ν(φ) = μ(φ) / ρ(φ)                       kinematic viscosity
//   τ(φ) = 0.5 + ν(φ)/cs²                    relaxation time
//   ω(φ) = 1/τ(φ)  (clamped [0.1, 1.95])
// Writes omega to NS cell's OMEGA field (per-cell, used by getOmegaf())
template <typename PFCELL, typename NSCELL>
struct OmegaUpdate {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;

  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    const T phi = pf_cell.template get<typename PFCELL::GenericRho>();
    const T rho_l = pf_cell.template get<ff::RHO_L<T>>();
    const T rho_h = pf_cell.template get<ff::RHO_H<T>>();
    const T eta_l = pf_cell.template get<ff::ETA_L<T>>();
    const T eta_h = pf_cell.template get<ff::ETA_H<T>>();

    const T mu = eta_l + phi * (eta_h - eta_l);
    const T rho = rho_l + phi * (rho_h - rho_l);
    const T nu = (rho > T{0}) ? (mu / rho) : T{0};

    const T tau = T{0.5} + nu / LatSet::cs2;
    T omega = T{1} / tau;

    // clamp for numerical stability
    if (omega > T{1.95}) omega = T{1.95};
    if (omega < T{0.1})  omega = T{0.1};

    ns_cell.template get<OMEGA<T>>() = omega;
  }
};

}  // namespace fflbm
