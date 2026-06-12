#pragma once

#include "data_struct/cell.h"
#include "utils/alias.h"

namespace moment {

// ===================================================================
//  Section 1: Phase-Field Macrosopic (f_α → φ)
//
//  φ = Σ_α f_α,  write back to GenericRho (PHI) field
// ===================================================================

template <typename CELL>
struct PhaseFieldPhi {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T phi = T{0};
    for (unsigned int i = 0; i < LatSet::q; ++i) phi += cell[i];
    cell.template get<typename CELL::GenericRho>() = phi;
  }
};

// ===================================================================
//  Section 2: Velocity-Based Macroscopic (g_α → p, u*)
//
//  p   = ρ·c_s²·Σ_α g_α
//  u*  = Σ_α g_α·e_α              
// ===================================================================

template <typename CELL>
struct VelocityMomenta {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T rho = cell.template get<typename CELL::GenericRho>();
    T sum = T{0};
    Vector<T, LatSet::d> u_star = T{};
    auto& F = cell.template get<FORCE<T, LatSet::d>>();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      sum += cell[i];
      u_star += latset::c<LatSet>(i) * cell[i];
    }
    cell.template get<VELOCITY<T,LatSet::d>>() = u_star + F * T{0.5} / rho;
    cell.template get<PRESSURE<T>>() = rho * LatSet::cs2 * sum;
  }
};

// ===================================================================
//  Section 3: Physical Property Interpolation (φ → ρ, η)
//
//  ρ = ρ_l + φ·(ρ_h - ρ_l)      Eq.(49)
//  η = η_l + φ·(η_h - η_l)      Eq.(50)
//
//  PF cell provides φ (GenericRho) and consts (RHO_L, RHO_H, ETA_L, ETA_H)
//  NS cell receives ρ (GenericRho) and η (PHYSICAL_ETA)
// ===================================================================

template <typename PFCELL, typename NSCELL>
struct PropertyInterpolation {
  using T = typename PFCELL::FloatType;

  __any__ static inline void apply(const PFCELL& pf_cell, NSCELL& ns_cell) {
    const T phi   = pf_cell.template get<typename PFCELL::GenericRho>();
    const T rho_l = pf_cell.template get<RHO_L<T>>();
    const T rho_h = pf_cell.template get<RHO_H<T>>();
    const T eta_l = pf_cell.template get<ETA_L<T>>();
    const T eta_h = pf_cell.template get<ETA_H<T>>();

    ns_cell.template get<typename NSCELL::GenericRho>() = rho_l + phi * (rho_h - rho_l);
    ns_cell.template get<PHYSICAL_ETA<T>>() = eta_l + phi * (eta_h - eta_l);
  }
};

// ===================================================================
//  Section 4: Legacy moment types (backward compat for boundary code)
// ===================================================================

template <typename T, typename LatSet>
struct Rho {
  static void apply(const BasicPopCell<T, LatSet>& cell, T& rho) {
    rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) rho += cell[i];
  }
};

template <typename T, typename LatSet>
struct Velocity {
  static void apply(const BasicPopCell<T, LatSet>& cell, Vector<T, LatSet::d>& u) {
    u.clear();
    T rho = T(0);
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u += latset::c<LatSet>(i) * cell[i];
    }
    u = u / rho;
  }
};

template <typename T, typename LatSet>
struct RhoVelocity {
  static void apply(const BasicPopCell<T, LatSet>& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T(0);
    u.clear();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rho += cell[i];
      u += latset::c<LatSet>(i) * cell[i];
    }
    u = u / rho;
  }
};

}  // namespace moment
