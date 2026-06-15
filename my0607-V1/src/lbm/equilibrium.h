// equilibrium.h
// Equilibrium distribution function operators (struct + apply pattern)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   Eq.18: PhaseFieldEquilibrium   (f_eq for phase-field)
//   Eq.31: NSEquilibrium           (g_eq for velocity-based NS)
//   Eq.40: MagEquilibrium          (h_eq for D2Q5 magnetic solver)

#pragma once

#include <array>

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// PhaseFieldEquilibrium (Eq.18)
// f_alpha^eq = omega_alpha * phi * (1 + e_alpha·u / cs^2)
// ================================================================
template <typename CELLTYPE>
struct PhaseFieldEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& feq) {
    T phi = cell.template get<PHI<T>>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    for (int k = 0; k < LatSet::q; ++k) {
      T eu = latset::c<LatSet>(k) * u;
      feq[k] = latset::w<LatSet>(k) * phi * (T{1} + LatSet::InvCs2 * eu);
    }
  }
};

// ================================================================
// NSEquilibrium (Eq.31) — velocity-based Navier-Stokes
// g_alpha^eq = omega_alpha * [ p/(rho*cs^2) + e_alpha·u/cs^2
//             + (e_alpha·u)^2/(2*cs^4) - u^2/(2*cs^2) ]
// KEY: zero-moment is p/(rho*cs^2), NOT rho
// ================================================================
template <typename CELLTYPE>
struct NSEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& geq) {
    T p = cell.template get<PRESSURE<T>>();
    T rho = cell.template get<RHO<T>>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    T p_over_rho_cs2 = p / (rho * LatSet::cs2);
    T uSqr = u * u;
    for (int k = 0; k < LatSet::q; ++k) {
      T eu = latset::c<LatSet>(k) * u;
      geq[k] = latset::w<LatSet>(k)
             * (p_over_rho_cs2 + LatSet::InvCs2 * eu
             + LatSet::InvCs4 * eu * eu * T{0.5} - LatSet::InvCs2 * uSqr * T{0.5});
    }
  }
};

// ================================================================
// MagEquilibrium (Eq.40) — D2Q5 magnetic solver
// h_alpha^eq = omega_alpha * psi   (omega_alpha = 1/5 for all alpha)
// ================================================================
template <typename CELLTYPE>
struct MagEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& heq) {
    T psi = cell.template get<PSI<T>>();
    for (int k = 0; k < LatSet::q; ++k) {
      heq[k] = latset::w<LatSet>(k) * psi;
    }
  }
};
