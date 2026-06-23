// moment.h
// Macroscopic moment recovery operators (struct + apply pattern)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   Eq.24: PhaseFieldMoment  (phi = sum(f))        -> writes PHI
//   Eq.35: VelocityMoment    (u = sum(g*e))         -> writes VELOCITY
//   Eq.36: PressureMoment    (p = rho*cs^2*sum(g))  -> writes PRESSURE
//   Eq.44: MagPotentialMoment(psi = sum(h))         -> writes PSI

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

namespace moment {

// ================================================================
// PhaseFieldMoment (Eq.24): phi = sum_alpha f_alpha
// ================================================================
template <typename CELLTYPE>
struct PhaseFieldMoment {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T phi = T{0};
    for (int k = 0; k < LatSet::q; ++k) {
      phi += cell[k];
    }
    cell.template get<PHI<T>>() = phi;
  }
};

// ================================================================
// VelocityMoment (Eq.35): u = sum_alpha(g_alpha * e_alpha)
// Note: dt/(2*rho)*F correction applied by collision operator
// ================================================================
template <typename CELLTYPE>
struct VelocityMoment {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    Vector<T, LatSet::d> u{};
    for (int k = 0; k < LatSet::q; ++k) {
      u += cell[k] * latset::c<LatSet>(k);
    }
    cell.template get<VELOCITY<T, LatSet::d>>() = u;
  }
};

// ================================================================
// PressureMoment (Eq.36): p = rho * cs^2 * sum_alpha(g_alpha)
// ================================================================
template <typename CELLTYPE>
struct PressureMoment {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T sum = T{0};
    for (int k = 0; k < LatSet::q; ++k) {
      sum += cell[k];
    }
    T rho = cell.template get<RHO<T>>();
    cell.template get<PRESSURE<T>>() = rho * LatSet::cs2 * sum;
  }
};

// ================================================================
// MagPotentialMoment (Eq.44): psi = sum_alpha h_alpha  (D2Q5)
// ================================================================
template <typename CELLTYPE>
struct MagPotentialMoment {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T psi = T{0};
    for (int k = 0; k < LatSet::q; ++k) {
      psi += cell[k];
    }
    cell.template get<PSI<T>>() = psi;
  }
};

// ================================================================
// NSMomentum (Eqs.35-36): combined VELOCITY + PRESSURE from populations
// Runs as a SEPARATE task BEFORE NSMRT collision, so NSMRT's equilibrium
// and GuoForce read fresh VELOCITY/PRESSURE from pre-collision g[k].
// ================================================================
template <typename CELLTYPE>
struct NSMomentum {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    // Eq.35: u = sum(g * e) + dt/(2*rho)*F_total
    T rho = cell.template get<RHO<T>>();
    const auto& F_total = cell.template get<FORCE<T, LatSet::d>>();
    Vector<T, LatSet::d> u{};
    for (int k = 0; k < LatSet::q; ++k) {
      u += cell[k] * latset::c<LatSet>(k);
    }
    u += T{0.5} / rho * F_total;
    cell.template get<VELOCITY<T, LatSet::d>>() = u;

    // Eq.36: p = rho * cs^2 * sum(g)
    T sum_g = T{0};
    for (int k = 0; k < LatSet::q; ++k) sum_g += cell[k];
    cell.template get<PRESSURE<T>>() = rho * LatSet::cs2 * sum_g;
  }
};

}  // namespace moment
