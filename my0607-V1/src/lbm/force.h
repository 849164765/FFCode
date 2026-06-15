// force.h
// Guo force scheme — converts total force vector to discrete lattice force
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025), Eq.32
//   G_beta = omega_beta * [ e_beta·F/cs^2 + (e_beta·u)(e_beta·F)/cs^4 - u·F/cs^2 ]
//
// Usage: coupling operators (T14) fill FORCE field with F_total,
//        GuoForce reads it and adds G_beta to populations.

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// GuoForce (Eq.32)
// Reads:  FORCE (F_total, pre-filled by coupling), VELOCITY
// Writes: POP (+= G_beta to each population)
// ================================================================
template <typename CELLTYPE>
struct GuoForce {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // BGK path: compute G_beta and add directly to populations
  __any__ static inline void apply(CELL& cell) {
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    for (int k = 0; k < LatSet::q; ++k) {
      cell[k] += _compute(k, F, u);
    }
  }

  // MRT path: compute G_beta and write to output array
  // Caller then transforms G to moment space: m += (I - S/2) * M9 * G
  __any__ static inline void GetForce(CELL& cell, std::array<T, LatSet::q>& G) {
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    for (int k = 0; k < LatSet::q; ++k) {
      G[k] = _compute(k, F, u);
    }
  }

 private:
  __any__ static inline T _compute(int k, const Vector<T, LatSet::d>& F,
                                   const Vector<T, LatSet::d>& u) {
    const T uF = u * F;
    const T ekF = latset::c<LatSet>(k) * F;
    const T ekU = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k)
         * (LatSet::InvCs2 * ekF + LatSet::InvCs4 * ekU * ekF - LatSet::InvCs2 * uF);
  }
};
