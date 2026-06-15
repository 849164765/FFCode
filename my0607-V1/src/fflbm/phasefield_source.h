// phasefield_source.h
// Phase-field source term operator (Eq.19)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   F_beta = omega_beta * (e_beta · n) * 4 * phi * (1 - phi) / W
//
// This is the Allen-Cahn interface sharpening source term.
// Writes to populations (BGK) or output array (MRT).

#pragma once

#include <array>

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// PhaseFieldSource (Eq.19)
// Reads:  PHI, PHASENORMAL, INTERFACEWIDTH
// Writes: POP (+= F_beta) [BGK]  or  output array [MRT]
// ================================================================
template <typename CELLTYPE>
struct PhaseFieldSource {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // BGK path: compute F_beta and add directly to populations
  __any__ static inline void apply(CELL& cell) {
    T phi = cell.template get<PHI<T>>();
    const auto& n = cell.template get<PHASENORMAL<T, LatSet::d>>();
    T invW = T{1} / cell.template get<INTERFACEWIDTH<T>>();
    T factor = T{4} * phi * (T{1} - phi) * invW;
    for (int k = 0; k < LatSet::q; ++k) {
      cell[k] += _compute(k, n, factor);
    }
  }

  // MRT path: compute F_beta and write to output array
  // Caller transforms to moment space: m += (I - S/2) * M9 * F
  __any__ static inline void GetSource(CELL& cell,
                                       std::array<T, LatSet::q>& F) {
    T phi = cell.template get<PHI<T>>();
    const auto& n = cell.template get<PHASENORMAL<T, LatSet::d>>();
    T invW = T{1} / cell.template get<INTERFACEWIDTH<T>>();
    T factor = T{4} * phi * (T{1} - phi) * invW;
    for (int k = 0; k < LatSet::q; ++k) {
      F[k] = _compute(k, n, factor);
    }
  }

 private:
  __any__ static inline T _compute(int k,
                                   const Vector<T, LatSet::d>& n,
                                   T factor) {
    return latset::w<LatSet>(k) * (latset::c<LatSet>(k) * n) * factor;
  }
};
