// force.h
// Guo forcing scheme for Lattice Boltzmann Method
//
// Applies a body force F to post-collision populations using the
// second-order accurate Guo scheme (Guo et al., 2002):
//
//   G_k = (1 - 0.5*omega) * w_k * [ (c_k - u)·F / cs² + (c_k·u)(c_k·F) / cs⁴ ]
//
// Usage:
//   using Gravity = force::Guo<CONSTFORCE<T, D>>;
//   using Entry = tmp::Key_TypePair<BulkFlag, Gravity>;
//   // ...
//   lattice.template ApplyCellDynamics<TaskSelector>(FlagFM);

#pragma once

#include "lbm/lattice_set.h"

namespace force {

template <typename ForceField>
struct Guo {
  template <typename CELL>
  static void apply(CELL& cell) {
    using T = typename CELL::FloatType;
    using LatSet = typename CELL::LatticeSet;
    static constexpr unsigned int D = LatSet::d;
    static constexpr unsigned int Q = LatSet::q;

    const T omega = cell.getOmega();
    const auto& u = cell.template get<VELOCITY<T, D>>();
    const auto& F = cell.template get<ForceField>();

    const T factor = T(1) - T(0.5) * omega;
    const T invCs2 = LatSet::InvCs2;
    const T invCs4 = LatSet::InvCs4;

    T u_dot_F = u * F;

    for (unsigned int k = 0; k < Q; ++k) {
      const auto& ck = latset::c<LatSet>(k);
      T ck_dot_F = ck * F;
      T ck_dot_u = ck * u;
      T term1 = (ck_dot_F - u_dot_F) * invCs2;
      T term2 = ck_dot_u * ck_dot_F * invCs4;
      T Gk = factor * latset::w<LatSet>(k) * (term1 + term2);
      cell[k] += Gk;
    }
  }
};

}  // namespace force
