// equilibrium.ur.h
// Equilibrium distribution functions (included by cell.h, bounce_back_boundary.h)
// Provides new Cell-based API (equilibrium::FirstOrder/SecondOrder) from collision.h
// and backward-compatible Equilibrium<T,LatSet> for legacy boundary code.

#pragma once

#include "lbm/collision.h"

// Backward-compatible interface for bounce_back_boundary.h
template <typename T, typename LatSet>
struct Equilibrium {
  static constexpr unsigned int D = LatSet::d;

  static T Order1(unsigned int k, const Vector<T, D>& u, T rho) {
    T cu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k) * rho * (T(1) + LatSet::InvCs2 * cu);
  }
  static T Order2(unsigned int k, const Vector<T, D>& u, T rho, T /*u2*/) {
    T cu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k) * rho
           * (T(1) + T(3) * cu + T(4.5) * cu * cu - T(1.5) * (u * u));
  }
};
