#pragma once

#include "data_struct/Vector.h"
#include "lbm/lattice_set.h"

// legacy backward-compatible equilibrium (needed by boundary conditions)
template <typename T, typename LatSet>
struct Equilibrium {
  __any__ static inline T Order1(int k, const Vector<T, LatSet::d>& u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order2(int k, const Vector<T, LatSet::d>& u, T rho, T u2) {
    T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T(1) + LatSet::InvCs2 * uc + uc * uc * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }
};

// backward-compatible equilibrium in namespace equilibrium
namespace equilibrium {

template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d>& u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }
};

template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d>& u, T rho, T u2) {
    T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T{0.5});
  }
};

}  // namespace equilibrium
