// moment.h — macroscopic moment computation (used by boundary conditions)
#pragma once
#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"

namespace moment {
template <typename T, typename LatSet>
struct Rho {
  static void apply(BasicPopCell<T,LatSet>& cell, T& rho) {
    rho = T(0);
    for (unsigned int k=0; k<LatSet::q; ++k) rho += cell[k];
  }
};
template <typename T, typename LatSet>
struct Velocity {
  static void apply(BasicPopCell<T,LatSet>& cell, Vector<T,LatSet::d>& vel) {
    T rho = T(0); Vector<T,LatSet::d> t{};
    for (unsigned int k=0; k<LatSet::q; ++k) { rho+=cell[k]; t=t+cell[k]*latset::c<LatSet>(k); }
    vel = (rho>T(0)) ? t*(T(1)/rho) : Vector<T,LatSet::d>{};
  }
};
}
