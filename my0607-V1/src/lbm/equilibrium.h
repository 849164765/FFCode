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

template <typename CELLTYPE>
struct PhaseFieldEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(std::array<T, LatSet::q>& feq, T phi,
                                   const Vector<T, LatSet::d>& u) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T edotu = u * latset::c<LatSet>(k);
      feq[k] = latset::w<LatSet>(k) * phi * (T{1} + LatSet::InvCs2 * edotu);
    }
  }
};

template <typename CELLTYPE>
struct NSFieldEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline void apply(std::array<T, LatSet::q>& feq, T rho, T p,
                                   const Vector<T, LatSet::d>& u) {
    
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T edotu = u * latset::c<LatSet>(k);
      feq[k] = latset::w<LatSet>(k) * ( p * LatSet::InvCs2 / rho + LatSet::InvCs2 * edotu + 
               edotu * edotu * T{0.5} * LatSet::InvCs4 - u.getnorm2() * T{0.5} * LatSet::InvCs2);
    }
  }
};


}  // namespace equilibrium
