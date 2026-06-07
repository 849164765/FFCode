// equlibrium.h

#pragma once
// lattice boltzmann method implementations

// #include "data_struct/Vector.h"
#include "lbm/lattice_set.h"

template <typename T, typename LatSet, typename TypePack>
class Cell;

namespace cudev {
  template <typename T, typename LatSet, typename TypePack>
  class Cell;
}  // namespace cudev

// calc equilibrium distribution function
// sum(feq_i) = rho, for both first and second order
template <typename T, typename LatSet>
struct Equilibrium {
  __any__ static inline T Order1(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order1_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order2(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T(1) + LatSet::InvCs2 * uc + uc * uc * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }

  // __any__ static inline T Order2_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho, T
  // u2) {
  //   T uc = u * latset::c<LatSet>(k);
  //   return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * uc + uc * uc * T(0.5) *
  //   LatSet::InvCs4 -
  //                          LatSet::InvCs2 * u2 * T(0.5));
  // }

  __any__ static void FirstOrder_Incompresible(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }
  __any__ static void SecondOrder(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
};

namespace equilibrium {

template <typename CELL>
struct SecondOrderImpl {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    const T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      const T uc = u * latset::c<LatSet>(k);
      feq[k] = latset::w<LatSet>(k) * rho *
      (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
    }
  }
};

template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d> &u, const T rho, const T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T{0.5});
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, rho, u);
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    SecondOrderImpl<CELL>::apply(feq, rho, u);
  }
};


template <typename CELL>
struct FirstOrderImpl {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
    }
  }
};

template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d> &u, const T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
  
    apply(feq, rho, u);

  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
      FirstOrderImpl<CELL>::apply(feq, rho, u);
  }
};

// Init feq
// common use: Init<moment::useFieldrhoU<CELL>, equilibrium::SecondOrder<CELL>>
template <typename MomentaScheme, typename EquilibriumScheme>
struct Init {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // macroscopic variables, usually from field value
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = feq[i];
    }
  }

};

template <typename CELLTYPE>
struct PhaseFieldEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  
  __any__ static inline void apply(std::array<T, LatSet::q>& feq, T phi, 
                                   const Vector<T, LatSet::d>& u) {
    T usqr = u[0]*u[0] + u[1]*u[1];
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T edotu = T{};
      for (unsigned int l = 0; l < LatSet::d; ++l) {
        edotu += latset::c<LatSet>(k)[l] * u[l];
      }
      feq[k] = latset::w<LatSet>(k) * phi * (T{1} + LatSet::InvCs2 * edotu);
    }
  }
  
  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& feq) {
    T phi = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      phi += cell[i];
    }
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, phi, u);
  }
};



}  // namespace equilibrium