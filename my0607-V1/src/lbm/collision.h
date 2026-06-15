// collision.h
// Lattice Boltzmann collision operators and equilibrium distribution functions

#pragma once

#include "lbm/lattice_set.h"
#include "utils/alias.h"

namespace equilibrium {

// First-order equilibrium for Allen-Cahn phase field:
//   f_eq(k) = w_k * phi * (1 + c_k·u / cs²)
template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int D = LatSet::d;

  static T get(CELL& cell, unsigned int k) {
    const T phi = cell.template get<PHI<T>>();
    const auto& u = cell.template get<VELOCITY<T, D>>();
    T cu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k) * phi * (T(1) + LatSet::InvCs2 * cu);
  }
};

// Second-order equilibrium for Navier-Stokes flow field (D2Q9 coefficients):
//   g_eq(k) = w_k * rho * (1 + 3(c·u) + 4.5(c·u)² - 1.5u²)
template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int D = LatSet::d;

  static T get(CELL& cell, unsigned int k) {
    const T rho = cell.template get<RHO<T>>();
    const auto& u = cell.template get<VELOCITY<T, D>>();
    T cu = latset::c<LatSet>(k) * u;
    T u2 = u * u;
    return latset::w<LatSet>(k) * rho
           * (T(1) + T(3) * cu + T(4.5) * cu * cu - T(1.5) * u2);
  }
};

}  // namespace equilibrium

namespace collision {

// BGK collision operator with optional source term (Guo forcing scheme).
//
// Template parameters:
//   Equilibrium  — equilibrium distribution functor (FirstOrder or SecondOrder)
//   SourceField  — field type for the source direction vector
//                   (use void when hasSource=false)
//   hasSource    — enable (true) or disable (false) source term
//
// Without source:
//   f_k' = f_k - omega * (f_k - f_k^eq)
//
// With source:
//   S_k = (1 - 0.5*omega) * w_k * (c_k · n)
//   f_k' = f_k - omega * (f_k - f_k^eq) + S_k
template <typename Equilibrium, typename SourceField, bool hasSource>
struct BGKSource {
  using T = typename Equilibrium::T;
  using LatSet = typename Equilibrium::LatSet;
  static constexpr unsigned int Q = LatSet::q;

  template <typename CELL>
  static void apply(CELL& cell) {
    const T omega = cell.getOmega();

    if constexpr (hasSource) {
      const T omega_half = T(0.5) * omega;
      const auto& source = cell.template get<SourceField>();
      for (unsigned int k = 0; k < Q; ++k) {
        T feq = Equilibrium::get(cell, k);
        T f_old = cell[k];
        T Sk = (T(1) - omega_half)
               * latset::w<LatSet>(k)
               * (latset::c<LatSet>(k) * source);
        cell[k] = f_old - omega * (f_old - feq) + Sk;
      }
    } else {
      for (unsigned int k = 0; k < Q; ++k) {
        T feq = Equilibrium::get(cell, k);
        cell[k] = cell[k] - omega * (cell[k] - feq);
      }
    }
  }
};

}  // namespace collision
