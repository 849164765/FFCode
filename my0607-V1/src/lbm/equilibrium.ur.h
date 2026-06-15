// equilibrium.ur.h
// Equilibrium distribution function helpers and constants
// Used by: bounce_back_boundary.h, cell.h, block_lattice.hh

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"

namespace equilibrium {

// ================================================================
// Constants: cs^2 = 1/3, cs^4 = 1/9  (for all lattices in this codebase)
// ================================================================
template <typename T>
constexpr T cs2 = T(1) / T(3);

template <typename T>
constexpr T cs4 = T(1) / T(9);

template <typename T>
constexpr T invCs2 = T(3);

// ================================================================
// D2Q9 omega weights  (omega_alpha, Eq.18 / Eq.31)
// omega0=4/9, omega1-4=1/9, omega5-8=1/36
// ================================================================
template <typename T>
constexpr T omega_D2Q9[9] = {
  T(4)/T(9),  T(1)/T(9),  T(1)/T(9),  T(1)/T(9),  T(1)/T(9),
  T(1)/T(36), T(1)/T(36), T(1)/T(36), T(1)/T(36)
};

// ================================================================
// D2Q5 omega weights  (Eq.40: omega_alpha = 1/5)
// ================================================================
template <typename T>
constexpr T omega_D2Q5[5] = {
  T(1)/T(5), T(1)/T(5), T(1)/T(5), T(1)/T(5), T(1)/T(5)
};

// ================================================================
// FirstOrder equilibrium (phase-field / advection-diffusion type)
// Eq.18: f_eq(k) = w_k * rho * (1 + e_k·u / cs^2)   with cs^2=1/3
// ================================================================
template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static T get(int k, const Vector<T, LatSet::d>& u, T rho) {
    T eu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k) * rho * (T(1) + LatSet::InvCs2 * eu);
  }
};

// ================================================================
// SecondOrder equilibrium (standard D2Q9)
// f_eq(k) = w_k * rho * [1 + e_k·u/cs^2 + (e_k·u)^2/(2*cs^4) - u^2/(2*cs^2)]
// ================================================================
template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr T invCs2 = LatSet::InvCs2;
  static constexpr T invCs4 = LatSet::InvCs4;

  static T get(int k, const Vector<T, LatSet::d>& u, T rho, T uSqr) {
    T eu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k) * rho
         * (T(1) + invCs2 * eu + invCs4 * eu * eu / T(2) - invCs2 * uSqr / T(2));
  }

  // Convenience: fill entire feq array
  static void apply(std::array<T, LatSet::q>& feq, T rho,
                    const Vector<T, LatSet::d>& u) {
    T uSqr = u * u;
    for (int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, rho, uSqr);
    }
  }
};

// ================================================================
// VelocityBasedNS equilibrium (Guo et al. Eq.31)
// g_eq(k) = w_k * [p/(rho*cs^2) + e_k·u/cs^2 + (e_k·u)^2/(2*cs^4) - u^2/(2*cs^2)]
// Note: zero-moment is p/(rho*cs^2), NOT rho
// ================================================================
template <typename CELL>
struct VelocityBasedNS {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr T cs2 = LatSet::cs2;
  static constexpr T invCs2 = LatSet::InvCs2;
  static constexpr T invCs4 = LatSet::InvCs4;

  static T get(int k, const Vector<T, LatSet::d>& u, T p_over_rho, T uSqr) {
    T eu = latset::c<LatSet>(k) * u;
    return latset::w<LatSet>(k)
         * (p_over_rho * invCs2 + invCs2 * eu + invCs4 * eu * eu / T(2) - invCs2 * uSqr / T(2));
  }
};

}  // namespace equilibrium
