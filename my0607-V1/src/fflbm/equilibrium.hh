// equilibrium.hh — equilibrium distribution functions
//
// SecondOrderEq<CELL> : Navier-Stokes compressible
//   feq_k = w_k·ρ·(1 + c_k·u/cs² + (c_k·u)²/(2cs⁴) - u²/(2cs²))
//
// FirstOrderEq<CELL>  : phase-field (order parameter)
//   feq_k = w_k·φ·(1 + c_k·u/cs²)

#pragma once

namespace fflbm {

template <typename CELL>
struct SecondOrderEq {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q>& feq,
                            const T rho, const Vector<T, LatSet::d>& u) {
    const T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      const T uc = u * latset::c<LatSet>(k);
      feq[k] = latset::w<LatSet>(k) * rho *
               (T{1} + LatSet::InvCs2 * uc +
                uc * uc * T{0.5} * LatSet::InvCs4 -
                LatSet::InvCs2 * u2 * T{0.5});
    }
  }

  // read rho and u from cell fields
  __any__ static void apply(std::array<T, LatSet::q>& feq, const CELL& cell) {
    const T rho = cell.template get<GenericRho>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, rho, u);
  }
};

template <typename CELL>
struct FirstOrderEq {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q>& feq,
                            const T phi, const Vector<T, LatSet::d>& u) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = latset::w<LatSet>(k) * phi *
               (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
    }
  }

  // read phi and u from cell fields
  __any__ static void apply(std::array<T, LatSet::q>& feq, const CELL& cell) {
    const T phi = cell.template get<GenericRho>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, phi, u);
  }
};

// ---------- Pressure-based equilibrium (g-LBE) ----------
// g_i^eq = w_i·[p̃ + c_i·u/cs² + (c_i·u)²/(2cs⁴) - u²/(2cs²)]
// where p̃ = p/cs² (pseudo-density), u is force-corrected velocity.
// Σg_i^eq = p̃, Σc_i·g_i^eq = u.
template <typename CELL>
struct PressureEq {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;

  __any__ static void apply(std::array<T, LatSet::q>& geq,
                            const T p_cs2, const Vector<T, LatSet::d>& u) {
    const T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      const T uc = u * latset::c<LatSet>(k);
      geq[k] = latset::w<LatSet>(k) *
               (p_cs2 + LatSet::InvCs2 * uc +
                uc * uc * T{0.5} * LatSet::InvCs4 -
                u2 * T{0.5} * LatSet::InvCs2);
    }
  }
};

}  // namespace fflbm
