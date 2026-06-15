// force.hh — Guo force scheme and buoyancy coupling
//
// GuoForce<T, LatSet>
//   F_i = w_i·F·[(c_i-u)/cs² + (c_i·u)·c_i/cs⁴]  (Guo 2002)
//
// Buoyancy<PFCELL, NSCELL>
//   F_b[1] = (ρ(φ) - ρ_h)·g
//   ρ(φ) = ρ_l + φ·(ρ_h - ρ_l)

#pragma once

#include "ff/ff2d.h"

namespace fflbm {

template <typename T, typename LatSet>
struct GuoForce {
  __any__ static void compute(std::array<T, LatSet::q>& Fi,
                              const Vector<T, LatSet::d>& u,
                              const Vector<T, LatSet::d>& F) {
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const auto& ci = latset::c<LatSet>(i);
      const T ci_dot_u = ci * u;
      Fi[i] = latset::w<LatSet>(i) * F *
              ((ci - u) * LatSet::InvCs2 +
               ci_dot_u * LatSet::InvCs4 * ci);
    }
  }
};

template <typename PFCELL, typename NSCELL>
struct Buoyancy {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;

  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell) {
    const T phi = pf_cell.template get<typename PFCELL::GenericRho>();
    const T rho_l = pf_cell.template get<ff::RHO_L<T>>();
    const T rho_h = pf_cell.template get<ff::RHO_H<T>>();
    const T g = pf_cell.template get<ff::GRAVITY<T>>();
    const T rho = rho_l + phi * (rho_h - rho_l);
    auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
    ns_force[1] += (rho - rho_h) * g;
  }
};

}  // namespace fflbm
