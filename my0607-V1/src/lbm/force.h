#pragma once

#include "data_struct/Vector.h"
#include "lbm/lattice_set.h"
#include "utils/alias.h"

namespace force {

// F_v : viscous force
template <typename PFCELL, typename NSCELL>
struct VisForce {
  using T = typename NSCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;

  __any__ static void apply(PFCELL& pfcell, NSCELL& nscell) {
    const Vector<T, LatSet::d>& u = nscell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& grad = pfcell.template get<GRAD<T, LatSet::d>>();
    T dux = T{}, duy = T{}, dvx = T{}, dvy = T{};

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T wi = latset::w<LatSet>(i);
      const auto& ci = latset::c<LatSet>(i);
      auto ui = nscell.getNeighbor(i).template get<VELOCITY<T, LatSet::d>>()[0];
      auto vi = nscell.getNeighbor(i).template get<VELOCITY<T, LatSet::d>>()[1];
      dux += wi * ui * ci[0];
      duy += wi * ui * ci[1];
      dvx += wi * vi * ci[0];
      dvy += wi * vi * ci[1];
    }

    dux /= LatSet::cs2; duy /= LatSet::cs2;
    dvx /= LatSet::cs2; dvy /= LatSet::cs2;

    T omega = nscell.template get<OMEGA<T>>();
    T nu = LatSet::cs2 * (T(1) / omega - T(0.5));
    T rho_h = pfcell.template get<RHO_H<T>>();
    T rho_l = pfcell.template get<RHO_L<T>>();

    auto& F = nscell.template get<FORCE<T, LatSet::d>>();
    F[0] += nu * (rho_h - rho_l) * (T{2} * (dux * grad[0]) + (duy + dvx) * grad[1]);
    F[1] += nu * (rho_h - rho_l) * (T{2} * (dvy * grad[1]) + (duy + dvx) * grad[0]);
  }
};

// F_p : pressure force
template <typename PFCELL, typename NSCELL>
struct PressForce {
  using T = typename NSCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;

  __any__ static void apply(PFCELL& pfcell, NSCELL& nscell) {
    T press = nscell.template get<PRESSURE<T>>();
    T rho = nscell.template get<typename NSCELL::GenericRho>();
    T rho_h = pfcell.template get<RHO_H<T>>();
    T rho_l = pfcell.template get<RHO_L<T>>();
    const auto& grad = pfcell.template get<GRAD<T, LatSet::d>>();

    auto& F = nscell.template get<FORCE<T, LatSet::d>>();
    F -= (press * (rho_h - rho_l) / rho) * grad;
  }
};

// F_s : surface tension force
template <typename PFCELL, typename NSCELL>
struct SurfaceTension {
  using T = typename NSCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;

  __any__ static void apply(PFCELL& pfcell, NSCELL& nscell) {
    const auto& grad = pfcell.template get<GRAD<T, LatSet::d>>();
    T laplacian = T{0};
    T phi = pfcell.template get<typename PFCELL::GenericRho>();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T phi_i = pfcell.getNeighbor(i).template get<typename PFCELL::GenericRho>();
      T wi = latset::w<LatSet>(i);
      laplacian += wi * (phi_i - phi);
    }
    laplacian *= T{2} / LatSet::cs2;

    T beta  = pfcell.template get<BETA<T>>();
    T kappa = pfcell.template get<KAPPA<T>>();

    T chem = T{4} * beta * (phi - T{0}) * (phi - T{1}) * (phi - T{0.5}) - kappa * laplacian;
    nscell.template get<FORCE<T, LatSet::d>>() += chem * grad;
  }
};

// F_b : buoyancy/gravity force
template <typename PFCELL, typename NSCELL>
struct Gravity {
  using T = typename NSCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;

  __any__ static void apply(PFCELL& pfcell, NSCELL& nscell) {
    T phi = pfcell.template get<typename PFCELL::GenericRho>();
    T rho_l = pfcell.template get<RHO_L<T>>();
    T rho_h = pfcell.template get<RHO_H<T>>();
    T g = nscell.template get<CONSTFORCE<T, LatSet::d>>()[1];  // y-component of gravity
    g = -0.01;
    T rho = rho_l + phi * (rho_h - rho_l);
    nscell.template get<FORCE<T, LatSet::d>>()[1] += (rho - rho_h) * g;
  }
};

}  // namespace force
