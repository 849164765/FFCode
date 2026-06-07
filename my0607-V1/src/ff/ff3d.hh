#pragma once

#include "ff/ff2d.h"

namespace ff {

template <typename CELL>
__any__ void FF3D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad;
  grad[0] = T{0};
  grad[1] = T{0};
  grad[2] = T{0};

  for (unsigned int k = 0; k < LatSet::q; ++k) {
    const auto& neighbor_pos = cell.getNeighbor(k);
    const auto& neighbor_neg = cell.getNeighbor(LatSet::opp[k]);
    T phi_pos = neighbor_pos.template get<GenericRho>();
    T phi_neg = neighbor_neg.template get<GenericRho>();

    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);

    grad[0] += wk * (phi_pos - phi_neg) * ck[0] / (2 * LatSet::cs2);
    grad[1] += wk * (phi_pos - phi_neg) * ck[1] / (2 * LatSet::cs2);
    grad[2] += wk * (phi_pos - phi_neg) * ck[2] / (2 * LatSet::cs2);
  }
  cell.template get<GRAD<T, LatSet::d>>() = grad;

  Vector<T, LatSet::d> n;
  T grad_mag = grad.getnorm();
  T epsilon = T(1e-10);
  if (grad_mag > epsilon) {
    n[0] = grad[0] / (grad_mag + epsilon);
    n[1] = grad[1] / (grad_mag + epsilon);
    n[2] = grad[2] / (grad_mag + epsilon);
  } else {
    n[0] = T{0};
    n[1] = T{0};
    n[2] = T{0};
  }
  cell.template get<NORMAL<T, LatSet::d>>() = n;
}

}  // namespace ff
