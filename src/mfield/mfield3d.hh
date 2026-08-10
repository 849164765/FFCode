#pragma once

#include "mfield/mfield3d.h"

namespace mfield {

template <typename PFCELL, typename MFCELL>
__any__ void MFUpdateCoeffs3D<PFCELL, MFCELL>::apply(PFCELL& pf_cell,
                                                     MFCELL& mf_cell) {
  using T = typename PFCELL::FloatType;

  const T phi = pf_cell.template get<PHI<T>>();
  const T mu_l = mf_cell.template get<MU_L<T>>();
  const T mu_h = mf_cell.template get<MU_H<T>>();
  const T chi_l = mf_cell.template get<CHI_L<T>>();
  const T chi_h = mf_cell.template get<CHI_H<T>>();
  const T mu = mu_l + phi * (mu_h - mu_l);
  const T chi = chi_l + phi * (chi_h - chi_l);

  // K changes only the convergence rate of the fixed-point solve.  The
  // continuum equation ∇·(mu∇psi)=0 is unchanged.
  const T solver_k = mf_cell.template get<PSI_K<T>>();
  const T tau_psi = T{0.5} + solver_k * mu;
  T omega_psi = T{1} / tau_psi;
  if (omega_psi > T{1.95}) omega_psi = T{1.95};
  if (omega_psi < T{0.01}) omega_psi = T{0.01};

  mf_cell.template get<MU_PERCELL<T>>() = mu;
  mf_cell.template get<CHI_PERCELL<T>>() = chi;
  mf_cell.template get<OMEGA_PSI<T>>() = omega_psi;
}

template <typename CELL>
__any__ void MFComputeH3D<CELL>::apply(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // D3Q7 has Q2 = sum_k w_k c_{k,a}^2 = 1/4, which differs from cs2=1/3.
  T q2 = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    q2 += wk * ck[0] * ck[0];
  }

  Vector<T, 3> grad{0, 0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const T psi_k = cell.getNeighbor(k).template get<PSI<T>>();
    const T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    for (unsigned int d = 0; d < 3; ++d) grad[d] += wk * ck[d] * psi_k;
  }
  for (unsigned int d = 0; d < 3; ++d) grad[d] /= q2;

  const T hx = -grad[0];
  const T hy = -grad[1];
  const T hz = -grad[2];
  cell.template get<HX<T>>() = hx;
  cell.template get<HY<T>>() = hy;
  cell.template get<HZ<T>>() = hz;
  cell.template get<HMAG<T>>() = std::sqrt(hx * hx + hy * hy + hz * hz);
}

template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce3D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;

  const T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  const T chi_l = mf_cell.template get<CHI_L<T>>();
  const T chi_h = mf_cell.template get<CHI_H<T>>();
  const T chi = chi_l + phi * (chi_h - chi_l);
  const T hmag = mf_cell.template get<HMAG<T>>();

  T q2 = T{0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    q2 += wk * ck[0] * ck[0];
  }

  Vector<T, 3> grad_hmag{0, 0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    const T hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    const T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    for (unsigned int d = 0; d < 3; ++d) {
      grad_hmag[d] += wk * ck[d] * hmag_k;
    }
  }
  for (unsigned int d = 0; d < 3; ++d) grad_hmag[d] /= q2;

  auto& force = ns_cell.template get<FORCE<T, 3>>();
  for (unsigned int d = 0; d < 3; ++d) {
    force[d] += chi * hmag * grad_hmag[d];
  }
}

}  // namespace mfield
