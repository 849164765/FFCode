#pragma once

#include "mfield/mfield2d.h"

namespace mfield {

// ---- MFUpdateCoeffs2D ----
// Update per-cell mu, chi, omega_psi from phi
// mu  = mu_l  + phi * (mu_h  - mu_l)
// chi = chi_l + phi * (chi_h - chi_l)
// tau_psi = 0.5 + epsilon * mu / cs²,  omega_psi = 1 / tau_psi (clamped)
//
// 论文式(37): (1/ε) ∂ψ/∂t = ∇·(μ∇ψ)
// ε 为伪时间加速参数: 大值 → 有效扩散系数 D_eff = ε*μ 增大 → 加速收敛到稳态
// 论文式(42): tau = 0.5 + 2.5*ε*μ/cs² (D2Q5 w=0.2, cs²=4/5)
// 本代码 D2Q5 w={1/3,1/6,...}, cs²=1/3, 对应 tau = 0.5 + ε*μ/cs²
// ε=1 时 D=μ (收敛慢, RT算例需~100万步); ε=10 时 D=10μ (~1万步收敛)
template <typename PFCELL, typename MFCELL>
__any__ void MFUpdateCoeffs2D<PFCELL, MFCELL>::apply(PFCELL& pf_cell,
                                                     MFCELL& mf_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;

  T phi   = pf_cell.template get<PHI<T>>();
  T mu_l  = mf_cell.template get<MU_L<T>>();
  T mu_h  = mf_cell.template get<MU_H<T>>();
  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T epsilon = mf_cell.template get<MF_EPSILON<T>>();

  T mu  = mu_l  + phi * (mu_h  - mu_l);
  T chi = chi_l + phi * (chi_h - chi_l);

  // tau = 0.5 + ε*μ / cs²  (有效扩散系数 D_eff = ε*μ)
  T tau_psi = T{0.5} + epsilon * mu / LatSet::cs2;
  T omega_psi = T{1} / tau_psi;
  if (omega_psi > T{1.95}) omega_psi = T{1.95};
  if (omega_psi < T{0.01}) omega_psi = T{0.01};

  mf_cell.template get<MU_PERCELL<T>>()   = mu;
  mf_cell.template get<CHI_PERCELL<T>>()  = chi;
  mf_cell.template get<OMEGA_PSI<T>>()    = omega_psi;
}

// ---- MFComputeH2D ----
// H = -∇ψ using D2Q5 4-direction isotropic gradient
// ∇ψ = Σ w_k * c_k * ψ_k / cs²
// D2Q5: c_k ∈ {(0,0),(1,0),(-1,0),(0,1),(0,-1)}
template <typename CELL>
__any__ void MFComputeH2D<CELL>::apply(CELL& cell) {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  T grad_x = T{0};
  T grad_y = T{0};

  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T psi_k = cell.getNeighbor(k).template get<PSI<T>>();
    T wk    = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_x += wk * ck[0] * psi_k;
    grad_y += wk * ck[1] * psi_k;
  }
  grad_x /= LatSet::cs2;
  grad_y /= LatSet::cs2;

  T Hx = -grad_x;
  T Hy = -grad_y;
  T Hmag = std::sqrt(Hx * Hx + Hy * Hy);

  cell.template get<HX<T>>()   = Hx;
  cell.template get<HY<T>>()   = Hy;
  cell.template get<HMAG<T>>() = Hmag;
}

// ---- MFMagneticForce2D ----
// Magnetic body force (Guo et al., Phys. Fluids 37, 022148 (2025), Eq. (8)):
//   F_m = (μ₀χ/2) ∇(|H|²) = μ₀χ|H|∇|H|   (with μ₀ ≡ 1)
//
// This form assumes linear magnetization (χ constant within each phase,
// interpolated by φ across the interface). The ∇χ (magnetostrictive) term
// is NOT included, following the reference.
template <typename PFCELL, typename MFCELL, typename NSCELL>
__any__ void MFMagneticForce2D<PFCELL, MFCELL, NSCELL>::apply(
    PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5 for neighbor access

  T chi_l = mf_cell.template get<CHI_L<T>>();
  T chi_h = mf_cell.template get<CHI_H<T>>();
  T Hmag  = mf_cell.template get<HMAG<T>>();
  T phi   = pf_cell.template get<typename PFCELL::GenericRho>();
  T chi   = chi_l + phi * (chi_h - chi_l);

  // ∇|H| via D2Q5 4-direction isotropic gradient
  Vector<T, 2> grad_Hmag{0, 0};
  for (unsigned int k = 1; k < LatSet::q; ++k) {
    T Hmag_k = mf_cell.getNeighbor(k).template get<HMAG<T>>();
    T wk = latset::w<LatSet>(k);
    const auto& ck = latset::c<LatSet>(k);
    grad_Hmag[0] += wk * ck[0] * Hmag_k;
    grad_Hmag[1] += wk * ck[1] * Hmag_k;
  }
  grad_Hmag[0] /= LatSet::cs2;
  grad_Hmag[1] /= LatSet::cs2;

  // F_mag = χ·|H|·∇|H|  (Guo 2025, Eq. (8), μ₀ ≡ 1)
  T Fmag_x = chi * Hmag * grad_Hmag[0];
  T Fmag_y = chi * Hmag * grad_Hmag[1];

  auto& ns_force = ns_cell.template get<FORCE<T, 2>>();
  ns_force[0] += Fmag_x;
  ns_force[1] += Fmag_y;
}

}  // namespace mfield
