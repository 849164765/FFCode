#pragma once

#include "ff/ff2d.h"

namespace ff {

template <typename CELL>
__any__ void FF3D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad;
  grad[0] = T{0};
  grad[1] = T{0};
  grad[2] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T phi_i = cell.getNeighbor(i).template get<GenericRho>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad[0] += wi * ci[0] * phi_i;
    grad[1] += wi * ci[1] * phi_i;
    grad[2] += wi * ci[2] * phi_i;
  }
  for (unsigned int d = 0; d < LatSet::d; ++d) {
    grad[d] /= LatSet::cs2;
  }
  cell.template get<GRAD<T, LatSet::d>>() = grad;

  Vector<T, LatSet::d> n;
  const T cutoff = T(0.02);
  const T eps = cutoff * cutoff;
  T mag2 = grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2];
  T inv_mag = T(1) / std::sqrt(mag2 + eps);
  for (unsigned int d = 0; d < LatSet::d; ++d) {
    n[d] = grad[d] * inv_mag;
  }
  cell.template get<NORMAL<T, LatSet::d>>() = n;
}

// ---- FFLaplacian3D ----
// ∇²φ = (2/(cs²*Δt²)) * Σ_{α=0}^{Q-1} w_α * [φ(x+e_α) - φ(x)]
// For lattice units (Δt=1, Δx=1): factor = 2/cs²
// Same formula as FFLaplacian2D, generic over LatSet::q.
template <typename CELL>
__any__ void FFLaplacian3D<CELL>::apply(CELL& cell) {
  T phi_self = cell.template get<GenericRho>();
  T laplacian = T{0};

  for (unsigned int k = 0; k < LatSet::q; ++k) {
    T phi_k = cell.getNeighbor(k).template get<GenericRho>();
    T wk = latset::w<LatSet>(k);
    // sum w_k * (phi_neighbor - phi_self)
    laplacian += wk * (phi_k - phi_self);
  }
  // factor 2/cs² for lattice units
  laplacian *= T{2} / LatSet::cs2;

  cell.template get<LAPLACIAN<T>>() = laplacian;
}

// ---- FFChemPotential3D ----
// λ = 4β * φ*(φ-1)*(φ-0.5) - κ * ∇²φ
// Same formula as FFChemPotential2D, no directional dependence.
template <typename CELL>
__any__ void FFChemPotential3D<CELL>::apply(CELL& cell) {
  T phi = cell.template get<GenericRho>();
  T laplacian = cell.template get<LAPLACIAN<T>>();
  T beta = cell.template get<BETA<T>>();
  T kappa = cell.template get<KAPPA<T>>();

  // φ*(φ-1)*(φ-0.5)
  T double_well = phi * (phi - T{1}) * (phi - T{0.5});
  T chem_potential = T{4} * beta * double_well - kappa * laplacian;

  cell.template get<CHEMICALPOTENTIAL<T>>() = chem_potential;
}

// ---- FFSurfaceTension3D ----
// F_s = λ * ∇φ  → added to ns_cell FORCE
template <typename PFCELL, typename NSCELL>
__any__ void FFSurfaceTension3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T chem_potential = pf_cell.template get<CHEMICALPOTENTIAL<T>>();
  const Vector<T, LatSet::d>& grad = pf_cell.template get<GRAD<T, LatSet::d>>();

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  for (unsigned int d = 0; d < LatSet::d; ++d) {
    ns_force[d] += chem_potential * grad[d];
  }
}

// ---- FFGravityForce3D ----
// Gravity acts in the z-direction (vertical, Nz is the domain height)
template <typename PFCELL, typename NSCELL>
__any__ void FFGravityForce3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T g = pf_cell.template get<GRAVITY<T>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  // z-direction is index 2 for D3Q19
  ns_cell.template get<FORCE<T, LatSet::d>>()[2] += -rho * g;
}

// ---- FFPreForce3D ----
// F_p = -(p/3) * DeltaRho * grad_phi
template <typename PFCELL, typename NSCELL>
__any__ void FFPreForce3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T p = ns_cell.template get<PRESSURE<T>>();
  T delta_rho = pf_cell.template get<DELTARHO<T>>();
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();

  T coeff = -p / T{3} * delta_rho;
  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  for (unsigned int d = 0; d < LatSet::d; ++d) {
    ns_force[d] += coeff * grad_phi[d];
  }
}

// ---- FFRhoOmegaUpdate3D ----
// WARNING: Cell::getOmega() returns the block-level scalar (from converter),
// NOT the per-cell OMEGA<T> field. So per-cell omega updates have NO effect on
// the collision. For variable viscosity, modify the collision operator to call
// getOmegaf() instead of getOmega(). See src/data_struct/cell.h:168.
//
// ρ = ρ_l + φ*(ρ_h - ρ_l)
// ν = η/ρ, τ = 0.5 + ν/cs², omega = 1/τ
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaUpdate3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    ns_cell.template get<DENSITY<T>>() = rho;
  } else {
    ns_cell.template get<typename NSCELL::GenericRho>() = rho;
  }

  T eta = eta_l + phi * (eta_h - eta_l);
  T nu = eta / rho;
  T tau = T{0.5} + nu / LatSet::cs2;
  T omega = T{1} / tau;
  // 上限从 1.95 提到 1.99: 允许低粘度流体 (KH 溶剂 eta~2.6e-4, RT 铁磁流体 eta~5e-3/rho=3)
  // 达到其本征 omega (1.98-1.997)。旧上限 1.95 会把这两类流体的有效运动粘度钳到
  // nu~0.00427 (Re_eff~300), 远低于论文意图 (Re=5000 / RT 铁磁流体 nu=0.00167),
  // 导致 KH 演化速率与 RT 尖钉下降速率偏离论文。1.99 -> tau=0.5025 仍稳定 (MRT 分离
  // 了体/幽灵矩 1.1-1.2, 剪切矩才绑 omega)。已核对 bubbleMag3d(omega<=1.93) 与
  // rosenMag3d(omega<=1.92) 不受影响。
  if (omega > T{1.99}) omega = T{1.99};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<OMEGA<T>>() = omega;
}

// ---- FFViscoForce3D (Palabos/Latt Simplified D3Q19 Basis) ----

// F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
// C_ab = Σ_k c_ka * c_kb * mgneq_k
// mgneq = InvM · S1 · (m - m_eq)   (first-pass MRT relaxation)
// S1 = diag(0, ω, ω, 0, s_q, 0, s_q, 0, s_q, ω, ω, ω, ω, ω, ω, ω, ω, ω, ω)
//
// Adds F_v to NS FORCE (which already contains F_s + F_b + F_p)

template <typename PFCELL, typename NSCELL>
__any__ void FFViscoForce3DM<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  // Read per-cell omega
  T omega{};
  if constexpr (ns_cell.template hasField<OMEGA<T>>()) {
    omega = ns_cell.template get<OMEGA<T>>();
  } else {
    omega = ns_cell.getOmega();
  }

  // Read actual density: DENSITY field (incompressible) or sum(f) (compressible)
  T f[LatSet::q];
  for (unsigned int k = 0; k < LatSet::q; ++k) f[k] = ns_cell[k];
  T rho;
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    rho = ns_cell.template get<DENSITY<T>>();
  } else {
    rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] 
        + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];
  }

  constexpr bool isIncompressible = NSCELL::template hasField<DENSITY<T>>();
  T ux, uy, uz;
  if constexpr (isIncompressible) {
    ux = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
    uy = f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
    uz = f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  } else {
    ux = (f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16]) / rho;
    uy = (f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18]) / rho;
    uz = (f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18]) / rho;
  }

  // Half-force correction using current FORCE (= F_s + F_b + F_p)
  const auto& F_in = ns_cell.template get<FORCE<T, LatSet::d>>();
  T halfInvRho = T{0.5} / rho;
  T ucx = ux + halfInvRho * F_in[0];
  T ucy = uy + halfInvRho * F_in[1];
  T ucz = uz + halfInvRho * F_in[2];

  // 严格遵循 Premnath (2007) Eq. 77 的物理绑定关系
  // 以及 Zu & He (2013) 对应力重构的渐近一致性要求
  T s_q = T{8} * (T{2} - omega) / (T{8} - omega); // 2D魔法公式在3D可作为近似，或直接设为1.0~1.2

  const T rtvec1[LatSet::q] {
    T{0},     // 0: rho (守恒)
    omega,    // 1: e (物理)
    omega,    // 2: epsilon (物理)
    T{0},     // 3: jx (守恒)
    s_q,      // 4: qx (能量通量)
    T{0},     // 5: jy (守恒)
    s_q,      // 6: qy 
    T{0},     // 7: jz (守恒)
    s_q,      // 8: qz 
    omega,    // 9: 剪切 (对应文献 s_10，必须等于 omega 以保证 Eq. 68 成立！)
    omega,    // 10: 幽灵矩 (保持物理一致性，绝不能是 1.5)
    omega,    // 11: 剪切 (对应文献 s_11，必须等于 omega)
    omega,    // 12: 幽灵矩 
    omega,    // 13: 剪切 (对应文献 s_12，必须等于 omega 以保证 Eq. 70 成立！)
    omega,    // 14: 剪切 (对应文献 s_13)
    omega,    // 15: 剪切 (对应文献 s_14)
    omega,    // 16: 幽灵矩
    omega,    // 17: 幽灵矩
    omega     // 18: 幽灵矩
  };

  // Moments from populations (M·f, unrolled D3Q19 matching library M matrix)
  T zeroth = rho; 
  T m_raw[LatSet::q];
  m_raw[0] = zeroth;
  m_raw[1] = -f[0] + f[7] + f[8] + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];
  m_raw[2] = f[0] - T{2}*(f[1]+f[2]+f[3]+f[4]+f[5]+f[6]) + f[7]+f[8]+f[9]+f[10]+f[11]+f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
  m_raw[3] = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
  m_raw[4] = -T{2}*f[1] + T{2}*f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
  m_raw[5] = f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
  m_raw[6] = -T{2}*f[3] + T{2}*f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
  m_raw[7] = f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  m_raw[8] = -T{2}*f[5] + T{2}*f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  m_raw[9] = T{2}*(f[1]+f[2]) - (f[3]+f[4]+f[5]+f[6]) + (f[7]+f[8]+f[9]+f[10]) - T{2}*(f[11]+f[12]) + (f[13]+f[14]+f[15]+f[16]) - T{2}*(f[17]+f[18]);
  m_raw[10] = -T{2}*(f[1]+f[2]) + (f[3]+f[4]+f[5]+f[6]) + (f[7]+f[8]+f[9]+f[10]) - T{2}*(f[11]+f[12]) + (f[13]+f[14]+f[15]+f[16]) - T{2}*(f[17]+f[18]);
  m_raw[11] = f[3]+f[4] - (f[5]+f[6]) + f[7]+f[8] - (f[9]+f[10]) + f[13]+f[14] - (f[15]+f[16]);
  m_raw[12] = -f[3]-f[4] + f[5]+f[6] + f[7]+f[8] - (f[9]+f[10]) + f[13]+f[14] - (f[15]+f[16]);
  m_raw[13] = f[7]+f[8] - f[13]-f[14];
  m_raw[14] = f[9]+f[10] - f[15]-f[16];
  m_raw[15] = f[11]+f[12] - f[17]-f[18];
  m_raw[16] = f[7]-f[8] - f[9]+f[10] + f[13]-f[14] - f[15]+f[16];
  m_raw[17] = f[7]-f[8] - f[11]+f[12] - f[13]+f[14] - f[17]+f[18];
  m_raw[18] = f[9]-f[10] - f[11]+f[12] - f[15]+f[16] + f[17]-f[18];

  // Equilibrium moments: Palabos/Latt Simplified D3Q19 Basis
  // NOTE: This basis naturally maintains Galilean invariance without O(u^3) corrections!
  // Many higher-order moments are EXACTLY ZERO in equilibrium.
  T ucx2 = ucx * ucx;
  T ucy2 = ucy * ucy;
  T ucz2 = ucz * ucz;
  T uc2 = ucx2 + ucy2 + ucz2;
  T m_eq[LatSet::q];
  
  if constexpr (isIncompressible) {
    m_eq[0] = zeroth;
    m_eq[1] = zeroth * uc2;
    m_eq[2] = T{0};
    m_eq[3] = ucx;
    m_eq[4] = T{0};
    m_eq[5] = ucy;
    m_eq[6] = T{0};
    m_eq[7] = ucz;
    m_eq[8] = T{0};
    m_eq[9] = T{2} * ucx2 - ucy2 - ucz2;
    m_eq[10] = T{0};
    m_eq[11] = ucy2 - ucz2;
    m_eq[12] = T{0};
    m_eq[13] = ucx * ucy;
    m_eq[14] = ucx * ucz;
    m_eq[15] = ucy * ucz;
    m_eq[16] = T{0}; m_eq[17] = T{0}; m_eq[18] = T{0};
  } else {
    m_eq[0] = rho;
    m_eq[1] = rho * uc2;
    m_eq[2] = T{0};
    m_eq[3] = rho * ucx;
    m_eq[4] = T{0};
    m_eq[5] = rho * ucy;
    m_eq[6] = T{0};
    m_eq[7] = rho * ucz;
    m_eq[8] = T{0};
    m_eq[9] = rho * (T{2} * ucx2 - ucy2 - ucz2);
    m_eq[10] = T{0};
    m_eq[11] = rho * (ucy2 - ucz2);
    m_eq[12] = T{0};
    m_eq[13] = rho * ucx * ucy;
    m_eq[14] = rho * ucx * ucz;
    m_eq[15] = rho * ucy * ucz;
    m_eq[16] = T{0}; m_eq[17] = T{0}; m_eq[18] = T{0};
  }

  // Deviation in moment space
  T dm[LatSet::q];
  for (unsigned int i = 0; i < LatSet::q; ++i) dm[i] = m_raw[i] - m_eq[i];

  // Non-equilibrium population: mgneq_i = Σ_j InvM(i,j) * S1(j) * dm_j
  T mgneq[LatSet::q] {};
  for (unsigned int i = 0; i < LatSet::q; ++i) {
    mgneq[i] = T{};
    for (unsigned int j = 0; j < LatSet::q; ++j) {
      mgneq[i] += mrt::InvM<LatSet>(i, j) * rtvec1[j] * dm[j];
    }
  }

  // C_ab = Σ_k c_ka * c_kb * mgneq_k (unrolled D3Q19 stress tensor)
  T Cxx = mgneq[1] + mgneq[2] + mgneq[7] + mgneq[8] + mgneq[9] + mgneq[10] + mgneq[13] + mgneq[14] + mgneq[15] + mgneq[16];
  T Cyy = mgneq[3] + mgneq[4] + mgneq[7] + mgneq[8] + mgneq[11] + mgneq[12] + mgneq[13] + mgneq[14] + mgneq[17] + mgneq[18];
  T Czz = mgneq[5] + mgneq[6] + mgneq[9] + mgneq[10] + mgneq[11] + mgneq[12] + mgneq[15] + mgneq[16] + mgneq[17] + mgneq[18];
  T Cxy = mgneq[7] + mgneq[8] - mgneq[13] - mgneq[14];
  T Cxz = mgneq[9] + mgneq[10] - mgneq[15] - mgneq[16];
  T Cyz = mgneq[11] + mgneq[12] - mgneq[17] - mgneq[18];

  // Read grad_phi and DeltaRho from PF cell
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T delta_rho = pf_cell.template get<DELTARHO<T>>();

  // mu = nu * rho = (tau-0.5)*cs²*rho = (1/omega - 0.5)*cs²*rho
  T invOmega = T{1} / omega;
  T nu = (invOmega - T{0.5}) * LatSet::cs2;
  T mu = nu * rho;

  // F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
  T prefactor = -T{3} * mu * delta_rho / rho;
  T Fv_x = prefactor * (Cxx * grad_phi[0] + Cxy * grad_phi[1] + Cxz * grad_phi[2]);
  T Fv_y = prefactor * (Cxy * grad_phi[0] + Cyy * grad_phi[1] + Cyz * grad_phi[2]);
  T Fv_z = prefactor * (Cxz * grad_phi[0] + Cyz * grad_phi[1] + Czz * grad_phi[2]);

  // Add F_v to NS FORCE
  ns_cell.template get<FORCE<T, LatSet::d>>()[0] += Fv_x;
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += Fv_y;
  ns_cell.template get<FORCE<T, LatSet::d>>()[2] += Fv_z;
}
// ---- FFViscoForce3D (SRT / BGK Version) ----

// F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
// C_ab = Σ_k c_ka * c_kb * f_neq_k
// f_neq_k = f_k - f_eq_k  (Standard BGK non-equilibrium)
//
// Adds F_v to NS FORCE (which already contains F_s + F_b + F_p)

template <typename PFCELL, typename NSCELL>
__any__ void FFViscoForce3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  // Read per-cell omega
  T omega{};
  if constexpr (ns_cell.template hasField<OMEGA<T>>()) {
    omega = ns_cell.template get<OMEGA<T>>();
  } else {
    omega = ns_cell.getOmega();
  }

  // Read actual density and populations
  T f[LatSet::q];
  for (unsigned int k = 0; k < LatSet::q; ++k) f[k] = ns_cell[k];
  T rho;
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    rho = ns_cell.template get<DENSITY<T>>();
  } else {
    rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] 
        + f[9] + f[10] + f[11] + f[12] + f[13] + f[14] + f[15] + f[16] + f[17] + f[18];
  }

  constexpr bool isIncompressible = NSCELL::template hasField<DENSITY<T>>();
  T ux, uy, uz;
  if constexpr (isIncompressible) {
    ux = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
    uy = f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
    uz = f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  } else {
    ux = (f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16]) / rho;
    uy = (f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18]) / rho;
    uz = (f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18]) / rho;
  }

  // Half-force correction using current FORCE
  const auto& F_in = ns_cell.template get<FORCE<T, LatSet::d>>();
  T halfInvRho = T{0.5} / rho;
  T ucx = ux + halfInvRho * F_in[0];
  T ucy = uy + halfInvRho * F_in[1];
  T ucz = uz + halfInvRho * F_in[2];

  // ---- SRT (BGK) Non-equilibrium distribution function ----
  // f_neq_k = f_k - f_eq_k
  // Standard D3Q19 equilibrium: f_eq = w * rho * (1 + 3(c·u) + 4.5(c·u)^2 - 1.5u^2)
  T f_neq[LatSet::q];
  T usqr = ucx * ucx + ucy * ucy + ucz * ucz;
  
  // Precompute common terms for equilibrium
  T w0 = T{1} / T{3};
  T w1 = T{1} / T{18};
  T w2 = T{1} / T{36};
  T base_eq = T{1} - T{1.5} * usqr;

  // 0: Rest node (0,0,0)
  f_neq[0] = f[0] - w0 * rho * base_eq;

  // 1,2: Face nodes (+-1, 0, 0) -> c·u = +-ux
  T cu_x = T{3} * ucx;
  f_neq[1] = f[1] - w1 * rho * (base_eq + cu_x + T{4.5} * ucx * ucx);
  f_neq[2] = f[2] - w1 * rho * (base_eq - cu_x + T{4.5} * ucx * ucx);

  // 3,4: Face nodes (0, +-1, 0) -> c·u = +-uy
  T cu_y = T{3} * ucy;
  f_neq[3] = f[3] - w1 * rho * (base_eq + cu_y + T{4.5} * ucy * ucy);
  f_neq[4] = f[4] - w1 * rho * (base_eq - cu_y + T{4.5} * ucy * ucy);

  // 5,6: Face nodes (0, 0, +-1) -> c·u = +-uz
  T cu_z = T{3} * ucz;
  f_neq[5] = f[5] - w1 * rho * (base_eq + cu_z + T{4.5} * ucz * ucz);
  f_neq[6] = f[6] - w1 * rho * (base_eq - cu_z + T{4.5} * ucz * ucz);

  // Edge nodes (diagonals)
  // Helper macro/lambda for edge nodes to keep code clean
  auto calc_edge = [&](T cx, T cy, T cz) {
    T cu = cx * ucx + cy * ucy + cz * ucz;
    return base_eq + T{3} * cu + T{4.5} * cu * cu;
  };

  // 7,8: (1,1,0), (-1,-1,0)
  f_neq[7]  = f[7]  - w2 * rho * calc_edge( 1,  1,  0);
  f_neq[8]  = f[8]  - w2 * rho * calc_edge(-1, -1,  0);

  // 9,10: (1,0,1), (-1,0,-1)
  f_neq[9]  = f[9]  - w2 * rho * calc_edge( 1,  0,  1);
  f_neq[10] = f[10] - w2 * rho * calc_edge(-1,  0, -1);

  // 11,12: (0,1,1), (0,-1,-1)
  f_neq[11] = f[11] - w2 * rho * calc_edge( 0,  1,  1);
  f_neq[12] = f[12] - w2 * rho * calc_edge( 0, -1, -1);

  // 13,14: (1,-1,0), (-1,1,0)
  f_neq[13] = f[13] - w2 * rho * calc_edge( 1, -1,  0);
  f_neq[14] = f[14] - w2 * rho * calc_edge(-1,  1,  0);

  // 15,16: (1,0,-1), (-1,0,1)
  f_neq[15] = f[15] - w2 * rho * calc_edge( 1,  0, -1);
  f_neq[16] = f[16] - w2 * rho * calc_edge(-1,  0,  1);

  // 17,18: (0,1,-1), (0,-1,1)
  f_neq[17] = f[17] - w2 * rho * calc_edge( 0,  1, -1);
  f_neq[18] = f[18] - w2 * rho * calc_edge( 0, -1,  1);


  // ---- Stress Tensor C_ab = Σ_k c_ka * c_kb * f_neq_k ----
  // Unrolled based on the specific D3Q19 velocity set mapping
  T Cxx = f_neq[1] + f_neq[2] + f_neq[7] + f_neq[8] + f_neq[9] + f_neq[10] 
        + f_neq[13] + f_neq[14] + f_neq[15] + f_neq[16];
        
  T Cyy = f_neq[3] + f_neq[4] + f_neq[7] + f_neq[8] + f_neq[11] + f_neq[12] 
        + f_neq[13] + f_neq[14] + f_neq[17] + f_neq[18];
        
  T Czz = f_neq[5] + f_neq[6] + f_neq[9] + f_neq[10] + f_neq[11] + f_neq[12] 
        + f_neq[15] + f_neq[16] + f_neq[17] + f_neq[18];

  T Cxy = f_neq[7] + f_neq[8] - f_neq[13] - f_neq[14];
  T Cxz = f_neq[9] + f_neq[10] - f_neq[15] - f_neq[16];
  T Cyz = f_neq[11] + f_neq[12] - f_neq[17] - f_neq[18];

  // Read grad_phi and DeltaRho from PF cell
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T delta_rho = pf_cell.template get<DELTARHO<T>>();

  // mu = nu * rho = (tau-0.5)*cs²*rho = (1/omega - 0.5)*cs²*rho
  T invOmega = T{1} / omega;
  T nu = (invOmega - T{0.5}) * LatSet::cs2;
  T mu = nu * rho;

  // F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
  T prefactor = -T{3} * mu * delta_rho / rho;
  T Fv_x = prefactor * (Cxx * grad_phi[0] + Cxy * grad_phi[1] + Cxz * grad_phi[2]);
  T Fv_y = prefactor * (Cxy * grad_phi[0] + Cyy * grad_phi[1] + Cyz * grad_phi[2]);
  T Fv_z = prefactor * (Cxz * grad_phi[0] + Cyz * grad_phi[1] + Czz * grad_phi[2]);

  // Add F_v to NS FORCE
  ns_cell.template get<FORCE<T, LatSet::d>>()[0] += Fv_x;
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += Fv_y;
  ns_cell.template get<FORCE<T, LatSet::d>>()[2] += Fv_z;
}

}  // namespace ff
