#pragma once

#include "ff/ff2d.h"

namespace ff {

template <typename CELL>
__any__ void FF3D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad;
  grad[0] = T{0};
  grad[1] = T{0};
  grad[2] = T{0};

  const T phi_self = cell.template get<GenericRho>();
  const std::size_t N = cell.getN();
  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T phi_i;
    // boundary check: out-of-range neighbor index (ghost/overlap cells)
    // falls back to self phi (zero-gradient contribution)
    if (cell.getNeighborId(i) < N) {
      phi_i = cell.getNeighbor(i).template get<GenericRho>();
    } else {
      phi_i = phi_self;
    }
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

  const std::size_t N = cell.getN();
  for (unsigned int k = 0; k < LatSet::q; ++k) {
    T phi_k;
    // boundary check: out-of-range neighbor index (ghost/overlap cells)
    // falls back to self phi (zero-laplacian contribution)
    if (cell.getNeighborId(k) < N) {
      phi_k = cell.getNeighbor(k).template get<GenericRho>();
    } else {
      phi_k = phi_self;
    }
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
// NOTE: p = Σf (PRESSURE field) forms a positive feedback loop:
//   p → F_p → collision (Guo force) → non-equilibrium pops → streaming → p grows.
// Under high density ratio (Δρ up to 999), the loop gain exceeds 1 and p grows
// exponentially, leading to NaN. Clamp p to break the loop.
// Evidence: p=2.83e-3 stable (step 80), p=0.0262 unstable (step 81) → p_max=0.01.
template <typename PFCELL, typename NSCELL>
__any__ void FFPreForce3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T p = ns_cell.template get<PRESSURE<T>>();
  // NaN/Inf guard and clamp to break positive feedback loop.
  if (!std::isfinite(p)) p = T{0};
  constexpr T p_max = T{0.01};
  if (p < -p_max) p = -p_max;
  if (p > p_max) p = p_max;
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
  // Clamp rho to [rho_l, rho_h] to prevent out-of-range values from
  // corrupted phi amplifying downstream divisions (e.g. F/rho in MRTForce).
  if (!std::isfinite(rho) || rho < rho_l) rho = rho_l;
  if (rho > rho_h) rho = rho_h;
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    ns_cell.template get<DENSITY<T>>() = rho;
  } else {
    ns_cell.template get<typename NSCELL::GenericRho>() = rho;
  }

  T eta = eta_l + phi * (eta_h - eta_l);
  T nu = eta / rho;
  T tau = T{0.5} + nu / LatSet::cs2;
  T omega = T{1} / tau;
  if (omega > T{1.95}) omega = T{1.95};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<OMEGA<T>>() = omega;
}

// ---- FFViscoForce3D ----
// F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
// C_ab = Σ_k c_ka * c_kb * mgneq_k
// mgneq = InvM · S1 · (m - m_eq)
// S1: conserved moments 0, energy modes omega, q-modes s_q, stress modes omega.
// Hand-written D3Q19 version, analogous to FFViscoForce2D.
template <typename PFCELL, typename NSCELL>
__any__ void FFViscoForce3D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  static_assert(LatSet::q == 19 && LatSet::d == 3,
                "FFViscoForce3D only supports D3Q19");

  // Read per-cell omega
  T omega{};
  if constexpr (ns_cell.template hasField<OMEGA<T>>()) {
    omega = ns_cell.template get<OMEGA<T>>();
  } else {
    omega = ns_cell.getOmega();
  }

  // Read populations and density
  T f[19];
  for (unsigned int k = 0; k < 19; ++k) f[k] = ns_cell[k];

  // zeroth = Σf (zeroth moment of populations, used for equilibrium consistency)
  // rho   = physical density from DENSITY field (from phi, used for prefactor only)
  // NOTE: If zeroth were set to DENSITY, m_raw[0]=Σf would mismatch m_eq[0]=rho,
  // producing huge dm[0] = Σf - DENSITY (e.g. 1 - 1000 = -999) and exploding F_v.
  T zeroth = T{0};
  for (unsigned int k = 0; k < 19; ++k) zeroth += f[k];
  T rho{};
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    rho = ns_cell.template get<DENSITY<T>>();
  } else {
    rho = zeroth;
  }

  // Compute velocity: sum(c*f) / rho for compressible, sum(c*f) for incompressible
  // D3Q19 ordering (from lattice_set.h):
  // 0:rest, 1:+x, 2:-x, 3:+y, 4:-y, 5:+z, 6:-z,
  // 7:+x+y, 8:-x-y, 9:+x+z, 10:-x-z, 11:+y+z, 12:-y-z,
  // 13:+x-y, 14:-x+y, 15:+x-z, 16:-x+z, 17:+y-z, 18:-y+z
  constexpr bool isIncompressible = NSCELL::template hasField<DENSITY<T>>();
  T ux = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
  T uy = f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
  T uz = f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  if constexpr (!isIncompressible) {
    ux /= rho;
    uy /= rho;
    uz /= rho;
  }

  // Half-force correction using current FORCE (= F_s + F_b + F_p)
  const auto& F_in = ns_cell.template get<FORCE<T, LatSet::d>>();
  T halfInvRho = T{0.5} / rho;
  T ucx = ux + halfInvRho * F_in[0];
  T ucy = uy + halfInvRho * F_in[1];
  T ucz = uz + halfInvRho * F_in[2];

  // ---- First-pass MRT relaxation vector (D3Q19) ----
  T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
  const T rtvec1[19] {
    T{0},     // 0:  rho (conserved)
    omega,    // 1:  e
    omega,    // 2:  epsilon
    T{0},     // 3:  jx (conserved)
    s_q,      // 4:  qx
    T{0},     // 5:  jy (conserved)
    s_q,      // 6:  qy
    T{0},     // 7:  jz (conserved)
    s_q,      // 8:  qz
    omega,    // 9:  stress
    omega,    // 10: stress
    omega,    // 11: shear
    omega,    // 12: higher-order
    omega,    // 13: shear
    omega,    // 14: stress
    omega,    // 15: shear
    omega,    // 16: higher-order
    omega,    // 17: higher-order
    omega     // 18: higher-order
  };

  // Equilibrium populations (standard D3Q19)
  T feq[19];
  T ucx2 = ucx * ucx;
  T ucy2 = ucy * ucy;
  T ucz2 = ucz * ucz;
  T uc2 = ucx2 + ucy2 + ucz2;
  T rho_eq = isIncompressible ? zeroth : rho;
  // rest direction
  feq[0] = latset::w<LatSet>(0) * (rho_eq - T{1.5} * uc2);
  for (unsigned int k = 1; k < 19; ++k) {
    const auto& ck = latset::c<LatSet>(k);
    T cu = ck[0] * ucx + ck[1] * ucy + ck[2] * ucz;
    T wk = latset::w<LatSet>(k);
    feq[k] = wk * (rho_eq + T{3} * cu + T{4.5} * cu * cu - T{1.5} * uc2);
  }

  // ---- Moments from populations (handwritten M·f for D3Q19) ----
  T m_raw[19];
  T sum1_6 = f[1] + f[2] + f[3] + f[4] + f[5] + f[6];
  T sum7_10 = f[7] + f[8] + f[9] + f[10];
  T sum11_14 = f[11] + f[12] + f[13] + f[14];
  T sum15_18 = f[15] + f[16] + f[17] + f[18];
  T sum7_18 = sum7_10 + sum11_14 + sum15_18;
  m_raw[0] = f[0] + sum1_6 + sum7_18;
  m_raw[1] = -f[0] + sum7_18;
  m_raw[2] = f[0] - T{2} * sum1_6 + sum7_18;
  m_raw[3] = f[1] - f[2] + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
  m_raw[4] = T{-2} * (f[1] - f[2]) + f[7] - f[8] + f[9] - f[10] + f[13] - f[14] + f[15] - f[16];
  m_raw[5] = f[3] - f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
  m_raw[6] = -T{2} * f[3] + T{2} * f[4] + f[7] - f[8] + f[11] - f[12] - f[13] + f[14] + f[17] - f[18];
  m_raw[7] = f[5] - f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  m_raw[8] = -T{2} * f[5] + T{2} * f[6] + f[9] - f[10] + f[11] - f[12] - f[15] + f[16] - f[17] + f[18];
  m_raw[9] = f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8] - f[9] - f[10] - f[11] - f[12] + f[13] + f[14] + f[15] + f[16] - f[17] - f[18];
  m_raw[10] = T{-2} * (f[1] + f[2]) + (f[3] + f[4] + f[5] + f[6])
            + (f[7] + f[8] + f[9] + f[10])
            - T{2} * (f[11] + f[12])
            + (f[13] + f[14] + f[15] + f[16])
            - T{2} * (f[17] + f[18]);
  m_raw[11] = (f[3] + f[4]) - (f[5] + f[6])
            + (f[7] + f[8]) - (f[9] + f[10])
            + (f[13] + f[14]) - (f[15] + f[16]);
  m_raw[12] = -(f[3] + f[4]) + (f[5] + f[6])
            + (f[7] + f[8]) - (f[9] + f[10])
            + (f[13] + f[14]) - (f[15] + f[16]);
  m_raw[13] = (f[7] + f[8]) - (f[13] + f[14]);
  m_raw[14] = (f[9] + f[10]) - (f[15] + f[16]);
  m_raw[15] = (f[11] + f[12]) - (f[17] + f[18]);
  m_raw[16] = (f[7] - f[8]) - (f[9] - f[10])
            + (f[13] - f[14]) - (f[15] - f[16]);
  m_raw[17] = (f[7] - f[8]) - (f[11] - f[12])
            - (f[13] - f[14]) - (f[17] - f[18]);
  m_raw[18] = (f[9] - f[10]) - (f[11] - f[12])
            - (f[15] - f[16]) + (f[17] - f[18]);

  // Equilibrium moments (same M applied to feq)
  T m_eq[19];
  T esum1_6 = feq[1] + feq[2] + feq[3] + feq[4] + feq[5] + feq[6];
  T esum7_10 = feq[7] + feq[8] + feq[9] + feq[10];
  T esum11_14 = feq[11] + feq[12] + feq[13] + feq[14];
  T esum15_18 = feq[15] + feq[16] + feq[17] + feq[18];
  T esum7_18 = esum7_10 + esum11_14 + esum15_18;
  m_eq[0] = feq[0] + esum1_6 + esum7_18;
  m_eq[1] = -feq[0] + esum7_18;
  m_eq[2] = feq[0] - T{2} * esum1_6 + esum7_18;
  m_eq[3] = feq[1] - feq[2] + feq[7] - feq[8] + feq[9] - feq[10] + feq[13] - feq[14] + feq[15] - feq[16];
  m_eq[4] = T{-2} * (feq[1] - feq[2]) + feq[7] - feq[8] + feq[9] - feq[10] + feq[13] - feq[14] + feq[15] - feq[16];
  m_eq[5] = feq[3] - feq[4] + feq[7] - feq[8] + feq[11] - feq[12] - feq[13] + feq[14] + feq[17] - feq[18];
  m_eq[6] = -T{2} * feq[3] + T{2} * feq[4] + feq[7] - feq[8] + feq[11] - feq[12] - feq[13] + feq[14] + feq[17] - feq[18];
  m_eq[7] = feq[5] - feq[6] + feq[9] - feq[10] + feq[11] - feq[12] - feq[15] + feq[16] - feq[17] + feq[18];
  m_eq[8] = -T{2} * feq[5] + T{2} * feq[6] + feq[9] - feq[10] + feq[11] - feq[12] - feq[15] + feq[16] - feq[17] + feq[18];
  m_eq[9] = feq[1] + feq[2] + feq[3] + feq[4] + feq[5] + feq[6] + feq[7] + feq[8] - feq[9] - feq[10] - feq[11] - feq[12] + feq[13] + feq[14] + feq[15] + feq[16] - feq[17] - feq[18];
  m_eq[10] = T{-2} * (feq[1] + feq[2]) + (feq[3] + feq[4] + feq[5] + feq[6])
           + (feq[7] + feq[8] + feq[9] + feq[10])
           - T{2} * (feq[11] + feq[12])
           + (feq[13] + feq[14] + feq[15] + feq[16])
           - T{2} * (feq[17] + feq[18]);
  m_eq[11] = (feq[3] + feq[4]) - (feq[5] + feq[6])
           + (feq[7] + feq[8]) - (feq[9] + feq[10])
           + (feq[13] + feq[14]) - (feq[15] + feq[16]);
  m_eq[12] = -(feq[3] + feq[4]) + (feq[5] + feq[6])
           + (feq[7] + feq[8]) - (feq[9] + feq[10])
           + (feq[13] + feq[14]) - (feq[15] + feq[16]);
  m_eq[13] = (feq[7] + feq[8]) - (feq[13] + feq[14]);
  m_eq[14] = (feq[9] + feq[10]) - (feq[15] + feq[16]);
  m_eq[15] = (feq[11] + feq[12]) - (feq[17] + feq[18]);
  m_eq[16] = (feq[7] - feq[8]) - (feq[9] - feq[10])
           + (feq[13] - feq[14]) - (feq[15] - feq[16]);
  m_eq[17] = (feq[7] - feq[8]) - (feq[11] - feq[12])
           - (feq[13] - feq[14]) - (feq[17] - feq[18]);
  m_eq[18] = (feq[9] - feq[10]) - (feq[11] - feq[12])
           - (feq[15] - feq[16]) + (feq[17] - feq[18]);

  // Deviation in moment space
  T dm[19];
  for (unsigned int i = 0; i < 19; ++i) dm[i] = m_raw[i] - m_eq[i];

  // C_ab = Σ_k c_ka * c_kb * (InvM · (rtvec1 * dm))_k
  // Fully expanded for D3Q19; s_q terms and j=2,16,17,18 terms vanish by symmetry.
  T Cxx = omega * (T{1} / T{3} * dm[1] + T{1} / T{3} * dm[9]);

  T Cyy = omega * (T{1} / T{3} * dm[1] - T{1} / T{6} * dm[9] + T{1} / T{2} * dm[11]);

  T Czz = omega * (T{1} / T{3} * dm[1] - T{1} / T{6} * dm[9] - T{1} / T{2} * dm[11]);

  T Cxy = omega * dm[13];

  T Cxz = omega * dm[14];

  T Cyz = omega * dm[15];

  // mu = nu * rho = (1/omega - 0.5) * cs² * rho
  T invOmega = T{1} / omega;
  T nu = (invOmega - T{0.5}) * LatSet::cs2;
  T mu = nu * rho;

  // F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
  T delta_rho = pf_cell.template get<DELTARHO<T>>();
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T prefactor = -T{3} * mu * delta_rho / rho;

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += prefactor * (Cxx * grad_phi[0] + Cxy * grad_phi[1] + Cxz * grad_phi[2]);
  ns_force[1] += prefactor * (Cxy * grad_phi[0] + Cyy * grad_phi[1] + Cyz * grad_phi[2]);
  ns_force[2] += prefactor * (Cxz * grad_phi[0] + Cyz * grad_phi[1] + Czz * grad_phi[2]);
}

}  // namespace ff
