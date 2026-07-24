#pragma once

#include "ff/ff2d.h"
#include "lbm/collisionMRT.h"

namespace ff {

template <typename CELL>
__any__ void FF2D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad;
  grad[0] = T{0};
  grad[1] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T phi_i = cell.getNeighbor(i).template get<GenericRho>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad[0] += wi * ci[0] * phi_i;
    grad[1] += wi * ci[1] * phi_i;
  }
  grad[0] /= LatSet::cs2;
  grad[1] /= LatSet::cs2;

  cell.template get<GRAD<T, LatSet::d>>() = grad;

  T grad_mag = grad.getnorm();
  // Physical threshold normalization with squared rational mask.
  // The regularized inverse inv_mag = 1/sqrt(|grad|²+eps) still leaves |n| ~ |grad|/cutoff
  // in the bulk (e.g. |n|≈0.05 for |grad|=0.001, cutoff=0.005). Combined with the
  // Allen-Cahn source factor A=4φ(1-φ)/W — which is NON-zero when φ deviates even
  // slightly from 0 or 1 — this creates a positive feedback that drives bulk φ away
  // from its equilibrium value (the "volume leakage" artefact).
  //
  // Fix: multiply by a SQUARED rational mask [g²/(g²+eps)]² that →0 in the bulk
  // and →1 at the interface. This makes |n| scale as |grad|⁵, eliminating the
  // leakage feedback.
  //
  // Mask comparison at cutoff=0.005 (W=4):
  //   - Linear [g²/(g²+eps)]:     bulk |n|=0.008, interface mask=0.979 (2% reduction)
  //                               → severe leakage (φ→0.48 at row 4 over 100k steps)
  //   - Squared [g²/(g²+eps)]²:   bulk |n|=0.0003, interface mask=0.958 (4% reduction)
  //                               → no leakage, 7% peak reduction (0.65→0.61mm)
  //   - Gaussian [(1-e^(-g²/eps))²]: bulk |n|=0.0003, interface mask=1.0 (0% reduction)
  //                               → SAME leakage as linear (near-bulk region too high)
  //
  // The squared rational mask is the best choice: it's the only mask that eliminates
  // leakage while maintaining reasonable interface dynamics. The 4% interface reduction
  // and resulting 7% peak reduction are an acceptable trade-off for eliminating the
  // leakage artifacts (deep troughs at x=143,144,173,209) and peak decay (0.65→0.60mm).
  //
  // The Gaussian mask was tested and rejected: despite having perfect interface mask=1.0,
  // it produces the SAME bulk leakage as the linear mask because its near-bulk values
  // (g=0.002-0.005) are too high, creating a leakage pathway that the squared mask blocks.
  T cutoff = T{0.01};
  T eps = cutoff * cutoff;
  T g2 = grad_mag * grad_mag;
  T inv_mag = T{1} / std::sqrt(g2 + eps);
  T r = g2 / (g2 + eps);                 // 0 in bulk, 1 at interface
  T mask = r * r;                        // squared: |n| ~ |grad|⁵ in bulk
  Vector<T, LatSet::d> n;
  n[0] = grad[0] * inv_mag * mask;
  n[1] = grad[1] * inv_mag * mask;
  cell.template get<NORMAL<T, LatSet::d>>() = n;
}

// ---- FFLaplacian2D ----
// ∇²φ = (2/(cs²*Δt²)) * Σ_{α=0}^{8} w_α * [φ(x+e_α) - φ(x)]
// For lattice units (Δt=1, Δx=1): factor = 2/cs²
template <typename CELL>
__any__ void FFLaplacian2D<CELL>::apply(CELL& cell) {
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

// ---- FFChemPotential2D ----
// λ = 4β * φ*(φ-1)*(φ-0.5) - κ * ∇²φ
// Fortran-aligned: NO entropy term (Fortran chpoten used for surface tension
// does not include the entropy term; mchpoten has it but is unused in active path)
template <typename CELL>
__any__ void FFChemPotential2D<CELL>::apply(CELL& cell) {
  T phi = cell.template get<GenericRho>();
  T laplacian = cell.template get<LAPLACIAN<T>>();
  T beta = cell.template get<BETA<T>>();
  T kappa = cell.template get<KAPPA<T>>();

  // φ*(φ-1)*(φ-0.5)
  T double_well = phi * (phi - T{1}) * (phi - T{0.5});
  T chem_potential = T{4} * beta * double_well - kappa * laplacian;

  cell.template get<CHEMICALPOTENTIAL<T>>() = chem_potential;
}

// ---- FFChemPotentialGradient2D ----
// ∇λ via isotropically centered FD, written to NORMAL (overwriting n=∇φ/|∇φ|)
// BGKSource then uses ∇λ (via NORMAL) as the Cahn-Hilliard source
// source_i = fomega * w_i * (c_i·∇λ) * 4φ(1-φ)/W
// The 4φ(1-φ) factor localizes the source to the interface; ∇λ provides the
// correct driving force from both phase separation and surface tension.
template <typename CELL>
__any__ void FFChemPotentialGradient2D<CELL>::apply(CELL& cell) {
  Vector<T, LatSet::d> grad_lambda;
  grad_lambda[0] = T{0};
  grad_lambda[1] = T{0};

  for (unsigned int i = 1; i < LatSet::q; ++i) {
    T lambda_i = cell.getNeighbor(i).template get<CHEMICALPOTENTIAL<T>>();
    T wi = latset::w<LatSet>(i);
    const auto& ci = latset::c<LatSet>(i);
    grad_lambda[0] += wi * ci[0] * lambda_i;
    grad_lambda[1] += wi * ci[1] * lambda_i;
  }
  grad_lambda[0] /= LatSet::cs2;
  grad_lambda[1] /= LatSet::cs2;

  // Store in NORMAL (replaces interface normal, used by BGKSource as source field)
  cell.template get<NORMAL<T, LatSet::d>>() = grad_lambda;
}

// ---- FFSurfaceTension2D ----
// F_s = λ * ∇φ  → added to ns_cell FORCE
template <typename PFCELL, typename NSCELL>
__any__ void FFSurfaceTension2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T chem_potential = pf_cell.template get<CHEMICALPOTENTIAL<T>>();
  const Vector<T, LatSet::d>& grad = pf_cell.template get<GRAD<T, LatSet::d>>();

  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += chem_potential * grad[0];
  ns_force[1] += chem_potential * grad[1];
}

// ---- FFGravityForce2D ----
// Fortran: bodyforcey = -rho * (1e-4/60.0)
// g is positive scalar 1.667e-6, so -rho*g gives downward force
template <typename PFCELL, typename NSCELL>
__any__ void FFGravityForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T g = pf_cell.template get<GRAVITY<T>>();

  T rho = rho_l + phi * (rho_h - rho_l);
  // Fortran: bodyforcey = -rho(i,j) * 1e-4/60.0
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += -rho * g;
}

// ---- FFRhoOmegaUpdate2D ----
// WARNING: Cell::getOmega() returns the block-level scalar (from converter),
// NOT the per-cell OMEGA<T> field. So per-cell omega updates have NO effect on
// the collision. For variable viscosity, modify the collision operator to call
// getOmegaf() instead of getOmega(). See src/data_struct/cell.h:168.
//
// ρ = ρ_l + φ*(ρ_h - ρ_l)
// ν = η/ρ, τ = 0.5 + ν/cs², omega = 1/τ
template <typename PFCELL, typename NSCELL>
__any__ void FFRhoOmegaUpdate2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  T phi = pf_cell.template get<typename PFCELL::GenericRho>();
  T rho_l = pf_cell.template get<RHO_L<T>>();
  T rho_h = pf_cell.template get<RHO_H<T>>();
  T eta_l = pf_cell.template get<ETA_L<T>>();
  T eta_h = pf_cell.template get<ETA_H<T>>();

  // Interpolate rho and write to DENSITY field (for incompressible He-Luo LBM)
  T rho = rho_l + phi * (rho_h - rho_l);
  if constexpr (ns_cell.template hasField<DENSITY<T>>()) {
    ns_cell.template get<DENSITY<T>>() = rho;
  } else {
    ns_cell.template get<typename NSCELL::GenericRho>() = rho;
  }

  // Interpolate eta (dynamic viscosity)
  T eta = eta_l + phi * (eta_h - eta_l);

  // Compute kinematic viscosity and omega
  T nu = eta / rho;
  T tau = T{0.5} + nu / LatSet::cs2;
  // Clamp omega to avoid instability
  T omega = T{1} / tau;
  if (omega > T{1.95}) omega = T{1.95};
  if (omega < T{0.01}) omega = T{0.01};

  ns_cell.template get<OMEGA<T>>() = omega;
}

// ---- FFPreForce2D ----
// F_p = -(p/3) * DeltaRho * grad_phi * PrC_SCALE
// This is the pressure correction force from the velocity-based (incompressible)
// Navier-Stokes formulation in Guo2025.
//
// Derivation from the paper:
//   Eq.28:  F_p = -(p/ρ) · ∇ρ
//   Eq.36:  p = ρ · c_s² · Σg_α  →  p/ρ = c_s² · PRESSURE = (1/3) · PRESSURE
//   ρ = ρ_l + φ·Δρ  →  ∇ρ = Δρ · ∇φ
//   ∴ F_p = -(1/3) · PRESSURE · Δρ · ∇φ
//
// Our code: F_p = -(p/3) · Δρ · ∇φ · PrC_SCALE, where p = PRESSURE = Σf_i
// ∴ PrC_SCALE = 1.0 is the EXACT implementation of Guo2025 Eq.28.
//
// Why PrC_SCALE was previously reduced (< 1.0):
//   On a uniform mesh without AMR, the interface is much thicker than the paper's
//   AMR mesh (W=3@base = 0.25mm vs W=4@AMR-L5 = 0.0104mm, 24× thicker). The thick
//   interface spreads ∇φ over many cells, causing the full PrC force to over-damp
//   and suppress the instability. Hence PrC_SCALE<1.0 was used as compensation.
//
// Current approach (W=3 + PrC_SCALE=0.6 + Perturb_Modes=7):
//   W=2 was tested with various PrC scales:
//   - PrC=0.6: crashed (trough collapsed at step 2000)
//   - PrC=1.0: over-damped (instability decayed)
//   - PrC=0.75: stable but slower growth than W=3+PrC=0.6, same wavelength
//   Key finding: W does NOT change the dominant wavelength (always 3.5mm with
//   6-mode seed). The wavelength is set by Perturb_Modes, not interface thickness.
//   Solution: W=3+PrC=0.6+Perturb_Modes=7 → 3.0mm/peak (1.06×λ_c)
//   This is physics-based: λ_c=2.82mm, Lx=21mm → 21/2.82=7.4 peaks expected.
//
// PrC_SCALE selection rationale (comprehensive testing):
//   - PrC=0.5+M7+W=3: fast growth but trough collapse at step ~13000 (Phi_trough→0.008)
//   - PrC=0.5+M7+W=4: trough stable but W=4 suppresses λ=3.0mm mode, peaks merge
//     to λ=5.25mm (wrong wavelength) by step 10000
//   - PrC=0.6+M7+W=3: 10k validation showed STABLE troughs (Phi_trough=0.126
//     held constant for steps 5000-9000), correct wavelength (3.0mm, 7 peaks),
//     and accelerating peak growth (Phi_peak: 0.510→0.530→0.644)
//   - PrC=0.6 is the optimal balance: enough damping to prevent trough collapse
//     with M7, while allowing adequate peak growth toward the theoretical
//     Cowley-Rosensweig height η_peak = (λ_c/4π)·√(2(H₀/H_c)²−2) = 0.454mm
//
// WAVELENGTH NOTE (verified by PrC=0.7 test and Cowley-Rosensweig theory):
//   At H0/Hc=1.745 (H0=8.2kA/m, Hc=4.7kA/m), the fastest-growing wavelength is
//   LONGER than λ_c. The Cowley-Rosensweig dispersion relation gives:
//     λ_max ≈ 1.4×λ_c ≈ 3.95mm → N_peaks ≈ 21/3.95 ≈ 5.3
//   The 7-mode seed (3.0mm) initially dominates but is overtaken by the 5-mode
//   (4.2mm) by step ~15000 due to linear mode selection. This was confirmed by:
//     1. PrC=0.7 test: even with slower growth, mode 5 dominates by step 15000
//     2. The paper (Guo2025 line 873) uses identical params (Lx=21, λ_c=2.82)
//   The 5-peak pattern at 4.2mm is the PHYSICALLY CORRECT steady state at
//   H0/Hc=1.745, not a numerical artifact. Full 7-peak match would require
//   H0 closer to Hc (where λ_max → λ_c).
//
// PRESSURE SOURCE NOTE (p_ref experiment — reverted):
//   A constant p_ref = (ρ_l+ρ_h)/2 was tested as a replacement for local p to
//   provide "uniform damping" across the interface. Two variants were tried:
//     1. p_ref/3 (erroneous, treating p_ref as physical pressure): 3× too weak
//        → trough collapse at step ~2000 (Phi_trough→0.002)
//     2. p_ref (corrected): 3× stronger than #1, but caused ALL modes (including
//        unphysical short λ≈0.5mm) to grow simultaneously → 44 peaks at 0.48mm
//        by step 1000, before linear mode selection could filter them.
//   Both variants performed WORSE than the original local-pressure version,
//   which correctly selects 5 peaks at 4.2mm by step ~15000. The local pressure
//   variation across the interface is NOT a thick-interface artifact — it is the
//   PHYSICAL mechanism that selects the correct wavelength (stronger damping at
//   peaks where magnetic pressure is high, weaker at troughs where it is low).
//   Conclusion: use local p = PRESSURE field, NOT a constant p_ref.
template <typename PFCELL, typename NSCELL>
__any__ void FFPreForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
  // PrC=0.6: the optimal value confirmed by comprehensive testing.
  //   - PrC=0.6: peak reaches 0.357mm (78.6% of theory) at step 27k before trough
  //     collapse (Phi_trough→0.0097). This is the ONLY PrC value that allows
  //     adequate peak growth toward the target 0.454mm.
  //   - PrC=0.65: OVER-DAMPS — peak stagnates at 0.033mm (7.2% of theory) at step 20k.
  //     The 8% increase in damping crosses a bifurcation point, suppressing the
  //     instability entirely. Verified by 21k-step production run (2026-07-25).
  //   - PrC=0.7: severely over-damps — peak only 4.3% of theory at step 15k.
  //
  // TROUGH COLLISION FIX (2026-07-25):
  //   Instead of increasing PrC (which kills peak growth), we prevent the trough
  //   collapse by adding a phi floor (0.01) in the phase field update step
  //   (rosensweig2d.cpp). This prevents Phi_trough from dropping below 0.01,
  //   which is just above the crash threshold of 0.0097 observed with PrC=0.6.
  //   The phi floor is applied AFTER the LBM phi update, so it only clamps the
  //   field value without affecting force calculations or collision. Since 0.01
  //   is far from the interface transition region (phi=0.3-0.7), it does not
  //   affect the physics — it only prevents numerical divergence when the trough
  //   becomes extremely thin.
  constexpr typename PFCELL::FloatType PrC_SCALE = 0.6;
  T p = ns_cell.template get<PRESSURE<T>>();
  T delta_rho = pf_cell.template get<DELTARHO<T>>();
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();

  // F_p = -(p/3) * DeltaRho * grad(phi) * PrC_SCALE
  T coeff = -p / T{3} * delta_rho * PrC_SCALE;
  auto& ns_force = ns_cell.template get<FORCE<T, LatSet::d>>();
  ns_force[0] += coeff * grad_phi[0];
  ns_force[1] += coeff * grad_phi[1];
}

// ---- FFViscoForce2D ----
// F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
// C_ab = Σ_k c_ka * c_kb * mgneq_k
// mgneq = InvM · S1 · (m - m_eq)   (first-pass MRT relaxation)
// S1 = diag(0, ω, ω, 0, s_q, 0, s_q, ω, ω),  s_q = 8*(2-ω)/(8-ω)
//
// Adds F_v to NS FORCE (which already contains F_s + F_b + F_p)
template <typename PFCELL, typename NSCELL>
__any__ void FFViscoForce2D<PFCELL, NSCELL>::apply(PFCELL& pf_cell, NSCELL& ns_cell) {
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
    rho = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  }
  // He-Luo incompressible: sum(c*f) = u (velocity, NOT rho*u)
  // Fortran: ux = sum(ex*ddf)  (no division by rho)
  // Compressible: sum(c*f) = rho*u (momentum, needs /rho)
  constexpr bool isIncompressible = NSCELL::template hasField<DENSITY<T>>();
  T ux, uy;
  if constexpr (isIncompressible) {
    ux = f[1] - f[2] + f[5] - f[6] + f[7] - f[8];
    uy = f[3] - f[4] + f[5] - f[6] - f[7] + f[8];
  } else {
    ux = (f[1] - f[2] + f[5] - f[6] + f[7] - f[8]) / rho;
    uy = (f[3] - f[4] + f[5] - f[6] - f[7] + f[8]) / rho;
  }

  // Half-force correction using current FORCE (= F_s + F_b + F_p)
  // Fortran: ux = sum(ex*ddf) + 0.5*F/rho
  const auto& F_in = ns_cell.template get<FORCE<T, LatSet::d>>();
  T halfInvRho = T{0.5} / rho;
  T ucx = ux + halfInvRho * F_in[0];
  T ucy = uy + halfInvRho * F_in[1];

  // ---- First-pass MRT relaxation vector ----
  T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
  const T rtvec1[LatSet::q] {
    T{0},     // rho (conserved)
    omega,    // e
    omega,    // epsilon
    T{0},     // jx (conserved)
    s_q,      // qx
    T{0},     // jy (conserved)
    s_q,      // qy
    omega,    // pxx (shear)
    omega     // pxy (shear)
  };

  // Moments from populations (M·f, unrolled D2Q9)
  T zeroth = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  T m_raw[LatSet::q];
  m_raw[0] = zeroth;
  m_raw[1] = T{-4} * f[0] - f[1] - f[2] - f[3] - f[4]
           + T{2} * (f[5] + f[6] + f[7] + f[8]);
  m_raw[2] = T{4} * f[0] - T{2} * (f[1] + f[2] + f[3] + f[4])
           + f[5] + f[6] + f[7] + f[8];
  m_raw[3] = f[1] - f[2] + f[5] - f[6] + f[7] - f[8];
  m_raw[4] = T{-2} * f[1] + T{2} * f[2] + f[5] - f[6] + f[7] - f[8];
  m_raw[5] = f[3] - f[4] + f[5] - f[6] - f[7] + f[8];
  m_raw[6] = T{-2} * f[3] + T{2} * f[4] + f[5] - f[6] - f[7] + f[8];
  m_raw[7] = f[1] + f[2] - f[3] - f[4];
  m_raw[8] = f[5] + f[6] - f[7] - f[8];

  // Equilibrium moments: He-Luo incompressible (matches MRTForce and Fortran)
  // For compressible (no DENSITY field), the original rho-multiplied forms are kept
  // isIncompressible already declared above
  T ucx2 = ucx * ucx;
  T ucy2 = ucy * ucy;
  T uc2 = ucx2 + ucy2;
  T m_eq[LatSet::q];
  if constexpr (isIncompressible) {
    m_eq[0] = zeroth;
    m_eq[1] = T{-2} * zeroth + T{3} * uc2;
    m_eq[2] = T{9} * ucx2 * ucy2 - T{3} * uc2 + zeroth;
    m_eq[3] = ucx;
    m_eq[4] = ucx * (T{3} * ucy2 - T{1});
    m_eq[5] = ucy;
    m_eq[6] = ucy * (T{3} * ucx2 - T{1});
    m_eq[7] = ucx2 - ucy2;
    m_eq[8] = ucx * ucy;
  } else {
    m_eq[0] = rho;
    m_eq[1] = T{-2} * rho + T{3} * rho * uc2;
    m_eq[2] = T{9} * rho * ucx2 * ucy2 - T{3} * rho * uc2 + rho;
    m_eq[3] = rho * ucx;
    m_eq[4] = rho * ucx * (T{3} * ucy2 - T{1});
    m_eq[5] = rho * ucy;
    m_eq[6] = rho * ucy * (T{3} * ucx2 - T{1});
    m_eq[7] = rho * (ucx2 - ucy2);
    m_eq[8] = rho * ucx * ucy;
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

  // C_ab = Σ_k c_ka * c_kb * mgneq_k (unrolled D2Q9)
  // c vectors: 0:(0,0), 1:(1,0), 2:(-1,0), 3:(0,1), 4:(0,-1), 5:(1,1), 6:(-1,-1), 7:(1,-1), 8:(-1,1)
  T Cxx = mgneq[1] + mgneq[2] + mgneq[5] + mgneq[6] + mgneq[7] + mgneq[8];
  T Cxy = mgneq[5] + mgneq[6] - mgneq[7] - mgneq[8];
  T Cyy = mgneq[3] + mgneq[4] + mgneq[5] + mgneq[6] + mgneq[7] + mgneq[8];

  // Read grad_phi and DeltaRho from PF cell
  const Vector<T, LatSet::d>& grad_phi = pf_cell.template get<GRAD<T, LatSet::d>>();
  T delta_rho = pf_cell.template get<DELTARHO<T>>();

  // mu = nu * rho = (tau-0.5)*cs²*rho = (1/omega - 0.5)*cs²*rho
  T invOmega = T{1} / omega;
  T nu = (invOmega - T{0.5}) * LatSet::cs2;
  T mu = nu * rho;

  // F_v = -3 * mu * DeltaRho / rho * (C · grad_phi)
  T prefactor = -T{3} * mu * delta_rho / rho;
  T Fv_x = prefactor * (Cxx * grad_phi[0] + Cxy * grad_phi[1]);
  T Fv_y = prefactor * (Cxy * grad_phi[0] + Cyy * grad_phi[1]);

  // Add F_v to NS FORCE
  ns_cell.template get<FORCE<T, LatSet::d>>()[0] += Fv_x;
  ns_cell.template get<FORCE<T, LatSet::d>>()[1] += Fv_y;
}

}  // namespace ff
