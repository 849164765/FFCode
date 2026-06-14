// phase_c.ur.h
// Phase C: Hydrodynamic g-LBE (Eqs. 13-32)
// from Fakhari et al., Phys. Rev. E 96, 053301 (2017)
//
// C1: ρ = ρ_L + (φ-φ_L)*(ρ_H-ρ_L)                        Eq.13
// C2: μ = μ_L + (φ-φ_L)*(μ_H-μ_L)                        Eq.24
// C3: τ = μ/(ρ·c_s²)                                      Eq.25
// C5: μ_φ = 4β·(φ-φ_L)(φ-φ_H)(φ-φ_0) - κ·∇²φ             Eq.5
// C6: F_s = μ_φ·∇φ                                        Eq.4
// C7: F_p = -p*·c_s²·∇ρ                                    Eq.19
// C10: F_α = δt·w_α·(e_α·F)/(ρ·c_s²)                     Eq.15
// C12: g-LBE MRT collision                                Eq.14+27
// C13: p* = Σ g_α,  u = Σ g_α·e_α + F·δt/(2ρ)           Eq.32

#pragma once

#include "data_struct/cell.h"
#include "lbm/lattice_set.h"
#include "ff/ff2d.h"          // RHO_L, RHO_H, ETA_L, ETA_H, BETA, KAPPA, GRAVITY, CHEMICALPOTENTIAL, LAPLACIAN


namespace phase_field {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;
#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// ===================================================================
// C1: ComputeDensityFromPhase — ρ = ρ_L + φ·(ρ_H - ρ_L)    Eq.13
// Reads:  PHI<T>, RHO_L<T>, RHO_H<T>
// Writes: RHO<T>
// ===================================================================
template <typename CELLTYPE>
struct ComputeDensityFromPhase {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T phi = cell.template get<PHI<T>>();
    const T rho_l = cell.template get<ff::RHO_L<T>>();
    const T rho_h = cell.template get<ff::RHO_H<T>>();
    const T rho = rho_l + phi * (rho_h - rho_l);
    cell.template get<RHO<T>>() = rho;
  }
};


// ===================================================================
// C2+C3: ComputeViscosityOmega — μ = μ_L + φ·(μ_H - μ_L),  τ = μ/(ρ·c_s²)
// Reads:  PHI<T>, ETA_L<T>, ETA_H<T>, RHO<T>
// Writes: OMEGA<T>
// ===================================================================
template <typename CELLTYPE>
struct ComputeViscosityOmega {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T phi = cell.template get<PHI<T>>();
    const T eta_l = cell.template get<ff::ETA_L<T>>();
    const T eta_h = cell.template get<ff::ETA_H<T>>();
    const T rho = cell.template get<RHO<T>>();
    const T eta = eta_l + phi * (eta_h - eta_l);
    const T tau = eta / (rho * LatSet::cs2);
    const T omega = T{1} / (tau + T{0.5});
    cell.template get<OMEGA<T>>() = omega;
  }
};


// ===================================================================
// C5: ComputeChemicalPotential — μ_φ = 4β·φ(φ-1)(φ-0.5) - κ·∇²φ
// Reads:  PHI<T>, BETA<T>, KAPPA<T>, LAPLACIAN<T>
// Writes: CHEMICALPOTENTIAL<T>
// ===================================================================
template <typename CELLTYPE>
struct ComputeChemicalPotential {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T phi = cell.template get<PHI<T>>();
    const T beta = cell.template get<ff::BETA<T>>();
    const T kappa = cell.template get<ff::KAPPA<T>>();
    const T laplacian = cell.template get<ff::LAPLACIAN<T>>();
    // φ(φ-1)(φ-0.5) = φ³ - 1.5φ² + 0.5φ
    const T double_well = phi * (phi - T{1}) * (phi - T{0.5});
    cell.template get<ff::CHEMICALPOTENTIAL<T>>() = T{4} * beta * double_well - kappa * laplacian;
  }
};


// ===================================================================
// C4+C6+C7: ComputeForces — surface tension + pressure + body forces
// ∇ρ = (ρ_H-ρ_L)·∇φ  (Eq.33)
// F_s = μ_φ·∇φ        (Eq.4)
// F_p = -p*·c_s²·∇ρ   (Eq.19)
// F = F_s + F_p + F_b  (Eq.18, F_b from GRAVITY constant)
//
// Reads:  GRAD<T,D>, CHEMICALPOTENTIAL<T>, PHI<T>, RHO_L<T>, RHO_H<T>,
//         GRAVITY<T>, pressure p* (via POP sum or RHO<T> as proxy)
// Writes: FORCE<T,D> (accumulates, use +=)
// ===================================================================
template <typename CELLTYPE>
struct ComputeForces {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const Vector<T, LatSet::d>& grad = cell.template get<GRAD<T, LatSet::d>>();
    const T chem_pot = cell.template get<ff::CHEMICALPOTENTIAL<T>>();
    const T rho_l = cell.template get<ff::RHO_L<T>>();
    const T rho_h = cell.template get<ff::RHO_H<T>>();
    const T p_star = cell.template get<RHO<T>>();  // pressure proxy
    const T g = cell.template get<ff::GRAVITY<T>>();

    // C6: Surface tension F_s = μ_φ * ∇φ
    // C7: Pressure force F_p = -p* * c_s² * (ρ_H-ρ_L) * ∇φ   (Eq.33)
    // C9: Total force F = F_s + F_p + F_b
    auto& force = cell.template get<FORCE<T, LatSet::d>>();
    const T drho = rho_h - rho_l;
    for (unsigned int d = 0; d < LatSet::d; ++d) {
      force[d] = chem_pot * grad[d]
               - p_star * LatSet::cs2 * drho * grad[d];
      // Body force in negative-y (gravity) direction
      if (d == 1) force[d] += g;
    }
  }
};


// ===================================================================
// C0: UpdatePressureVelocityG — computes p* and u from g-populations
// This runs BEFORE G_Collision to provide macroscopic fields for
// force computation and equilibrium.
//
// p* = sum(g_alpha)                                    Eq.32a
// u  = sum(g_alpha * e_alpha) + F * dt / (2 * p*)     Eq.32b
//      (p* used as rho proxy in half-force correction; valid for Ma<<1)
//
// Reads:  POP<T,q>, FORCE<T,D>
// Writes: RHO<T> (= p*), VELOCITY<T,D>
// ===================================================================
template <typename CELLTYPE>
struct UpdatePressureVelocityG {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();

    // p* = sum(g)
    T p_star = T{0};
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha)
      p_star += cell[alpha];

    // momentum = sum(g*e)
    Vector<T, LatSet::d> u;
    for (unsigned int d = 0; d < LatSet::d; ++d) u[d] = T{0};
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const auto& e_a = latset::c<LatSet>(alpha);
      for (unsigned int d = 0; d < LatSet::d; ++d)
        u[d] += cell[alpha] * e_a[d];
    }

    // Half-force correction: u += F * dt / (2 * p*)
    // p* as rho proxy (incompressible approximation, valid for moderate density ratios)
    const T half_over_p = T{0.5} / p_star;
    for (unsigned int d = 0; d < LatSet::d; ++d)
      u[d] += F[d] * half_over_p;

    cell.template get<RHO<T>>() = p_star;
    cell.template get<VELOCITY<T, LatSet::d>>() = u;
  }
};

// ===================================================================
// C10-C12: G_Collision — BGK collision for g-LBE
// g_eq_α = p*·w_α + (Γ_α - w_α)                            Eq.17
// F_α = w_α·(e_α·F)/(ρ·c_s²)                              Eq.15
// ḡ_eq_α = g_eq_α - ½ F_α                                  Eq.16
// g_α_new = (1-ω)·g_α + ω·ḡ_eq_α + F_α                     Eq.14 BGK
// (MRT: g_α_new = g_α - M⁻¹ŜM(g_α - ḡ_eq_α) + F_α)
//
// Reads:  RHO<T>, VELOCITY<T,D>, FORCE<T,D>, OMEGA<T>
//         g_α populations (via cell[i])
// Writes: g_α populations (collided values)
// ===================================================================
template <typename CELLTYPE>
struct G_Collision {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const T p_star = cell.template get<RHO<T>>();
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& F = cell.template get<FORCE<T, LatSet::d>>();
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T rho = p_star;  // density = pressure proxy for incompressible

    // Compute F·u for force correction
    T Fu = T{0};
    for (unsigned int d = 0; d < LatSet::d; ++d) Fu += F[d] * u[d];

    // Precompute u²
    T u2 = T{0};
    for (unsigned int d = 0; d < LatSet::d; ++d) u2 += u[d] * u[d];

    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const auto& e_a = latset::c<LatSet>(alpha);
      const T w_a = latset::w<LatSet>(alpha);
      T c_u = T{0}, c_F = T{0};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        c_u += e_a[d] * u[d];
        c_F += e_a[d] * F[d];
      }

      // Eq.17: g_eq_α = p*·w_α + (Γ_α - w_α)
      // Γ_α = w_α·[1 + c·u/c_s² + (c·u)²/(2c_s⁴) - u²/(2c_s²)]
      const T Gamma = w_a * (T{1}
        + c_u * LatSet::InvCs2
        + c_u * c_u * LatSet::InvCs4 * T{0.5}
        - u2 * LatSet::InvCs2 * T{0.5});
      const T g_eq = p_star * w_a + (Gamma - w_a);

      // Eq.15: F_α = w_α·(e_α·F)/(ρ·c_s²)
      const T F_alpha = w_a * c_F / (rho * LatSet::cs2);

      // Eq.16: ḡ_eq = g_eq - ½ F_α
      const T gbar_eq = g_eq - T{0.5} * F_alpha;

      // Eq.14 BGK: g_new = (1-ω)·g_old + ω·ḡ_eq + F_α
      cell[alpha] = _omega * cell[alpha] + omega * gbar_eq + F_alpha;
    }
  }
};


// ===================================================================
// C13: UpdatePressureVelocity — p* = Σ g_α,  u = Σ g_α·e_α + F/(2ρ)
// Reads:  g_α populations (via cell[i]), FORCE<T,D>
// Writes: RHO<T> (= p*), VELOCITY<T,D> (= u)
// ===================================================================
template <typename CELLTYPE>
struct UpdatePressureVelocity {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    const Vector<T, LatSet::d>& F = cell.template get<FORCE<T, LatSet::d>>();
    const T rho = cell.template get<RHO<T>>();

    // p* = Σ g_α                             Eq.32a
    T p_star = T{0};
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha)
      p_star += cell[alpha];

    // u = Σ g_α·e_α + F·δt/(2ρ)              Eq.32b
    Vector<T, LatSet::d> u;
    for (unsigned int d = 0; d < LatSet::d; ++d) u[d] = T{0};
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const auto& e_a = latset::c<LatSet>(alpha);
      for (unsigned int d = 0; d < LatSet::d; ++d)
        u[d] += cell[alpha] * e_a[d];
    }
    const T half_dt_over_rho = T{0.5} / rho;
    for (unsigned int d = 0; d < LatSet::d; ++d)
      u[d] += F[d] * half_dt_over_rho;

    cell.template get<RHO<T>>() = p_star;
    cell.template get<VELOCITY<T, LatSet::d>>() = u;
  }
};

}  // namespace phase_field
