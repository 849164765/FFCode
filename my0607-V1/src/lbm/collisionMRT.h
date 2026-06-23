// collisionMRT.h
// MRT (multiple-relaxation-time) collision operators
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//   Eq.15: PhaseFieldMRT  (D2Q9 MRT for phase-field)
//   Eq.29: NSMRT           (D2Q9 MRT for velocity-based NS)
//   Eq.38: MagMRT          (D2Q5 MRT for magnetic solver)

#pragma once

#include <array>

#include "lbm/mrt_matrix.h"
#include "lbm/equilibrium.h"
#include "lbm/force.h"
#include "fflbm/phasefield_source.h"
#include "utils/alias.h"

// ================================================================
// PhaseFieldMRT (Eq.15) — D2Q9 MRT for phase-field
// f_tilde = f - M^{-1} S M (f - f_eq) + M^{-1}(I - S/2) M F  (F = F_beta from Eq.19)
// Relaxation: S = diag(1,1,1,s3,1,s5,1,1,1)
//   1/s3 = 1/s5 = M_phi/(cs^2*dt) + 0.5  (Eq.22)
// ================================================================
template <typename CELLTYPE>
struct PhaseFieldMRT {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T f[9];
    std::array<T, 9> feq{};
    T m[9], m_eq[9];

    // 1. Read populations
    for (int k = 0; k < 9; ++k) f[k] = cell[k];

    // 2. Equilibrium
    PhaseFieldEquilibrium<CELL>::apply(cell, feq);

    // 3. MRT forward transform: m = M9 * f, m_eq = M9 * feq
    mrt::M9Transform(m, f);
    mrt::M9Transform(m_eq, feq.data());

    // 4. Relaxation: S = [1,1,1,s3,1,s5,1,1,1]
    //    s3 = s5 = 1/tau_phi,  tau_phi = 0.5 + M_phi/cs^2  (Eq.22)
    T tau_pf = cell.template get<TAUPHI<T>>();
    T s_pf = T{1} / tau_pf;
    T S[9] = {T{1}, T{1}, T{1}, s_pf, T{1}, s_pf, T{1}, T{1}, T{1}};

    // 5. Source term in moment space: M_F = M9 * F_beta
    //    m += (I - S/2) * M_F   (Eq.15, source from Eq.19)
    std::array<T, 9> F_src{};
    PhaseFieldSource<CELL>::GetSource(cell, F_src);
    T M_F[9];
    mrt::M9Transform(M_F, F_src.data());
    for (int k = 0; k < 9; ++k) {
      m[k] += (T{1} - S[k] * T{0.5}) * M_F[k];
    }

    // 6. Collision in moment space: m -= S * (m - m_eq)
    for (int k = 0; k < 9; ++k) {
      m[k] -= S[k] * (m[k] - m_eq[k]);
    }

    // 7. MRT inverse transform: f = invM9 * m
    mrt::M9InverseTransform(f, m);

    // 8. Write back populations
    for (int k = 0; k < 9; ++k) cell[k] = f[k];
  }
};

// ================================================================
// NSMRT (Eq.29) — D2Q9 MRT for velocity-based Navier-Stokes
// g_tilde = g - M^{-1} S M (g - g_eq) + M^{-1}(I - S/2) M G
// Relaxation: S = diag(1,1,1,1,1,1,1,s7,s8)
//   1/s7 = 1/s8 = eta/(rho*cs^2*dt) + 0.5  (Eq.33)
// After collision: sets VELOCITY (Eq.35) and PRESSURE (Eq.36)
// ================================================================
// NSMRT (Eq.29): D2Q9 MRT collision for velocity-based NS
// Reads VELOCITY and PRESSURE from cell field (set by moment::NSMomentum
// as a SEPARATE task BEFORE this collision). Does NOT write VELOCITY/PRESSURE.
// ================================================================
template <typename CELLTYPE>
struct NSMRT {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T g[9];
    std::array<T, 9> geq{};
    T m[9], m_eq[9];

    // 1. Read populations and physical fields
    for (int k = 0; k < 9; ++k) g[k] = cell[k];
    T rho = cell.template get<RHO<T>>();
    T eta = cell.template get<VISCOSITY<T>>();

    // 2. Relaxation: s7 = s8 = 1/(eta/(rho*cs^2) + 0.5)  (Eq.33)
    T tau_nu = eta / (rho * LatSet::cs2) + T{0.5};
    T s_nu = T{1} / tau_nu;
    T S[9] = {T{1}, T{1}, T{1}, T{1}, T{1}, T{1}, T{1}, s_nu, s_nu};

    // 3. Equilibrium (reads VELOCITY/PRESSURE from field set by NSMomentum)
    NSEquilibrium<CELL>::apply(cell, geq);

    // 4. MRT forward transform
    mrt::M9Transform(m, g);
    mrt::M9Transform(m_eq, geq.data());

    // 5. Force in moment space: M_G = M9 * G, then m += (I - S/2) * M_G
    std::array<T, 9> G{};
    GuoForce<CELL>::GetForce(cell, G);
    T M_G[9];
    mrt::M9Transform(M_G, G.data());
    for (int k = 0; k < 9; ++k) {
      m[k] += (T{1} - S[k] * T{0.5}) * M_G[k];
    }

    // 6. Collision in moment space: m -= S * (m - m_eq)
    for (int k = 0; k < 9; ++k) {
      m[k] -= S[k] * (m[k] - m_eq[k]);
    }

    // 7. MRT inverse transform
    mrt::M9InverseTransform(g, m);

    // 8. Write back populations
    for (int k = 0; k < 9; ++k) cell[k] = g[k];
  }
};

// ================================================================
// MagMRT (Eq.38) — D2Q5 MRT for magnetic solver
// h_tilde = h - M_5^{-1} S M_5 (h - h_eq)
// Note: NO force source term (Eq.38 has no F term)
// Relaxation: S = diag(1, s1, s2, 1, 1)
//   1/s1 = 1/s2 = 2.5*epsilon*mu/(c^2*dt) + 0.5  (Eq.42)
// ================================================================
template <typename CELLTYPE>
struct MagMRT {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T h[9];
    std::array<T, 9> heq{};
    T m[9], m_eq[9];

    for (int k = 0; k < 9; ++k) h[k] = cell[k];

    // 2. Equilibrium: h_alpha^eq = w_alpha * psi
    MagEquilibrium<CELL>::apply(cell, heq);

    // 3. MRT forward transform: m = M5 * h
    mrt::M9Transform(m, h);
    mrt::M9Transform(m_eq, heq.data());

    constexpr T s_h = T{1};
    T S[9] = {T{1}, T{1}, T{1}, s_h, T{1}, s_h, T{1}, T{1}, T{1}};

    for (int k = 0; k < 9; ++k)
      m[k] -= S[k] * (m[k] - m_eq[k]);

    mrt::M9InverseTransform(h, m);

    for (int k = 0; k < 9; ++k) cell[k] = h[k];
  }
};
