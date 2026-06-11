#pragma once

#include "lbm/moment.h"

// legacy namespace for backward compat
namespace legacy {

template <typename T, typename LatSet>
struct BGK {
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void apply(BasicPopCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(BasicPopCell<T, LatSet>& cell, const std::array<T, LatSet::q>& source) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * source[i];
    }
  }
};

}  // namespace legacy

// ===================================================================
//  Section 1: Phase-Field BGK Collision (Allen-Cahn)
//
//  f̃_α = f_α - ω_φ·(f_α - f_α^eq) + (1 - ω_φ/2)·F_α
//
//  f_α^eq = ω_α·φ·(1 + e_α·u / c_s²)
//  F_α    = ω_α·(e_α·n) · 4φ(1-φ)/W
// ===================================================================

template <typename CELL>
struct PhaseFieldBGK {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    // --- macroscopic: φ from field (pre-computed by PhaseFieldPhi) ---
    const T phi = cell.template get<typename CELL::GenericRho>();

    // --- get velocity and normal from field ---
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL<T, LatSet::d>>();
    const T W = cell.template get<INTERFACEWIDTH<T>>();

    // --- equilibrium: f_α^eq = ω_α·φ·(1 + e_α·u/c_s²) ---
    std::array<T, LatSet::q> feq;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T eu = u * latset::c<LatSet>(i);
      feq[i] = latset::w<LatSet>(i) * phi * (T{1} + LatSet::InvCs2 * eu);
    }

    // --- source: F_α = ω_α·(e_α·n) · 4φ(1-φ)/W ---
    const T coeff = T{4} * phi * (T{1} - phi) / W;
    std::array<T, LatSet::q> F;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = latset::c<LatSet>(i) * n;
      F[i] = latset::w<LatSet>(i) * en * coeff;
    }

    // --- BGK: f̃ = f - ω(f - feq) + (1-ω/2)F ---
    const T omega = cell.getOmega();
    const T halfOmega = T{0.5} * omega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = cell[i] - omega * (cell[i] - feq[i]) + (T{1} - halfOmega) * F[i];
    }
  }
};


// ===================================================================
//  Section 2: Velocity-Based BGK Collision (Navier-Stokes)
//
//  g̃_α = g_α - ω_g·(g_α - g_α^eq) + (1 - ω_g/2)·G_α
//
//  g_α^eq = ω_α·[ p/(ρc_s²) + e_α·u/c_s² + (e_α·u)²/(2c_s⁴) - u²/(2c_s²) ]
//
//  G_α = ω_α·{ e_α·F/c_s² + [(u⊗F+F⊗u) : (e_α⊗e_α - c_s²·I)] / (2c_s⁴) }
//
//  Macroscopic recovery:
//    p = ρ·c_s²·Σ g_α
//    u = Σ g_α·e_α + dt/(2ρ)·F^total
// ===================================================================

template <typename CELL, typename MomentScheme>
struct VelocityBGK {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    // --- macroscopic: p, u* from g_α ---
    T p;
    Vector<T, LatSet::d> u_star;
    MomentScheme::apply(cell, p, u_star);

    // --- full velocity: u = u* + dt/(2ρ)·F^total ---
    const T rho = cell.template get<typename CELL::GenericRho>();
    const Vector<T, LatSet::d>& F_total = cell.template get<FORCE<T, LatSet::d>>();
    const T coeff_u = T{1} / (T{2} * rho);  // dt = 1 in lattice units
    Vector<T, LatSet::d> u = u_star + F_total * coeff_u;

    // --- write velocity back to field ---
    cell.template get<VELOCITY<T, LatSet::d>>() = u;

    // --- equilibrium: g_α^eq = ω_α·[p/(ρc_s²) + e·u/c_s² + (e·u)²/(2c_s⁴) - u²/(2c_s²)] ---
    const T u2 = u.getnorm2();
    const T p_over_rho_cs2 = p / (rho * LatSet::cs2);

    std::array<T, LatSet::q> geq;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T eu = latset::c<LatSet>(i) * u;
      geq[i] = latset::w<LatSet>(i) * (p_over_rho_cs2 + LatSet::InvCs2 * eu
               + T{0.5} * eu * eu * LatSet::InvCs4 - T{0.5} * LatSet::InvCs2 * u2);
    }

    // --- force term: G_α = ω_α·{e·F/c_s² + [(uF+Fu):(ee-c_s²I)]/(2c_s⁴)} ---
    std::array<T, LatSet::q> G;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const auto& e = latset::c<LatSet>(i);

      T term1 = (e * F_total) * LatSet::InvCs2;

      T term2 = T{0};
      for (unsigned int a = 0; a < LatSet::d; ++a) {
        for (unsigned int b = 0; b < LatSet::d; ++b) {
          T eab = e[a] * e[b];
          T delta = (a == b) ? LatSet::cs2 : T{0};
          term2 += (u[a] * F_total[b] + F_total[a] * u[b]) * (eab - delta);
        }
      }
      term2 *= T{0.5} * LatSet::InvCs4;

      G[i] = latset::w<LatSet>(i) * (term1 + term2);
    }

    // --- BGK: g̃ = g - ω(g - geq) + (1-ω/2)G ---
    const T omega = cell.getOmega();
    const T halfOmega = T{0.5} * omega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = cell[i] - omega * (cell[i] - geq[i]) + (T{1} - halfOmega) * G[i];
    }
  }
};
