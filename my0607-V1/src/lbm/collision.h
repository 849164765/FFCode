#pragma once

#include "lbm/moment.h"
#include "ff/ff2d.h"

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

template <typename CELL, typename EquilibriumScheme>
struct PhaseFieldBGK {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    // --- macroscopic: φ from field ---
    const T phi = cell.template get<typename CELL::GenericRho>();

    // --- get velocity from field ---
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();

    // --- equilibrium ---
    std::array<T, LatSet::q> feq;
    EquilibriumScheme::apply(feq, phi, u);

    // --- source ---
    std::array<T, LatSet::q> F;
    ff::PhaseFieldSource<CELL>::apply(cell, F);

    // --- BGK ---
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

template <typename CELL, typename EquilibriumScheme>
struct VelocityBGK {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // --- macroscopic from fields (computed by VelocityMomenta task) ---
    T p = cell.template get<PRESSURE<T>>();
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    T rho = cell.template get<GenericRho>();

    // --- equilibrium ---
    std::array<T, LatSet::q> geq;
    EquilibriumScheme::apply(geq, rho, p, u);

    // --- force term ---
    std::array<T, LatSet::q> G;
    ff::NSFieldSource<CELL>::apply(cell, G);

    // --- BGK ---
    const T omega = cell.getOmega();
    const T halfOmega = T{0.5} * omega;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = cell[i] - omega * (cell[i] - geq[i]) + (T{1} - halfOmega) * G[i];
    }
  }
};
