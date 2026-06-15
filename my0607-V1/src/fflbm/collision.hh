// collision.hh — BGK collision operators
//
// BGKCollision<MomentaScheme, EquilibriumScheme>
//   Plain BGK:  f_i = ω·feq_i + (1-ω)·f_i
//
// BGKForceCollision<MomentaScheme, EquilibriumScheme>
//   BGK + Guo force:  f_i = ω·feq_i + (1-ω)·f_i + (1-ω/2)·F_i
//   Reads FORCE field from cell.
//
// BGKSourceCollision<EquilibriumScheme>
//   BGK + Allen-Cahn source:
//     source_i = (1-ω/2)·w_i·(c_i·n)·4φ(1-φ)/W
//     f_i = ω·feq_i + (1-ω)·f_i + source_i
//   Reads NORMAL field (n vector) and INTERFACEWIDTH from cell.

#pragma once

namespace fflbm {

// ---------- plain BGK ----------
template <typename MomentaScheme, typename EquilibriumScheme>
struct BGKCollision {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);

    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);

    const T omega = cell.getOmega();
    const T omega_bar = cell.get_Omega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + omega_bar * cell[i];
    }
  }
};

// ---------- BGK + Guo force ----------
template <typename MomentaScheme, typename EquilibriumScheme>
struct BGKForceCollision {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    // read force vector from cell field
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();

    // force-corrected macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);

    // Guo discrete force term:
    // F_i = w_i * F · [(c_i - u)/cs² + (c_i·u)/cs⁴ * c_i]
    std::array<T, LatSet::q> Fi{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const auto& ci = latset::c<LatSet>(i);
      const T ci_dot_u = ci * u;
      Fi[i] = latset::w<LatSet>(i) * F *
              ((ci - u) * LatSet::InvCs2 +
               ci_dot_u * LatSet::InvCs4 * ci);
    }

    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);

    const T omega = cell.getOmega();
    const T omega_bar = cell.get_Omega();
    const T fomega = cell.getfOmega();  // = 1 - omega/2

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + omega_bar * cell[i] + fomega * Fi[i];
    }
  }
};

// ---------- BGK + Allen-Cahn source ----------
// source_i = (1-ω/2)·w_i·(c_i·n)·4φ(1-φ)/W
// n = ∇φ/|∇φ| (interface normal) from NORMAL field (computed by PhiGradient)
template <typename EquilibriumScheme>
struct BGKSourceCollision {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    const T phi = cell.template get<typename CELL::GenericRho>();
    const auto& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const auto& n = cell.template get<NORMAL<T, LatSet::d>>();
    const T W = cell.template get<INTERFACEWIDTH<T>>();

    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, phi, u);

    const T omega = cell.getOmega();
    const T omega_bar = cell.get_Omega();
    const T fomega = cell.getfOmega();  // = 1 - omega/2

    // 4φ(1-φ)/W — interface localization factor
    const T phi_factor = T{4} * phi * (T{1} - phi) / W;

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const T en = latset::c<LatSet>(i) * n;
      const T source = fomega * latset::w<LatSet>(i) * en * phi_factor;
      cell[i] = omega * feq[i] + omega_bar * cell[i] + source;
    }
  }
};

// ---------- Pressure-based BGK + Guo force (per-cell omega) ----------
// g-LBE for incompressible NS:
//   g_i^* = ω·g_i^eq + (1-ω)·g_i + (1-ω/2)·F_i
//   p̃ = Σg_i,  u = Σc_i·g_i + F/(2ρ₀)
//   ω = per-cell omega from OMEGA field (cell.getOmegaf())
//   MomentaScheme: computes p̃ and u (e.g. PressureMoment<CELL>)
//   EquilibriumScheme: computes g_i^eq (e.g. PressureEq<CELL>)
template <typename MomentaScheme, typename EquilibriumScheme>
struct PressureBGKCollision {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // read force vector
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();

    // force-corrected moments: p̃, u
    T p_cs2{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, p_cs2, u);

    // write back to fields so PF lattice reads VELOCITY
    cell.template get<VELOCITY<T, LatSet::d>>() = u;
    cell.template get<GenericRho>() = p_cs2;

    // equilibrium
    std::array<T, LatSet::q> geq{};
    EquilibriumScheme::apply(geq, p_cs2, u);

    // Guo discrete force
    std::array<T, LatSet::q> Fi{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      const auto& ci = latset::c<LatSet>(i);
      const T ci_dot_u = ci * u;
      Fi[i] = latset::w<LatSet>(i) * F *
              ((ci - u) * LatSet::InvCs2 +
               ci_dot_u * LatSet::InvCs4 * ci);
    }

    // per-cell omega from OMEGA field (spatially varying)
    const T omega = cell.template get<OMEGA<T>>();
    const T omega_bar = T{1} - omega;
    const T fomega = omega_bar + T{0.5} * omega;  // = 1 - ω/2

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * geq[i] + omega_bar * cell[i] + fomega * Fi[i];
    }
  }
};

}  // namespace fflbm
