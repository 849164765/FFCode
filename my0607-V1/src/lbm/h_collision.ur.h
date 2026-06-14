// h_collision.ur.h
// Phase B: Interface-tracking h-LBE (Eqs. 6-12)
// from Fakhari et al., Phys. Rev. E 96, 053301 (2017)
//
// B1-B4: H_Collision — BGK collision for h distribution
//   h_α(x+e_α) = h_α(x) - (h_α - h̄^eq_α)/(τ_φ+½) + F^φ_α      Eq.6
//   F^φ_α = [1-4(φ-φ₀)²]/ξ · w_α · e_α·n                     Eq.7
//   h^eq_α = φ · Γ_α                                           Eq.10 via Γ
//   h̄^eq_α = h^eq_α - ½·F^φ_α                                  Eq.9
//
// B5: UpdatePhaseField — φ = Σ_α h_α                          Eq.12
//
// Note: Streaming is done separately by BlockLattice::Stream()

#pragma once

#include "data_struct/cell.h"
#include "lbm/lattice_set.h"
#include "lbm/equilibrium.ur.h"

namespace phase_field {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;
#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// ===================================================================
// H_Collision — h-LBE BGK collision with phase-field forcing
//
// Reads:  PHI<T>, NORMAL<T,D>, INTERFACEWIDTH<T>, VELOCITY<T,D>
//         h_α populations (via cell[i] = POP)
// Writes: h_α populations (collided values, ready for streaming)
//
// Omega for collision: cell.getOmega() / cell.get_Omega() / cell.getfOmega()
//   omega   = 1/(τ_φ+½)     → weight for equilibrium
//   _omega  = 1 - omega      → weight for old population
//   fomega  = τ_φ/(τ_φ+½)    → weight for forcing term
// ===================================================================

template <typename CELLTYPE>
struct H_Collision {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr T phi0 = T{0.5};  // φ₀

  __any__ static inline void apply(CELL& cell) {
    const T phi = cell.template get<PHI<T>>();
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL<T, LatSet::d>>();
    const T xi = cell.template get<INTERFACEWIDTH<T>>();
    const T omega_phi = cell.getOmega();
    const T _omega_phi = cell.get_Omega();
    const T fomega_phi = cell.getfOmega();

    // ---- Precompute Γ_α (Eq.10) for all directions ----
    std::array<T, LatSet::q> Gamma;
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const auto& e_alpha = latset::c<LatSet>(alpha);
      const T w_alpha = latset::w<LatSet>(alpha);
      T cu = T{0};
      T u2 = T{0};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cu += e_alpha[d] * u[d];
        u2 += u[d] * u[d];
      }
      Gamma[alpha] = w_alpha * (T{1}
        + cu * LatSet::InvCs2
        + cu * cu * LatSet::InvCs4 * T{0.5}
        - u2 * LatSet::InvCs2 * T{0.5});
    }

    // ---- B1: Compute forcing term F^φ_α (Eq.7) ----
    // F^φ_α = [1 - 4(φ-φ₀)²] / ξ · w_α · e_α·n
    // Note: 1 - 4(φ-φ₀)² = 4φ(1-φ)
    const T phi_factor = T{4} * phi * (T{1} - phi) / xi;

    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      const auto& e_alpha = latset::c<LatSet>(alpha);
      T en = T{0};
      for (unsigned int d = 0; d < LatSet::d; ++d) en += e_alpha[d] * n[d];
      const T F_phi_alpha = phi_factor * latset::w<LatSet>(alpha) * en;

      // ---- B2-B3: equilibrium and modified equilibrium ----
      const T h_eq_alpha = phi * Gamma[alpha];

      // ---- B4: BGK collision (Eq.6) ----
      // h_new = omega * h_eq + _omega * h_old + fomega * F^φ
      cell[alpha] = omega_phi * h_eq_alpha
                  + _omega_phi * cell[alpha]
                  + fomega_phi * F_phi_alpha;
    }
  }
};


// ===================================================================
// UpdatePhaseField — φ = Σ h_α  (Eq.12, after streaming)
// Reads:  h_α populations (via cell[i])
// Writes: PHI<T>
// ===================================================================

template <typename CELLTYPE>
struct UpdatePhaseField {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    T phi = T{0};
    for (unsigned int alpha = 0; alpha < LatSet::q; ++alpha) {
      phi += cell[alpha];
    }
    cell.template get<PHI<T>>() = phi;
  }
};

}  // namespace phase_field
