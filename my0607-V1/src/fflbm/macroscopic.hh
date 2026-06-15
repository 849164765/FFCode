// macroscopic.hh — compute macroscopic variables from distribution functions
//
// Operators:
//   RhoU<CELL>      : rho = Σf_i,  u = Σc_i·f_i / rho  (no external force)
//   ForceRhoU<CELL> : rho = Σf_i,  u = (Σc_i·f_i + F/2) / rho  (Guo force correction)

#pragma once

namespace fflbm {

template <typename CELL>
struct RhoU {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  // compute rho and u, optionally write back to cell fields
  __any__ static void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T{0};
    u.clear();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T f = cell[k];
      rho += f;
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * f;
      }
    }
    if (rho > T{0}) {
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] /= rho;
      }
    }
  }

  // convenience: read-only from cell, don't write back
  __any__ static void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u,
                            bool /*write_back*/) {
    apply(cell, rho, u);
  }
};

// Force-corrected moment: u = (Σc_i·f_i + F/2) / ρ
// This is the standard Guo (2002) force-corrected velocity for BGK.
// The force F is read from the FORCE<T,d> field of the cell.
template <typename CELL>
struct ForceRhoU {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell, T& rho, Vector<T, LatSet::d>& u) {
    rho = T{0};
    u.clear();

    // sum of populations
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T f = cell[k];
      rho += f;
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * f;
      }
    }

    // Guo half-force correction:  u = (Σc·f + F/2) / ρ
    const auto& force = cell.template get<FORCE<T, LatSet::d>>();
    T half_rho = (rho > T{0}) ? (T{0.5} / rho) : T{0};
    for (unsigned int d = 0; d < LatSet::d; ++d) {
      u[d] = u[d] / rho + force[d] * half_rho;
    }
  }

  // write-back variant — writes rho to GenericRho and u to VELOCITY
  __any__ static void applyWriteBack(CELL& cell) {
    T rho;
    Vector<T, LatSet::d> u;
    apply(cell, rho, u);
    cell.template get<GenericRho>() = rho;
    cell.template get<VELOCITY<T, LatSet::d>>() = u;
  }
};

// Pressure-based moment (g-LBE):
//   p̃ = p/cs² = Σg_i
//   u = Σc_i·g_i + F·Δt/(2ρ₀)  (Guo half-force correction)
// ρ₀ is the reference density (constant, typically average density).
// Set PressureMoment<CELL>::rho0 before simulation starts.
template <typename CELL>
struct PressureMoment {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static inline T rho0 = T{1.0};  // set before use

  __any__ static void apply(CELL& cell, T& p_cs2, Vector<T, LatSet::d>& u) {
    p_cs2 = T{0};
    u.clear();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T g = cell[k];
      p_cs2 += g;
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * g;
      }
    }
    // Guo half-force correction: u += F·Δt/(2·ρ₀)
    const auto& F = cell.template get<FORCE<T, LatSet::d>>();
    if (rho0 > T{0}) {
      T half = T{0.5} / rho0;
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += F[d] * half;
      }
    }
  }
};

}  // namespace fflbm
