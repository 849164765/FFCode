// coupling.h
// Three-solver coupling operators (dual-cell interface)
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025)
//
// Each operator reads from one solver's cell and writes to another's.
// Interface: apply(CELL0&, CELL1&) — supported by BlockLatManagerCoupling.

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// PFtoNS_Properties (Eqs.49-50): PF → NS material properties
//   rho = rho_l + phi * (rho_h - rho_l)
//   eta = eta_l + phi * (eta_h - eta_l)
//   rho_l, rho_h, eta_l, eta_h are Data<T,Base> constants registered on NS cell
// ================================================================
template <typename CELL_PF_TYPE, typename CELL_NS_TYPE>
struct PFtoNS_Properties {
  using CELL_PF = CELL_PF_TYPE;
  using CELL_NS = CELL_NS_TYPE;
  using T = typename CELL_PF::FloatType;

  __any__ static inline void apply(CELL_PF& pf, CELL_NS& ns) {
    T phi = pf.template get<PHI<T>>();

    T rho_l = ns.template get<RHOL<T>>();
    T rho_h = ns.template get<RHOH<T>>();
    T eta_l = ns.template get<ETAL<T>>();
    T eta_h = ns.template get<ETAH<T>>();

    // Eq.49: rho(phi)
    ns.template get<RHO<T>>() = rho_l + phi * (rho_h - rho_l);

    // Eq.50: eta(phi)
    ns.template get<VISCOSITY<T>>() = eta_l + phi * (eta_h - eta_l);
  }
};

// ================================================================
// PFtoNS_Forces (Eq.12): PF → NS surface tension force
//   F_s = lambda_phi * grad_phi
// Accumulates into NS FORCE field (caller zeroes FORCE before first call).
// ================================================================
template <typename CELL_PF_TYPE, typename CELL_NS_TYPE>
struct PFtoNS_Forces {
  using CELL_PF = CELL_PF_TYPE;
  using CELL_NS = CELL_NS_TYPE;
  using T = typename CELL_PF::FloatType;
  static constexpr unsigned int D = CELL_PF::LatticeSet::d;

  __any__ static inline void apply(CELL_PF& pf, CELL_NS& ns) {
    T lambda = pf.template get<PHASECHEMPOTENTIAL<T>>();
    const auto& grad = pf.template get<GRAD<T, D>>();
    ns.template get<FORCE<T, D>>() += lambda * grad;
  }
};

// ================================================================
// PFtoMag_Permeability (Eq.51): PF → Mag magnetic permeability
//   mu = mu_l + phi * (mu_h - mu_l)
//   mu_l, mu_h are Data<T,Base> constants registered on Mag cell
// ================================================================
template <typename CELL_PF_TYPE, typename CELL_MAG_TYPE>
struct PFtoMag_Permeability {
  using CELL_PF = CELL_PF_TYPE;
  using CELL_MAG = CELL_MAG_TYPE;
  using T = typename CELL_PF::FloatType;

  __any__ static inline void apply(CELL_PF& pf, CELL_MAG& mag) {
    T phi = pf.template get<PHI<T>>();

    T mu_l = mag.template get<MUL<T>>();
    T mu_h = mag.template get<MUH<T>>();

    mag.template get<MAGPERMEABILITY<T>>() = mu_l + phi * (mu_h - mu_l);
  }
};

// ================================================================
// MagtoNS_Force (Eq.8): Mag → NS magnetic force
//   F_total += F_m  (F_m computed by MagForce, stored on Mag cell FORCE)
// ================================================================
template <typename CELL_MAG_TYPE, typename CELL_NS_TYPE>
struct MagtoNS_Force {
  using CELL_MAG = CELL_MAG_TYPE;
  using CELL_NS = CELL_NS_TYPE;
  using T = typename CELL_MAG::FloatType;
  static constexpr unsigned int D = CELL_MAG::LatticeSet::d;

  __any__ static inline void apply(CELL_MAG& mag, CELL_NS& ns) {
    const auto& Fm = mag.template get<FORCE<T, D>>();
    ns.template get<FORCE<T, D>>() += Fm;
  }
};

// ================================================================
// NStoPF_Velocity: NS → PF velocity coupling
//   Copies VELOCITY from NS cell to PF cell for convection term
// ================================================================
template <typename CELL_NS_TYPE, typename CELL_PF_TYPE>
struct NStoPF_Velocity {
  using CELL_NS = CELL_NS_TYPE;
  using CELL_PF = CELL_PF_TYPE;
  using T = typename CELL_NS::FloatType;
  static constexpr unsigned int D = CELL_NS::LatticeSet::d;

  __any__ static inline void apply(CELL_NS& ns, CELL_PF& pf) {
    const auto& u = ns.template get<VELOCITY<T, D>>();
    pf.template get<VELOCITY<T, D>>() = u;
  }
};
