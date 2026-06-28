// equlibrium.h

#pragma once
// lattice boltzmann method implementations

// #include "data_struct/Vector.h"
#include "lbm/lattice_set.h"

template <typename T, typename LatSet, typename TypePack>
class Cell;

namespace cudev {
  template <typename T, typename LatSet, typename TypePack>
  class Cell;
}  // namespace cudev

// calc equilibrium distribution function
// sum(feq_i) = rho, for both first and second order
template <typename T, typename LatSet>
struct Equilibrium {
  __any__ static inline T Order1(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order1_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho) {
    return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }

  __any__ static inline T Order2(int k, const Vector<T, LatSet::d> &u, T rho, T u2) {
    T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T(1) + LatSet::InvCs2 * uc + uc * uc * T(0.5) * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T(0.5));
  }

  // __any__ static inline T Order2_Incompresible(int k, const Vector<T, LatSet::d> &u, T rho, T
  // u2) {
  //   T uc = u * latset::c<LatSet>(k);
  //   return latset::w<LatSet>(k) * (rho + LatSet::InvCs2 * uc + uc * uc * T(0.5) *
  //   LatSet::InvCs4 -
  //                          LatSet::InvCs2 * u2 * T(0.5));
  // }

  __any__ static void FirstOrder_Incompresible(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order1_Incompresible(k, u, rho);
    }
  }
  __any__ static void SecondOrder(std::array<T, LatSet::q> &feq, const Vector<T, LatSet::d> &u, T rho) {
    T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = Order2(k, u, rho, u2);
    }
  }
};

namespace equilibrium {

template <typename CELL>
struct SecondOrderImpl {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    const T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      const T uc = u * latset::c<LatSet>(k);
      feq[k] = latset::w<LatSet>(k) * rho *
      (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
    }
  }
};

template <typename CELL>
struct SecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d> &u, const T rho, const T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * rho *
           (T{1} + LatSet::InvCs2 * uc + uc * uc * T{0.5} * LatSet::InvCs4 -
            LatSet::InvCs2 * u2 * T{0.5});
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, rho, u);
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    SecondOrderImpl<CELL>::apply(feq, rho, u);
  }
};

// He-Luo incompressible equilibrium (for variable-density two-phase flow)
// f_i^eq = w_i * p + w_i * [InvCs2*(c_i·u) + 0.5*InvCs4*(c_i·u)^2 - 0.5*InvCs2*u^2]
// where p = sum(f_i) is pressure (not density).
// The actual density rho comes separately from the phase-field (phi interpolation).
template <typename CELL>
struct IncompressibleSecondOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d>& u, const T p, const T u2) {
    const T uc = u * latset::c<LatSet>(k);
    return latset::w<LatSet>(k) * (p + LatSet::InvCs2 * uc
           + uc * uc * T{0.5} * LatSet::InvCs4 - LatSet::InvCs2 * u2 * T{0.5});
  }

  __any__ static void apply(std::array<T, LatSet::q>& feq, const T p,
                             const Vector<T, LatSet::d>& u) {
    const T u2 = u.getnorm2();
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = get(k, u, p, u2);
    }
  }
};


template <typename CELL>
struct FirstOrderImpl {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
    }
  }
};

template <typename CELL>
struct FirstOrder {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using CELLTYPE = CELL;
  using GenericRho = typename CELL::GenericRho;

  __any__ static inline T get(unsigned int k, const Vector<T, LatSet::d> &u, const T rho) {
    return latset::w<LatSet>(k) * rho * (T{1} + LatSet::InvCs2 * (u * latset::c<LatSet>(k)));
  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const CELL &cell) {
    const T rho = cell.template get<GenericRho>();
    const Vector<T, LatSet::d> &u = cell.template get<VELOCITY<T, LatSet::d>>();
  
    apply(feq, rho, u);

  }
  __any__ static void apply(std::array<T, LatSet::q> &feq, const T rho, const Vector<T, LatSet::d> &u) {
      FirstOrderImpl<CELL>::apply(feq, rho, u);
  }
};

// Init feq
// common use: Init<moment::useFieldrhoU<CELL>, equilibrium::SecondOrder<CELL>>
template <typename MomentaScheme, typename EquilibriumScheme>
struct Init {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // macroscopic variables, usually from field value
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = feq[i];
    }
  }

};

template <typename CELLTYPE>
struct PhaseFieldEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  
  __any__ static inline void apply(std::array<T, LatSet::q>& feq, T phi, 
                                   const Vector<T, LatSet::d>& u) {
    T usqr = u[0]*u[0] + u[1]*u[1];
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      T edotu = T{};
      for (unsigned int l = 0; l < LatSet::d; ++l) {
        edotu += latset::c<LatSet>(k)[l] * u[l];
      }
      feq[k] = latset::w<LatSet>(k) * phi * (T{1} + LatSet::InvCs2 * edotu);
    }
  }
  
  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& feq) {
    T phi = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      phi += cell[i];
    }
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    apply(feq, phi, u);
  }
};


// Magnetic field equilibrium (paper Eq. 40):
//   h_α^eq = ω_α^mag * ψ,  where ψ = Σ h_α
// For D2Q5: all ω_α = 0.2, independent of the lattice quadrature weights.
// The magnetic equation is a pure diffusion equation —
// no advection term, no velocity dependence in the equilibrium.
template <typename CELLTYPE>
struct MagEquilibrium {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  // Compute equilibrium from scalar psi
  __any__ static inline void apply(std::array<T, LatSet::q>& feq, const T psi) {
    constexpr T w_eq = T{1} / T{LatSet::q};  // = 0.2 for D2Q5, 1/7 for D3Q7
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      feq[k] = w_eq * psi;
    }
  }

  // Compute equilibrium from cell populations (psi = Σ cell[k])
  __any__ static inline void apply(CELL& cell, std::array<T, LatSet::q>& feq) {
    T psi = T{};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      psi += cell[i];
    }
    apply(feq, psi);
  }
};

}  // namespace equilibrium