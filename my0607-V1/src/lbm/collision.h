// collision.h

#pragma once

#include "lbm/moment.ur.h"

namespace collision {

template <typename MomentaScheme, typename EquilibriumScheme>
struct BGK {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    MomentaScheme::apply(cell, rho, u);
    // moment::template rhoU<CELL, WriteToField>::apply(cell);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);
    // EquilibriumScheme::apply(feq, cell);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

};

template <typename MomentaScheme, typename EquilibriumScheme, typename ForceScheme>
struct BGKForce {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    const auto force = ForceScheme::getForce(cell);
    MomentaScheme::apply(cell, force, rho, u);
    // compute force term
    std::array<T, LatSet::q> fi{};
    ForceScheme::apply(u, force, fi);
    // equilibrium distribution function
    T p = cell.template get<PRESSURE<T>>();
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, p ,u);
    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    if constexpr (CELL::template hasField<PRESSURE<T>>()) {
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
      }
    } else {
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
      }
    }
  }
};

// a typical BGK collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
// force term
template <typename EquilibriumScheme, typename SOURCE, bool WriteToField = false>
struct BGKSource_Feq_Rho {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // update macroscopic variables
    T rho{};
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const auto source = cell.template get<SOURCE>();
    moment::template sourcerho<CELL, SOURCE, WriteToField>::apply(cell, source, rho);
    // equilibrium distribution function
    std::array<T, LatSet::q> feq{};
    EquilibriumScheme::apply(feq, rho, u);

    // BGK collision
    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * source * latset::w<LatSet>(i);
    }
  }
};

// full way bounce back, could be regarded as a mpdified collision process
// swap the populations in the opposite direction
template <typename CELLTYPE>
struct BounceBack {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  __any__ static void apply(CELL& cell) {
    for (unsigned int i = startdir; i < LatSet::q; i += 2) {
      const T temp = cell[i];
      const unsigned int iopp = i + 1;
      cell[i] = cell[iopp];
      cell[iopp] = temp;
    }
  }

};

  // static inline void apply(CELL &cell, unsigned int k) {
  //   cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) +
  //             2 * LatSet::InvCs2 * latset::w<LatSet>(k) * cell.template get<GenericRho>() *
  //               (cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(k));
  // }

// full way bounce back with moving wall, could be regarded as a modified collision process
// swap the populations in the opposite direction
template <typename CELLTYPE>
struct BounceBackMovingWall {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;
  static constexpr unsigned int startdir = LatSet::q % 2 == 0 ? 0 : 1;

  __any__ static void apply(CELL& cell) {
    const T rhox = 2 * LatSet::InvCs2 * cell.template get<GenericRho>();
    for (unsigned int i = startdir; i < LatSet::q; i += 2) {
      const T temp = cell[i];
      const unsigned int iopp = i + 1;
      const T uc = cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(i) * latset::w<LatSet>(i) * rhox;
      cell[i] = cell[iopp] + uc;
      cell[iopp] = temp - uc;
    }
  }

};

// PeriodicBoundary: no-op collision operator for periodic boundary cells
//
// Design rationale:
// - Periodic boundary cells should NOT be modified during the collision phase
// - Their distribution function values are determined by:
//   1. The streaming step (populations move to neighbors)
//   2. FixedPeriodicBoundaryManager::Apply() (copies data from opposite side)
// - If collision modifies these cells (e.g., local swap = BounceBack),
//   it would corrupt the distribution functions and cause artifacts
//
// Usage:
//   using PeriodicTask = tmp::Key_TypePair<PeriodicFlag, collision::PeriodicBoundary<CELL>>;
//   using TaskCollection = tmp::TupleWrapper<BulkTask, PeriodicTask>;
//
// Must be used together with FixedPeriodicBoundaryManager for actual data transfer

template <typename CELLTYPE>
struct PeriodicBoundary {
  using CELL = CELLTYPE;
  using LatSet = typename CELL::LatticeSet;
  using T = typename CELL::FloatType;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    moment::template rho<CELL, true>::apply(cell);
  }

};




template <typename EquilibriumScheme,typename INTERFACEWIDTH, bool WriteToField = false>
struct BGKSourceWithInterface {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell) {
    T phi = T{0};
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    moment::template rho<CELL,WriteToField>::apply(cell);

    Vector<T, LatSet::d> grad;
    Vector<T, LatSet::d> n;

    if (LatSet::d == 3){
      grad[0] = T{0};
      grad[1] = T{0};
      grad[2] = T{0};

      for (unsigned int k = 0; k < LatSet::q; ++k) {
        const auto& neighbor_pos = cell.getNeighbor(k);
        const auto& neighbor_neg = cell.getNeighbor(LatSet::opp[k]);
        T phi_pos = neighbor_pos.template get<GenericRho>();
        T phi_neg = neighbor_neg.template get<GenericRho>();

        T wk = latset::w<LatSet>(k);
        const auto& ck = latset::c<LatSet>(k);

        grad[0] += wk * (phi_pos - phi_neg) * ck[0] / (2 * LatSet::cs2);
        grad[1] += wk * (phi_pos - phi_neg) * ck[1] / (2 * LatSet::cs2);
        grad[2] += wk * (phi_pos - phi_neg) * ck[2] / (2 * LatSet::cs2);
      }

      cell.template get<GRAD<T, LatSet::d>>() = grad;

      T grad_mag = grad.getnorm();
      T epsilon = T(1e-10);
      if (grad_mag > epsilon) {
        n[0] = grad[0] / (grad_mag + epsilon);
        n[1] = grad[1] / (grad_mag + epsilon);
        n[2] = grad[2] / (grad_mag + epsilon);
      } else {
        n[0] = T{0};
        n[1] = T{0};
        n[2] = T{0};
      }
    }else if (LatSet::d == 2){
      grad[0] = T{0};
      grad[1] = T{0};

      for (unsigned int k = 0; k < LatSet::q; ++k) {
        const auto& neighbor_pos = cell.getNeighbor(k);
        const auto& neighbor_neg = cell.getNeighbor(LatSet::opp[k]);
        T phi_pos = neighbor_pos.template get<GenericRho>();
        T phi_neg = neighbor_neg.template get<GenericRho>();

        T wk = latset::w<LatSet>(k);
        const auto& ck = latset::c<LatSet>(k);

        grad[0] += wk * (phi_pos - phi_neg) * ck[0] / (2 * LatSet::cs2);
        grad[1] += wk * (phi_pos - phi_neg) * ck[1] / (2 * LatSet::cs2);
      }
      cell.template get<GRAD<T, LatSet::d>>() = grad;

      T grad_mag = grad.getnorm();
      T epsilon = T(1e-10);
      Vector<T, LatSet::d> n;
      if (grad_mag > epsilon) {
        n[0] = grad[0] / (grad_mag + epsilon);
        n[1] = grad[1] / (grad_mag + epsilon);
      } else {
        n[0] = T{0};
        n[1] = T{0};
      }
    }

    phi = cell.template get<GenericRho>();
    std::array<T, LatSet::q> feq;
    EquilibriumScheme::apply(feq, phi, u);
    T omega_phi = cell.getOmega();
    T _omega_phi = cell.get_Omega();
    T factor = cell.getfOmega();
    T interfacewidth = cell.template get<INTERFACEWIDTH>();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = T{};
      Vector<T, LatSet::d> c = latset::c<LatSet>(i);
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        en += c[d] * n[d];
      } 
      T source = factor * (latset::w<LatSet>(i) * en * (T{4} * phi *(T{1} - phi) )) / interfacewidth;
      cell[i] = omega_phi * feq[i] + _omega_phi * cell[i] + source;
    }
  }
};

template <typename EquilibriumScheme, typename NORMAL, bool WriteToField = false>
struct BGKSource {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell) {
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL>();

    T phi = cell.template get<GenericRho>();
    std::array<T, LatSet::q> feq;
    EquilibriumScheme::apply(feq, phi, u);
    T omega_phi = cell.getOmega();
    T _omega_phi = cell.get_Omega();
    T factor = cell.getfOmega();
    T interfacewidth = cell.template get<INTERFACEWIDTH<T>>();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = T{};
      Vector<T, LatSet::d> c = latset::c<LatSet>(i);
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        en += c[d] * n[d];
      } 
      T source = factor * (latset::w<LatSet>(i) * en * (T{4} * phi *(T{1} - phi) )) / interfacewidth;
      cell[i] = omega_phi * feq[i] + _omega_phi * cell[i] + source;
    }
  }
};

}  // namespace collision

// old version of BGK collision

namespace legacy {

template <typename T, typename LatSet>
struct BGK {
  // BGK collision operator
  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void apply(PopCell<T, LatSet>& cell) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i];
    }
  }

  // BGK collision operator with force
  // template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  // __any__ static void applyForce(PopCell<T, LatSet>& cell, const Vector<T, LatSet::d>& force) {
  //   std::array<T, LatSet::q> feq{};
  //   GetFeq(feq, cell.getVelocity(), cell.getRho());

  //   std::array<T, LatSet::q> fi{};
  //   force::ForcePop<T, LatSet>::compute(fi, cell.getVelocity(), force);

  //   const T omega = cell.getOmega();
  //   const T _omega = cell.get_Omega();
  //   const T fomega = cell.getfOmega();

  //   for (int i = 0; i < LatSet::q; ++i) {
  //     cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
  //   }
  // }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(PopCell<T, LatSet>& cell, const std::array<T, LatSet::q>& fi) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * fi[i];
    }
  }

  template <void (*GetFeq)(std::array<T, LatSet::q>&, const Vector<T, LatSet::d>&, T)>
  __any__ static void applySource(PopCell<T, LatSet>& cell, const T S) {
    std::array<T, LatSet::q> feq{};
    GetFeq(feq, cell.getVelocity(), cell.getRho());

    const T omega = cell.getOmega();
    const T _omega = cell.get_Omega();
    const T fomega = cell.getfOmega();

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      cell[i] = omega * feq[i] + _omega * cell[i] + fomega * S * latset::w<LatSet>(i);
    }
  }
};

}  // namespace collision