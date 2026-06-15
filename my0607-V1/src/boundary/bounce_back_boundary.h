// bounceback like boundary condition

#pragma once

#include "boundary/basic_boundary.h"
#include "lbm/equilibrium.ur.h"

// --------------------------------------------------------------------------------------
// ---- BlockBoundary operators (CELL-compatible) ------------
// --------------------------------------------------------------------------------------

namespace bounceback {

template <typename CELL>
struct normal {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k));
  }

  static inline void apply(CyclicArray<T>& arr, const CyclicArray<T>& arrk, std::size_t id) {
    arr[id] = arrk.getPrevious(id);
  }
};

template <typename CELL>
struct anti_simple {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * cell.template get<GenericRho>() * latset::w<LatSet>(k) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};
template <typename CELL>
struct anti_O1 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * equilibrium::FirstOrder<CELL>::get(
                    k, cell.template get<VELOCITY<T, LatSet::d>>(),
                    cell.template get<GenericRho>()) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};
template <typename CELL>
struct anti_O2 {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = 2 * equilibrium::SecondOrder<CELL>::get(
                    k, cell.template get<VELOCITY<T, LatSet::d>>(),
                    cell.template get<GenericRho>(),
                    cell.template get<VELOCITY<T, LatSet::d>>().getnorm2()) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};
template <typename CELL>
struct anti_pressure {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    // get the interpolated velocity
    const Vector<T, LatSet::d> uwall =
      cell.template get<VELOCITY<T, LatSet::d>>() +
      T{0.5} * (cell.template get<VELOCITY<T, LatSet::d>>() -
                cell.getNeighbor(latset::opp<LatSet>(k)).template get<VELOCITY<T, LatSet::d>>());
    cell[k] = 2 * cell.template get<GenericRho>() * latset::w<LatSet>(k) *
                (T{1} + std::pow((uwall * latset::c<LatSet>(k)), 2) * T{0.5} * LatSet::InvCs4 -
                 uwall.getnorm2() * T(0.5) * LatSet::InvCs2) -
              cell.getPrevious(latset::opp<LatSet>(k));
  }
};

template <typename CELL>
struct movingwall {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  static inline void apply(CELL &cell, unsigned int k) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) +
              2 * LatSet::InvCs2 * latset::w<LatSet>(k) * cell.template get<GenericRho>() *
                (cell.template get<VELOCITY<T, LatSet::d>>() * latset::c<LatSet>(k));
  }
  static inline void apply(CELL &cell, unsigned int k, const Vector<T, LatSet::d> &wall_velocity) {
    cell[k] = cell.getPrevious(latset::opp<LatSet>(k)) - 2 * LatSet::InvCs2 * latset::w<LatSet>(k) *
                                                   cell.template get<GenericRho>() *
                                                   (wall_velocity * latset::c<LatSet>(k));
  }
};

}  // namespace bounceback
