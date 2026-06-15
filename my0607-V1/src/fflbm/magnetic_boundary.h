// magnetic_boundary.h
// Neumann boundary condition for D2Q5 magnetic solver
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025), Eqs.45-48
//
// For each boundary cell, fills populations pointing into the domain
// using the non-equilibrium extrapolation / virtual node scheme.

#pragma once

#include "lbm/lattice_set.h"
#include "data_struct/Vector.h"
#include "utils/alias.h"

// ================================================================
// MagNeumannBC (Eqs.45-48)
// Only applies to boundary cells (FLAG != 0).
// For each outward direction k, fills opp(k) using interior neighbor:
//   h_k(x0) = h_k(x1)                       (conserved, Eq.45-46)
//   h_2(x0) = (1-1/c)*h_2(x1) + (1/c)*h_4(x1) + eps*mu*H_n/c  (Eq.47)
//   h_4(x0) = -(1/c)*h_2(x1) + (1+1/c)*h_4(x1) + eps*mu*H_n/c (Eq.48)
//
// D2Q5 c[]: [0]=(0,0), [1]=(1,0), [2]=(-1,0), [3]=(0,1), [4]=(0,-1)
// ================================================================
template <typename CELLTYPE>
struct MagNeumannBC {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static inline void apply(CELL& cell) {
    // Only applicable to boundary cells (FLAG != 0)
    if (cell.template get<FLAG>() == 0) return;

    // For each direction k != 0:
    // if neighbor in direction k is outside, fill opp(k) using Eqs.47-48
    for (int k = 1; k < LatSet::q; ++k) {
      auto neighbor = cell.getNeighbor(k);
      if (neighbor.template get<FLAG>() == 0) continue;

      int kin = latset::opp<LatSet>(k);  // direction pointing into domain
      auto cell_x1 = cell.getNeighbor(kin);

      // CFL = dx/dt = 1 in standard lattice units
      T cfl = T{1};

      // eps_mu_Hn = epsilon * mu * H_n
      //   epsilon = 1/mu0 = 1 (lattice units)
      //   mu = MAGPERMEABILITY (from cell)
      //   H_n = prescribed normal field (zero for insulating boundary)
      T mu = cell.template get<MAGPERMEABILITY<T>>();
      T H_n = T{0};  // TODO: read from config parameter
      T eps_mu_Hn = mu * H_n;

      // Eqs.47-48: virtual node populations
      int ktan = _tangential(kin);
      T h_in_x1 = cell_x1[kin];
      T h_tan_x1 = cell_x1[ktan];
      T invC = T{1} / cfl;

      // Eq.47: incoming direction
      cell[kin] = (T{1} - invC) * h_in_x1 + invC * h_tan_x1
                + eps_mu_Hn * invC;
      // Eq.48: tangential direction
      cell[ktan] = -invC * h_in_x1 + (T{1} + invC) * h_tan_x1
                  + eps_mu_Hn * invC;
      break;  // one boundary face per cell
    }
  }

 private:
  // D2Q5 c[]: [1]=(1,0), [2]=(-1,0), [3]=(0,1), [4]=(0,-1)
  // tangential: swap x↔y
  __any__ static constexpr int _tangential(int k) {
    switch (k) {
      case 1: return 3;   // (1,0) → (0,1)
      case 2: return 4;   // (-1,0) → (0,-1)
      case 3: return 1;   // (0,1) → (1,0)
      case 4: return 2;   // (0,-1) → (-1,0)
      default: return 0;
    }
  }
};
