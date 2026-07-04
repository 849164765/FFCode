// collisionMRT.h

#pragma once

#include "lbm/moment.ur.h"

// Forward declare magnetic field base types (defined in ff/ff2d.h)
namespace ff { struct MagOmegaBase; }

namespace mrtdata {
// mrt transformation matrix: mapping from population to moment space
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> M[Q][Q] = {};

// Inverse of M
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> InvM[Q][Q] = {};

// relaxation times
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s[Q] = {};

// relaxation times
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s2[Q] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ int shearIndexes = {};

// relevant indexes of r. t. for shear viscosity
template <unsigned int D, unsigned int Q>
__constexpr__ int shearViscIndexes[shearIndexes<D,Q>] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s_phase_source[Q] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s1[Q] = {};

// mrt transformation matrix
template <>
__constexpr__ Fraction<> M<2,9>[9][9] = {
{ 1,  1,  1,  1,  1,  1,  1,  1,  1},
{-4, -1, -1, -1, -1,  2,  2,  2,  2},
{ 4, -2, -2, -2, -2,  1,  1,  1,  1},
{ 0,  1, -1,  0,  0,  1, -1,  1, -1},
{ 0, -2,  2,  0,  0,  1, -1,  1, -1},
{ 0,  0,  0,  1, -1,  1, -1, -1,  1},
{ 0,  0,  0, -2,  2,  1, -1, -1,  1},
{ 0,  1,  1, -1, -1,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  1,  1, -1, -1}
};

template <>
__constexpr__ Fraction<> InvM<2,9>[9][9] = {
  {{1, 9}, {-1,  9}, { 1,  9},       0,        0,       0,        0,       0,       0},
  {{1, 9}, {-1, 36}, {-1, 18}, { 1, 6}, {-1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18}, {-1, 6}, { 1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, { 1, 6}, {-1,  6}, {-1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, {-1, 6}, { 1,  6}, {-1, 4},       0},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, { 1, 6}, { 1, 12},       0, { 1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, {-1, 6}, {-1, 12},       0, { 1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, {-1, 6}, {-1, 12},       0, {-1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, { 1, 6}, { 1, 12},       0, {-1, 4}}
};

template <>
__constexpr__ Fraction<> s<2,9>[9] = {
0, {11, 10}, {11, 10}, 0, {11, 10}, 0, {11, 10}, 0, 0
};

template <>
__constexpr__ Fraction<> s_phase_source<2,9>[9] = {
  0, {11, 10}, {11, 10}, 0, {11, 10}, 0, {11, 10}, 1, 1
};

template <>
__constexpr__ Fraction<> s1<2,9>[9] = {
  1, {14, 10}, {14, 10}, 1, {12, 10}, 1, {12, 10}, 1, 1
};

template <>
__constexpr__ int shearIndexes<2,9> = 2;

template <>
__constexpr__ int shearViscIndexes<2,9>[shearIndexes<2,9>] = { 7, 8};

// ===================================================================
// D2Q5 MRT — magnetic field potential solver
// Hu & Li (2018) Phys. Rev. E 98, 033301, §II.D
// Moments: m = (ψ, Jx, Jy, e, pxx)
//   ψ   = Σ h_i         (conserved)
//   Jx  = h_1 - h_2     (flux x)
//   Jy  = h_3 - h_4     (flux y)
//   e   = 4h_0 - Σ_{i>0} h_i  (energy)
//   pxx = (h_1+h_2) - (h_3+h_4) (normal stress diff)
// ===================================================================
template <>
__constexpr__ Fraction<> M<2,5>[5][5] = {
  { 1,  1,  1,  1,  1},
  { 0,  1, -1,  0,  0},
  { 0,  0,  0,  1, -1},
  { 4, -1, -1, -1, -1},
  { 0,  1,  1, -1, -1}
};

template <>
__constexpr__ Fraction<> InvM<2,5>[5][5] = {
  {{1, 5},       0,       0,  {1, 5},       0},
  {{1, 5},  {1, 2},       0, {-1,20},  {1, 4}},
  {{1, 5}, {-1, 2},       0, {-1,20},  {1, 4}},
  {{1, 5},       0,  {1, 2}, {-1,20}, {-1, 4}},
  {{1, 5},       0, {-1, 2}, {-1,20}, {-1, 4}}
};

template <>
__constexpr__ int shearIndexes<2,5> = 0;

template <>
__constexpr__ int shearViscIndexes<2,5>[0] = {};

template <>
__constexpr__ Fraction<> M<3,19>[19][19] = {
{ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{-1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{ 1, -2, -2, -2, -2, -2, -2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{ 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0},
{ 0, -2,  2,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0},
{ 0,  0,  0,  1, -1,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1},
{ 0,  0,  0, -2,  2,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1},
{ 0,  0,  0,  0,  0,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1},
{ 0,  0,  0,  0,  0, -2,  2,  0,  0,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1},
{ 0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1, -2, -2,  1,  1,  1,  1, -2, -2},
{ 0, -2, -2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2,  1,  1,  1,  1, -2, -2},
{ 0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,  0,  0,  1,  1, -1, -1,  0,  0},
{ 0,  0,  0, -1, -1,  1,  1,  1,  1, -1, -1,  0,  0,  1,  1, -1, -1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1},
{ 0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1, -1,  1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1, -1,  1,  0,  0, -1,  1},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0, -1,  1,  1, -1}
};

template <>
__constexpr__ Fraction<> InvM<3,19>[19][19] = {
  {{1, 3}, {-1, 2}, {1, 6},       0,        0,       0,       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18}, {1, 6}, {-1, 6},       0,       0,        0,        0,  {1,12}, {-1,12},        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18}, {-1,6},  {1, 6},       0,       0,        0,        0,  {1,12}, {-1,12},        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0, {1, 6}, {-1,6},        0,        0, {-1,24},  {1,24},  {1, 8}, {-1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0, {-1,6},  {1,6},        0,        0, {-1,24},  {1,24},  {1, 8}, {-1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0,      0,       0,  {1, 6}, {-1, 6}, {-1,24},  {1,24}, {-1, 8},  {1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0,      0,       0, {-1, 6},  {1, 6}, {-1,24},  {1,24}, {-1, 8},  {1, 8},        0,        0,        0,        0,        0,        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24}, {1,12}, {1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16},  {1, 4},        0,        0,  {1, 8},  {1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24}, {-1,12},{-1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16},  {1, 4},        0,        0, {-1, 8}, {-1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24},      0,       0,  {1,12},  {1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0,  {1, 4},        0, {-1, 8},        0,  {1, 8}},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24},      0,       0, {-1,12}, {-1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0,  {1, 4},        0,  {1, 8},        0, {-1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {1,12}, {1,24},  {1,12},  {1,24}, {-1,24}, {-1,24},        0,        0,        0,        0,  {1, 4},        0, {-1, 8}, {-1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {-1,12},{-1,24}, {-1,12}, {-1,24}, {-1,24}, {-1,24},        0,        0,        0,        0,  {1, 4},        0,  {1, 8},  {1, 8}},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24}, {-1,12},{-1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16}, {-1, 4},        0,        0,  {1, 8}, {-1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24}, {1,12}, {1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16}, {-1, 4},        0,        0, {-1, 8},  {1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24},      0,       0, {-1,12}, {-1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0, {-1, 4},        0, {-1, 8},        0, {-1, 8}},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24},      0,       0,  {1,12},  {1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0, {-1, 4},        0,  {1, 8},        0,  {1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {1,12}, {1,24}, {-1,12}, {-1,24}, {-1,24}, {-1,24},        0,        0,        0,        0, {-1, 4},        0, {-1, 8},  {1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {-1,12},{-1,24}, {1,12},  {1,24}, {-1,24}, {-1,24},        0,        0,        0,        0, {-1, 4},        0,  {1, 8}, {-1, 8}}
};

template <>
__constexpr__ Fraction<> s<3,19>[19] = {
0, {119, 100}, {14, 10}, 0, {12, 10}, 0, {12, 10}, 0, {12, 10}, 0, {14, 10}, 0, {14, 10}, 0, 0, 0, {99, 50}, {99, 50}, {99, 50}
};

template <>
__constexpr__ Fraction<> s1<3,19>[19] = {
1, {119, 100}, {14, 10}, 1, {12, 10}, 1, {12, 10}, 1, {12, 10}, 0, {14, 10}, 0, {14, 10}, 0, 0, 0, {99, 50}, {99, 50}, {99, 50}
};

template <>
__constexpr__ int shearIndexes<3,19> = 5;

template <>
__constexpr__ int shearViscIndexes<3,19>[shearIndexes<3,19>] = { 9, 11, 13, 14, 15};

} // namespace mrtdata

namespace mrt {

template <typename LatSet>
constexpr typename LatSet::FloatType M(unsigned int i, unsigned int j) {
  return mrtdata::M<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType InvM(unsigned int i, unsigned int j) {
  return mrtdata::InvM<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType s(unsigned int i) {
  return mrtdata::s<LatSet::d,LatSet::q>[i].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType s1(unsigned int i) {
  return mrtdata::s1<LatSet::d,LatSet::q>[i].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr int shearIndexes() {
  return mrtdata::shearIndexes<LatSet::d,LatSet::q>;
}

template <typename LatSet>
constexpr int shearViscIndexes(unsigned int i) {
  return mrtdata::shearViscIndexes<LatSet::d,LatSet::q>[i];
}

} // namespace mrt

namespace collision {

// a typical MRT collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
template <typename CELL, bool WriteToField = false>
struct MRT_Feq_RhoU {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = equilibrium::SecondOrder<CELL>;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    const T omega = cell.getOmega();
    // relaxation time vector
    T rtvec[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rtvec[i] = mrt::s<LatSet>(i);
    }
    for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i) {
      rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega;
    }
    // relaxation time matrix
    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhoU<CELL, WriteToField>::apply(cell, rho, u);

    // compute Momenta
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }
    // compute Equilibrium
    std::array<T, LatSet::q> feq{};
    equilibrium::SecondOrder<CELL>::apply(feq, rho, u);
    T momentaEq[LatSet::q] {};
    for(unsigned int i = 0; i < LatSet::q; ++i) {
      momentaEq[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
      }
    }
    // MRT collision
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * (momenta[j] - momentaEq[j]);
      }
      cell[i] -= momentaEqS;
    }
  }

};

// MRT source collision for phase-field equations.
//
// Template parameters:
//   UseCHRelaxation = false (default): interface sharpening mode.
//     j-moments use s=1 (full relaxation), suitable for NORMAL=∇φ/|∇φ|
//     (e.g., simpledrop2dMRT).
//   UseCHRelaxation = true: Cahn-Hilliard diffusion mode.
//     j-moments use s=omega_phi (matches BGKSource for correct
//     effective mobility), suitable for NORMAL=∇μ
//     (e.g., rti2dMRT, bubble).
template <typename EquilibriumScheme, typename NORMAL, bool WriteToField = false, bool UseCHRelaxation = false>
struct MRTSource {
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
    T interfacewidth = cell.template get<INTERFACEWIDTH<T>>();

    T rtvec[LatSet::q] {};
    if constexpr (UseCHRelaxation) {
      // Cahn-Hilliard: Fortran-aligned relaxation vector
      if constexpr (std::is_same_v<LatSet, D2Q9<T>> || std::is_same_v<LatSet, D2Q5<T>>) {
        // Fortran Sf = {1.0, 1.1, 1.1, 1/tau, 1/tau, 1/tau, 1/tau, 1.2, 1.2}
        rtvec[0] = T{1};                                    // rho/phi (conserved)
        rtvec[1] = T{11./10.};                               // e
        rtvec[2] = T{11./10.};                               // epsilon
        rtvec[3] = omega_phi;                                // jx (1/tau)
        rtvec[4] = omega_phi;                                // qx (1/tau)
        if constexpr (LatSet::q > 5) {
          rtvec[5] = omega_phi;                              // jy (1/tau)
          rtvec[6] = omega_phi;                              // qy (1/tau)
          rtvec[7] = T{12./10.};                             // pxx (1.2)
          rtvec[8] = T{12./10.};                             // pxy (1.2)
        }
      } else if constexpr (std::is_same_v<LatSet, D3Q19<T>>) {
        // D3Q19 CH relaxation vector aligned with mrtdata::M<3,19> row order:
        // 0:rho, 1:e, 2:epsilon, 3:jx, 4:qx, 5:jy, 6:qy, 7:jz, 8:qz,
        // 9-18: stress/shear and higher-order moments (9,11,13,14,15 are shear modes)
        rtvec[0] = T{1};            // rho/phi (conserved)
        rtvec[1] = T{11./10.};       // e
        rtvec[2] = T{11./10.};       // epsilon
        rtvec[3] = omega_phi;        // jx (1/tau_phi)
        rtvec[4] = omega_phi;        // qx (1/tau_phi)
        rtvec[5] = omega_phi;        // jy (1/tau_phi)
        rtvec[6] = omega_phi;        // qy (1/tau_phi)
        rtvec[7] = omega_phi;        // jz (1/tau_phi)
        rtvec[8] = omega_phi;        // qz (1/tau_phi)
        rtvec[9] = T{12./10.};       // pxx
        rtvec[10] = T{12./10.};      // pww
        rtvec[11] = T{12./10.};      // pxy
        rtvec[12] = T{12./10.};      // pyz
        rtvec[13] = T{12./10.};      // pxz
        rtvec[14] = T{12./10.};      // higher-order stress
        rtvec[15] = T{12./10.};      // higher-order stress
        rtvec[16] = T{12./10.};      // higher-order ghost
        rtvec[17] = T{12./10.};      // higher-order ghost
        rtvec[18] = T{12./10.};      // higher-order ghost
      } else {
        // Generic fallback: use omega_phi for all moments
        for (unsigned int i = 0; i < LatSet::q; ++i) {
          rtvec[i] = omega_phi;
        }
      }
    } else {
      // Interface sharpening: s1 + shear override (original behaviour)
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        rtvec[i] = mrt::s1<LatSet>(i);
      }
      for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i) {
        rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega_phi;
      }
    }

    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }

    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momentaEq[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
      }
    }

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * (momenta[j] - momentaEq[j]);
      }
      cell[i] -= momentaEqS;
    }

    T A = T{4} * phi * (T{1} - phi) / interfacewidth;

    T S_body[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = T{};
      Vector<T, LatSet::d> c = latset::c<LatSet>(i);
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        en += c[d] * n[d];
      }
      S_body[i] = latset::w<LatSet>(i) * en * A;
    }

    T source_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j) {
      source_m[j] = T{};
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        source_m[j] += mrt::M<LatSet>(j, k) * S_body[k];
      }
    }

    T source_m_corrected[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j) {
      source_m_corrected[j] = (T{1} - rtvec[j] / T{2}) * source_m[j];
    }

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        source += mrt::InvM<LatSet>(i, j) * source_m_corrected[j];
      }
      cell[i] += source;
    }
  }

};

// MRT collision with Guo force scheme for Navier-Stokes
// replaces BGKForce<moment::forcerhoU, equilibrium::SecondOrder, force::Force>
template <typename CELL, typename ForceField>
struct MRTForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // Per-cell omega (spatially varying relaxation) for variable viscosity
    T omega{};
    if constexpr (cell.template hasField<OMEGA<T>>()) {
      omega = cell.template get<OMEGA<T>>();
    } else {
      omega = cell.getOmega();
    }
    const auto& F = cell.template get<ForceField>();

    // Compile-time check: incompressible (He-Luo) vs compressible LBM
    // DENSITY field present => He-Luo incompressible (Fortran-aligned)
    constexpr bool isIncompressible = CELL::template hasField<DENSITY<T>>();

    // 1. Compute zeroth moment p=Σf (pressure) and first moment Σ(c·f)
    //    Incompressible: u = Σ(c·f) (velocity, NOT divided by rho)
    //    Compressible:   u = Σ(c·f)/rho (momentum divided by rho)
    T zeroth{};  // p = Σf (pressure)
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      zeroth += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
      }
    }

    // rho: incompressible reads from DENSITY field (phi-based, set by
    // FFRhoOmegaUpdate3D); compressible falls back to rho = Σf
    T rho{};
    if constexpr (isIncompressible) {
      rho = cell.template get<DENSITY<T>>();
    } else {
      rho = zeroth;
      for (unsigned int d = 0; d < LatSet::d; ++d) u[d] /= rho;
    }

    // 2. Half-force correction: u_c = u + F/(2*rho)
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d) uc[d] += F[d] / (T{2} * rho);

    // Write macroscopic fields
    // Incompressible: do NOT write GenericRho (=DENSITY), which is set from phi
    // by FFRhoOmegaUpdate3D and must not be overwritten with Σf.
    // Compressible: write rho = Σf to GenericRho.
    if constexpr (!isIncompressible) {
      cell.template get<GenericRho>() = rho;
    }
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;
    T rtvec[LatSet::q] {};
    if constexpr (std::is_same_v<LatSet, D3Q19<T>>) {
      // =====================================================================
      // 3D Multiphase MRT Relaxation Vector (Based on Premnath 2007 & Fakhari 2017)
      // 核心原则：解耦物理剪切粘度与数值耗散，镇压 3D 界面高频噪声与声波！
      // =====================================================================
      
      // 放弃 2D 魔法公式，固定为常数以提供稳定的界面耗散
      T s_q_fixed = T{1.1}; 
      T s_bulk = T{1.1};     // 控制体粘度，用于快速衰减界面声波
      T s_ghost = T{1.2};    // 幽灵矩强制耗散，充当高频噪声垃圾桶

      rtvec[0] = T{0};       // rho (守恒)
      rtvec[1] = s_bulk;     // e (能量 -> 控制体粘度，必须 > omega！)
      rtvec[2] = s_bulk;     // epsilon (控制体粘度)
      rtvec[3] = T{0};       // jx (守恒)
      rtvec[4] = s_q_fixed;  // qx (能量通量)
      rtvec[5] = T{0};       // jy (守恒)
      rtvec[6] = s_q_fixed;  // qy (能量通量)
      rtvec[7] = T{0};       // jz (守恒)
      rtvec[8] = s_q_fixed;  // qz (能量通量)
      
      // 剪切应力矩 (控制物理剪切粘度，必须严格绑定 omega)
      // 对应库里的 shearViscIndexes: 9, 11, 13, 14, 15
      rtvec[9]  = omega;     // 3cx^2 - r^2 (剪切)
      rtvec[11] = omega;     // cy^2 - cz^2 (剪切)
      rtvec[13] = omega;     // xy (剪切)
      rtvec[14] = omega;     // xz (剪切)
      rtvec[15] = omega;     // yz (剪切)

      // 幽灵矩 / 高阶非流体力学矩 (无宏观物理意义，强制拉满以耗散 3D 噪声)
      // 对应索引: 10, 12, 16, 17, 18
      rtvec[10] = s_ghost;   // Ghost moment (MAXIMUM DISSIPATION)
      rtvec[12] = s_ghost;   // Ghost moment
      rtvec[16] = s_ghost;   // Ghost moment
      rtvec[17] = s_ghost;   // Ghost moment
      rtvec[18] = s_ghost;   // Ghost moment
      // =====================================================================

    } else if constexpr (std::is_same_v<LatSet, D2Q9<T>>) {
      // 2D 保持你原有的逻辑（2D 矩空间紧凑，原有逻辑已足够稳定）
      T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
      rtvec[0] = T{0}; rtvec[1] = omega; rtvec[2] = omega;
      rtvec[3] = T{0}; rtvec[4] = s_q;
      rtvec[5] = T{0}; rtvec[6] = s_q;
      rtvec[7] = omega; rtvec[8] = omega;
    } else {
      for (unsigned int i = 0; i < LatSet::q; ++i) rtvec[i] = omega;
    }

    // 4. Moments from population
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];

    // 5. Equilibrium in moment space
    std::array<T, LatSet::q> feq{};
    if constexpr (isIncompressible) {
      // Incompressible (He-Luo): feq = w*(p + 3*cu + 4.5*cu² - 1.5*u²), p = Σf
      equilibrium::IncompressibleSecondOrder<CELL>::apply(feq, zeroth, uc);
    } else {
      // Compressible: feq = w*rho*(1 + 3*cu + 4.5*cu² - 1.5*u²)
      equilibrium::SecondOrder<CELL>::apply(feq, rho, uc);
    }
    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];

    // 6. Force source in population space (Guo force scheme)
    T uF = T{};
    for (unsigned int d = 0; d < LatSet::d; ++d) uF += uc[d] * F[d];
    T fi_force[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cF = T{}, cu = T{};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cF += latset::c<LatSet>(i)[d] * F[d];
        cu += latset::c<LatSet>(i)[d] * uc[d];
      }
      // Incompressible Guo force needs /rho factor (Fortran L824: DF = w*[...] * dt/rho, dt=1)
      // Compressible Guo force has no /rho factor
      if constexpr (isIncompressible) {
        fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                       + LatSet::InvCs4 * cu * cF) / rho;
      } else {
        fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                       + LatSet::InvCs4 * cu * cF);
      }
    }

    // 7. Force source in moment space
    T force_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j)
      for (unsigned int k = 0; k < LatSet::q; ++k)
        force_m[j] += mrt::M<LatSet>(j, k) * fi_force[k];

    // 8. Combined MRT collision + force source
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll += invM_ij * rtvec[j] * (momenta[j] - momentaEq[j]);
        source += invM_ij * (T{1} - rtvec[j] / T{2}) * force_m[j];
      }
      cell[i] = cell[i] - coll + source;
    }
  }
};

// ===================================================================
// MRTMag — MRT collision for magnetic potential equation (D2Q5).
// Hu & Li (2018) Phys. Rev. E 98, 033301, §II.D Eq.(40)-(44)
//
// Solves: ∂ψ/∂t = ∇·(εμ∇ψ) via pseudo-time relaxation.
// Relaxation: S^h = diag(1, ω_mag, ω_mag, 1, 1)
//   ψ (m₀):   conserved — s=1, m₀=m_eq₀=ψ identically
//   Jx,Jy:    relaxed by ω_mag = 1/τ_mag (diffusion rate)
//   e, pxx:   fully relaxed to equilibrium=0 with s=1
//
// Equilibrium: h_α^eq = w_eq·ψ  (uniform weights = 1/q)
// Paper default: τ_mag = 4.0 → ω_mag = 0.25  (§III)
// ===================================================================
template <typename CELL, bool WriteToField = false>
struct MRTMag {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = equilibrium::MagEquilibrium<CELL>;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // Per-cell omega (spatially varying via Eq.43: τ = 0.5 + εμ/(c_s²·Δt))
    // Falls back to block-level converter omega if no per-cell field
    T omega_mag;
    if constexpr (cell.template hasField<GenericField<GenericArray<T>, ff::MagOmegaBase>>()) {
      omega_mag = cell.template get<GenericField<GenericArray<T>, ff::MagOmegaBase>>();
    } else {
      omega_mag = cell.getOmega();
    }

    // Relaxation vector: S = (1, ω_mag, ω_mag, 1, 1)
    T rtvec[LatSet::q];
    rtvec[0] = T{1};
    rtvec[1] = omega_mag;
    rtvec[2] = omega_mag;
    for (unsigned int i = 3; i < LatSet::q; ++i) rtvec[i] = T{1};

    // InvM_S = InvM · diag(rtvec)
    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    // Forward transform: m = M · f
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }

    // Collision: f'[i] = f[i] - Σ_{j>0} InvM[i,j]·s_j·m[j]
    // m_eq[j] = 0 ∀ j>0 (uniform equilibrium → vanishing higher moments)
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      for (unsigned int j = 1; j < LatSet::q; ++j) {
        coll += InvM_S[i][j] * momenta[j];
      }
      cell[i] -= coll;
    }
  }
};

}  // namespace collision