// collisionMRT.h

#pragma once

#include "lbm/moment.ur.h"

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
      // Cahn-Hilliard: use omega_phi for all → numerically matches BGKSource
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        rtvec[i] = omega_phi;
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

  __any__ static void apply(CELL& cell) {
    const T omega = cell.getOmega();
    const auto& F = cell.template get<ForceField>();

    // 1. Compute rho and u_raw from populations
    T rho{};
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      rho += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
      }
    }
    for (unsigned int d = 0; d < LatSet::d; ++d) u[d] /= rho;

    // 2. Half-force correction: u_c = u + F/(2*rho)
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d) uc[d] += F[d] / (T{2} * rho);

    // Write macroscopic fields
    cell.template get<RHO<T>>() = rho;
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;

    // 3. D2Q9 MRT relaxation time vector
    //    全设为 omega，使 MRTForce 数值等价于 BGKForce
    //    稳定后可按需调大指定矩的 s_k 来利用 MRT 优势
    T rtvec[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) rtvec[i] = omega;

    // 4. Moments from population
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];

    // 5. Equilibrium in moment space
    std::array<T, LatSet::q> feq{};
    equilibrium::SecondOrder<CELL>::apply(feq, rho, uc);
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
      fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                     + LatSet::InvCs4 * cu * cF);
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

// MRT collision for velocity-based Navier-Stokes per Guo et al. 2025.
//
// Eq.29: g̃ = g - M⁻¹·S·M·(g-g_eq) + M⁻¹·(I-S/2)·M·G
// Eq.33: 1/s₇ = 1/s₈ = η/(ρc²) + 0.5  (per-cell shear omega)
// Conserved: s₀, s₃, s₅ = 0.  Non-conserved: standard D2Q9 stability values.
//
// Requires: OMEGA<T> field for per-cell ω(φ), RHO<T> for ρ(φ),
//           PRESSURE<T>, FORCE<T,D>, VELOCITY<T,D>.
template <typename CELL, typename ForceField>
struct MRTVelocityBased {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    const auto& F = cell.template get<ForceField>();

    // ---- 1. Velocity-based macro update (Eqs.35-36) ----
    T gsum = T{};
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      gsum += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d)
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
    }
    const T rho_phi = cell.template get<typename CELL::GenericRho>(); // ρ(φ)
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d)
      uc[d] += F[d] / (T{2} * rho_phi);          // Eq.35: u = Σe·g + F/(2ρ)
    T pressure = rho_phi * LatSet::cs2 * gsum;     // Eq.36: p = ρc²·Σg
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;
    cell.template get<PRESSURE<T>>() = pressure;

    // ---- 2. Relaxation rate vector (standard D2Q9 MRT + uniform omega) ----
    // Using block-level omega for shear to avoid force imbalance from
    // per-cell ω variation: (1-s₇/2) must be uniform for consistent F/ρ.
    // Per-cell ω(φ) (Eq.33) can be enabled after force-term renormalization.
    T omega_sh = cell.getOmega();

    T rtvec[LatSet::q] {};
    rtvec[0] = T{0};                              // s₀: conserved (density)
    rtvec[1] = T{1.2};                            // s₁: energy
    rtvec[2] = T{1.4};                            // s₂: energy square
    rtvec[3] = T{0};                              // s₃: conserved (momentum x)
    rtvec[4] = T{1.2};                            // s₄: energy flux x
    rtvec[5] = T{0};                              // s₅: conserved (momentum y)
    rtvec[6] = T{1.2};                            // s₆: energy flux y
    rtvec[7] = omega_sh;                          // s₇: shear (xx stress)
    rtvec[8] = omega_sh;                          // s₈: shear (xy stress)

    // ---- 3. Population moments m = M·g ----
    T m[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        m[i] += mrt::M<LatSet>(i, j) * cell[j];

    // ---- 4. Equilibrium moments m_eq = M·g_eq (Eq.31) ----
    std::array<T, LatSet::q> feq{};
    equilibrium::NSFieldEquilibrium<CELL>::apply(feq, rho_phi, pressure, uc);
    T mEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        mEq[i] += mrt::M<LatSet>(i, j) * feq[j];

    // ---- 5. Guo force in population space (Eq.32), divide by ρ (Eq.25) ----
    T uF = uc[0] * F[0] + uc[1] * F[1];
    T invRho = T{1} / rho_phi;
    T G_pop[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cF = T{}, cu = T{};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cF += latset::c<LatSet>(i)[d] * F[d];
        cu += latset::c<LatSet>(i)[d] * uc[d];
      }
      G_pop[i] = latset::w<LatSet>(i)
                * (LatSet::InvCs2 * (cF - uF) + LatSet::InvCs4 * cu * cF)
                * invRho;
    }
    // ---- 6. Force in moment space: F_m = M·G_pop ----
    T F_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j)
      for (unsigned int k = 0; k < LatSet::q; ++k)
        F_m[j] += mrt::M<LatSet>(j, k) * G_pop[k];

    // ---- 7. MRT collision + force: Eq.29 ----
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll = T{}, src = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll += invM_ij * rtvec[j] * (m[j] - mEq[j]);
        src  += invM_ij * (T{1} - rtvec[j] / T{2}) * F_m[j];
      }
      cell[i] = cell[i] - coll + src;
    }

    // ---- 8. Light pressure relaxation (zeroth moment only) ----
    // s₀=0 conserves Σg, but streaming across interface mixes 1/ρ from
    // neighboring cells with different φ, creating Σg ≠ 1/ρ(φ_local).
    // MRT separation allows this without disturbing shear modes.
    {
      T gsum2 = T{};
      for (unsigned int i = 0; i < LatSet::q; ++i) gsum2 += cell[i];
      T g_target = T{1} / rho_phi;
      T delta = gsum2 - g_target;
      constexpr T omega_p = T{0.1};
      for (unsigned int i = 0; i < LatSet::q; ++i)
        cell[i] -= omega_p * latset::w<LatSet>(i) * delta;
      cell.template get<PRESSURE<T>>() = rho_phi * LatSet::cs2 * (gsum2 - omega_p * delta);
    }
  }
};

// ===== Simplified variant: all s = omega (BGK-equivalent MRT) =====
template <typename CELL, typename ForceField>
struct MRTVelocityBasedBGK {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;

  __any__ static void apply(CELL& cell) {
    T omega;
    if constexpr (CELL::template hasField<OMEGA<T>>()) {
      omega = cell.template get<OMEGA<T>>();
    } else {
      omega = cell.getOmega();
    }
    const auto& F = cell.template get<ForceField>();

    // Macro update (velocity-based)
    T gsum = T{};
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      gsum += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d)
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
    }
    const T rho_phi = cell.template get<typename CELL::GenericRho>();
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d) uc[d] += F[d] / (T{2} * rho_phi);
    T pressure = rho_phi * LatSet::cs2 * gsum;
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;
    cell.template get<PRESSURE<T>>() = pressure;

    // MRT with all s = omega (equivalent to BGK in moment space)
    T rtvec[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) rtvec[i] = omega;

    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];

    std::array<T, LatSet::q> feq{};
    equilibrium::NSFieldEquilibrium<CELL>::apply(feq, rho_phi, pressure, uc);
    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];

    T invRho_bgk = T{1} / rho_phi;
    T uF_dot = uc[0] * F[0] + uc[1] * F[1];
    T fi_force[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cF = T{}, cu = T{};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cF += latset::c<LatSet>(i)[d] * F[d];
        cu += latset::c<LatSet>(i)[d] * uc[d];
      }
      fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF_dot)
                     + LatSet::InvCs4 * cu * cF) * invRho_bgk;
    }
    T force_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j)
      for (unsigned int k = 0; k < LatSet::q; ++k)
        force_m[j] += mrt::M<LatSet>(j, k) * fi_force[k];

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{}, source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll += invM_ij * rtvec[j] * (momenta[j] - momentaEq[j]);
        source += invM_ij * (T{1} - rtvec[j] / T{2}) * force_m[j];
      }
      cell[i] = cell[i] - coll + source;
    }

    // Pressure relaxation
    {
      T gsum2{};
      for (unsigned int i = 0; i < LatSet::q; ++i) gsum2 += cell[i];
      T gsum_target = T{1} / rho_phi;
      T delta = gsum2 - gsum_target;
      constexpr T omega_p = T{0.05};
      for (unsigned int i = 0; i < LatSet::q; ++i)
        cell[i] -= omega_p * latset::w<LatSet>(i) * delta;
      cell.template get<PRESSURE<T>>() = rho_phi * LatSet::cs2 * (gsum2 - omega_p * delta);
    }
  }
};

}