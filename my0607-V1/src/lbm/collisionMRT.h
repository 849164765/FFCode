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

template <typename EquilibriumScheme, typename NORMAL, bool WriteToField = false>
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
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rtvec[i] = mrt::s1<LatSet>(i);
    }
    for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i) {
      rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega_phi;
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

}