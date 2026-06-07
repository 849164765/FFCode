// collisionMRT.h

#pragma once


#include "lbm/moment.ur.h"
#include "lbm/collisionMRT.h"

#ifdef _UNROLLFOR

namespace collision {

#ifdef __CUDA_ARCH__
template <typename T, typename LatSet, typename TypePack>
using CELL = cudev::Cell<T, LatSet, TypePack>;

#else
template <typename T, typename LatSet, typename TypePack>
using CELL = Cell<T, LatSet, TypePack>;
#endif

// a typical MRT collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
template <typename T, typename TypePack, bool WriteToField>
struct MRT_Feq_RhoU<CELL<T, D2Q9<T> ,TypePack>, WriteToField> {
  using LatSet = D2Q9<T>;
  using equilibriumscheme = equilibrium::SecondOrder<CELL<T, D2Q9<T> ,TypePack>>;
  using GenericRho = typename CELL<T, D2Q9<T> ,TypePack>::GenericRho;

  __any__ static void apply(CELL<T, D2Q9<T> ,TypePack>& cell) {
    const T omega = cell.getOmega();
    // relaxation time vector
    const T rtvec[LatSet::q] {T{}, T{11./10.}, T{11./10.}, T{}, T{11./10.}, T{}, T{11./10.}, omega, omega};
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
    moment::template rhoU<CELL<T, D2Q9<T> ,TypePack>, WriteToField>::apply(cell, rho, u);

    // compute Momenta
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }
    // compute Equilibrium
    // std::array<T, LatSet::q> feq{};
    // equilibrium::SecondOrder<CELL<T, D2Q9<T> ,TypePack>>::apply(feq, rho, u);
    // T momentaEq[LatSet::q] {};
    // for(unsigned int i = 0; i < LatSet::q; ++i) {
    //   momentaEq[i] = T{};
    //   for (unsigned int j = 0; j < LatSet::q; ++j) {
    //     momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
    //   }
    // }
    // T momentaEq[LatSet::q] {};
    // const T ux2 = u[0] * u[0];
    // const T uy2 = u[1] * u[1];
    // const T rhoux = rho * u[0];
    // const T rhouy = rho * u[1];
    // const T cse0 = 3 * rho * (ux2 + uy2);
    // momentaEq[0] = rho;
    // momentaEq[1] = -2 * rho + cse0;
    // momentaEq[2] = 9 * rho * ux2 * uy2 - cse0 + rho;
    // momentaEq[3] = rhoux;
    // momentaEq[4] = rhoux * (3 * uy2 - 1);
    // momentaEq[5] = rhouy;
    // momentaEq[6] = rhouy * (3 * ux2 - 1);
    // momentaEq[7] = rho * (ux2 - uy2);
    // momentaEq[8] = rho * u[0] * u[1];

    // compute off-Equilibrium part: momenta - momentaEq
    T delmomenta[LatSet::q] {};
    const T ux2 = u[0] * u[0];
    const T uy2 = u[1] * u[1];
    const T rhoux = rho * u[0];
    const T rhouy = rho * u[1];
    const T cse0 = 3 * rho * (ux2 + uy2);
    delmomenta[0] = momenta[0] - rho;
    delmomenta[1] = momenta[1] + 2 * rho - cse0;
    delmomenta[2] = momenta[2] - 9 * rho * ux2 * uy2 + cse0 - rho;
    delmomenta[3] = momenta[3] - rhoux;
    delmomenta[4] = momenta[4] - rhoux * (3 * uy2 - 1);
    delmomenta[5] = momenta[5] - rhouy;
    delmomenta[6] = momenta[6] - rhouy * (3 * ux2 - 1);
    delmomenta[7] = momenta[7] - rho * (ux2 - uy2);
    delmomenta[8] = momenta[8] - rho * u[0] * u[1];

    // MRT collision
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * delmomenta[j];
      }
      cell[i] -= momentaEqS;
    }
  }

};

template <typename T, typename TypePack, bool WriteToField>
struct MRTSource<equilibrium::SecondOrder<CELL<T, D2Q9<T>, TypePack>>, NORMAL<T, 2>, WriteToField> {
  using LatSet = D2Q9<T>;
  using CELLTYPE = CELL<T, D2Q9<T>, TypePack>;
  using equilibriumscheme = equilibrium::SecondOrder<CELLTYPE>;
  using GenericRho = typename CELLTYPE::GenericRho;

  __any__ static void apply(CELLTYPE& cell) {
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL<T, 2>>();
    T phi = cell.template get<GenericRho>();

    T omega_phi = cell.getOmega();
    T interfacewidth = cell.template get<INTERFACEWIDTH<T>>();

    const T rtvec[LatSet::q] {T{1}, T{14./10.}, T{14./10.}, T{1}, T{12./10.}, T{1}, T{12./10.}, omega_phi, omega_phi};

    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    T momenta[LatSet::q] {};
    momenta[0] = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
    momenta[1] = T{-4} * cell[0] - cell[1] - cell[2] - cell[3] - cell[4] + T{2} * (cell[5] + cell[6] + cell[7] + cell[8]);
    momenta[2] = T{4} * cell[0] - T{2} * (cell[1] + cell[2] + cell[3] + cell[4]) + cell[5] + cell[6] + cell[7] + cell[8];
    momenta[3] = cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8];
    momenta[4] = T{-2} * cell[1] + T{2} * cell[2] + cell[5] - cell[6] + cell[7] - cell[8];
    momenta[5] = cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8];
    momenta[6] = T{-2} * cell[3] + T{2} * cell[4] + cell[5] - cell[6] - cell[7] + cell[8];
    momenta[7] = cell[1] + cell[2] - cell[3] - cell[4];
    momenta[8] = cell[5] + cell[6] - cell[7] - cell[8];

    std::array<T, LatSet::q> feq{};
    equilibriumscheme::apply(feq, phi, u);

    T momentaEq[LatSet::q] {};
    momentaEq[0] = feq[0] + feq[1] + feq[2] + feq[3] + feq[4] + feq[5] + feq[6] + feq[7] + feq[8];
    momentaEq[1] = T{-4} * feq[0] - feq[1] - feq[2] - feq[3] - feq[4] + T{2} * (feq[5] + feq[6] + feq[7] + feq[8]);
    momentaEq[2] = T{4} * feq[0] - T{2} * (feq[1] + feq[2] + feq[3] + feq[4]) + feq[5] + feq[6] + feq[7] + feq[8];
    momentaEq[3] = feq[1] - feq[2] + feq[5] - feq[6] + feq[7] - feq[8];
    momentaEq[4] = T{-2} * feq[1] + T{2} * feq[2] + feq[5] - feq[6] + feq[7] - feq[8];
    momentaEq[5] = feq[3] - feq[4] + feq[5] - feq[6] - feq[7] + feq[8];
    momentaEq[6] = T{-2} * feq[3] + T{2} * feq[4] + feq[5] - feq[6] - feq[7] + feq[8];
    momentaEq[7] = feq[1] + feq[2] - feq[3] - feq[4];
    momentaEq[8] = feq[5] + feq[6] - feq[7] - feq[8];

    T delmomenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      delmomenta[i] = momenta[i] - momentaEq[i];
    }

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * delmomenta[j];
      }
      cell[i] -= momentaEqS;
    }

    T A = T{4} * phi * (T{1} - phi) / interfacewidth;

    T S_body[LatSet::q] {};
    S_body[0] = T{};
    S_body[1] = T{1} / T{9} * n[0] * A;
    S_body[2] = T{1} / T{9} * (-n[0]) * A;
    S_body[3] = T{1} / T{9} * n[1] * A;
    S_body[4] = T{1} / T{9} * (-n[1]) * A;
    S_body[5] = T{1} / T{36} * (n[0] + n[1]) * A;
    S_body[6] = T{1} / T{36} * (-n[0] - n[1]) * A;
    S_body[7] = T{1} / T{36} * (n[0] - n[1]) * A;
    S_body[8] = T{1} / T{36} * (-n[0] + n[1]) * A;

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

}

#endif