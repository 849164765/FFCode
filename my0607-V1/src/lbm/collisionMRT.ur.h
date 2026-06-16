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

template <typename T, typename TypePack, bool WriteToField, bool UseCHRelaxation>
struct MRTSource<equilibrium::SecondOrder<CELL<T, D2Q9<T>, TypePack>>, NORMAL<T, 2>, WriteToField, UseCHRelaxation> {
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

    // relaxation-time vector — Fortran-aligned for UseCHRelaxation=true
    // Fortran Sf = {1.0, 1.1, 1.1, 1/tau, 1/tau, 1/tau, 1/tau, 1.2, 1.2}
    const T rtvec[LatSet::q] {
      T{1},                                              // 0: rho/phi (Fortran: s=1.0)
      (UseCHRelaxation ? T{11./10.} : T{14./10.}),       // 1: e
      (UseCHRelaxation ? T{11./10.} : T{14./10.}),       // 2: epsilon
      (UseCHRelaxation ? omega_phi  : T{1}),             // 3: jx (Fortran: 1/tau for CH)
      (UseCHRelaxation ? omega_phi  : T{12./10.}),       // 4: qx (Fortran: 1/tau for CH)
      (UseCHRelaxation ? omega_phi  : T{1}),             // 5: jy
      (UseCHRelaxation ? omega_phi  : T{12./10.}),       // 6: qy
      (UseCHRelaxation ? T{12./10.} : omega_phi),        // 7: pxx (Fortran: 1.2 for CH)
      (UseCHRelaxation ? T{12./10.} : omega_phi)         // 8: pxy
    };

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

// MRT collision with Guo force scheme for Navier-Stokes (D2Q9 specialization)
// Fortran-accurate MRT: proper relaxation vector with s_q for q-moments
template <typename T, typename TypePack, typename ForceField>
struct MRTForce<CELL<T, D2Q9<T>, TypePack>, ForceField> {
  using LatSet = D2Q9<T>;

  __any__ static void apply(CELL<T, D2Q9<T>, TypePack>& cell) {
    // Per-cell omega (spatially varying relaxation) for variable viscosity
    T omega{};
    if constexpr (cell.template hasField<OMEGA<T>>()) {
      omega = cell.template get<OMEGA<T>>();
    } else {
      omega = cell.getOmega();
    }
    const auto& F = cell.template get<ForceField>();

    // Fortran-accurate MRT relaxation vector
    // S = diag(0, ω, ω, 0, s_q, 0, s_q, ω, ω)
    T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
    const T rtvec[LatSet::q] {
      T{0}, omega, omega, T{0}, s_q, T{0}, s_q, omega, omega
    };

    // Compile-time check: incompressible (He-Luo) vs compressible LBM
    static constexpr bool isIncompressible =
        CELL<T, D2Q9<T>, TypePack>::template hasField<DENSITY<T>>();

    // ---- Common: zeroth moment (p for incomp, rho for comp) and raw momentum ----
    T zeroth = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
             + cell[5] + cell[6] + cell[7] + cell[8];
    T ux_raw = cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8];
    T uy_raw = cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8];

    // ---- Branched: density source and velocity half-force correction ----
    T dens{}, ucx{}, ucy{};
    if constexpr (isIncompressible) {
      dens = cell.template get<DENSITY<T>>();   // from phi interpolation
      ucx = (ux_raw + T{0.5} * F[0]) / dens;
      ucy = (uy_raw + T{0.5} * F[1]) / dens;
      cell.template get<PRESSURE<T>>() = zeroth;
    } else {
      dens = zeroth;                              // sum(f) ~ 1.0
      ucx = (ux_raw + T{0.5} * F[0]) / dens;
      ucy = (uy_raw + T{0.5} * F[1]) / dens;
      cell.template get<RHO<T>>() = zeroth;
    }
    cell.template get<VELOCITY<T, LatSet::d>>()[0] = ucx;
    cell.template get<VELOCITY<T, LatSet::d>>()[1] = ucy;

    // ---- Common: raw moments from populations (M·f, unrolled D2Q9) ----
    T m[LatSet::q];
    m[0] = zeroth;
    m[1] = T{-4} * cell[0] - cell[1] - cell[2] - cell[3] - cell[4]
         + T{2} * (cell[5] + cell[6] + cell[7] + cell[8]);
    m[2] = T{4} * cell[0] - T{2} * (cell[1] + cell[2] + cell[3] + cell[4])
         + cell[5] + cell[6] + cell[7] + cell[8];
    m[3] = cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8];
    m[4] = T{-2} * cell[1] + T{2} * cell[2] + cell[5] - cell[6] + cell[7] - cell[8];
    m[5] = cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8];
    m[6] = T{-2} * cell[3] + T{2} * cell[4] + cell[5] - cell[6] - cell[7] + cell[8];
    m[7] = cell[1] + cell[2] - cell[3] - cell[4];
    m[8] = cell[5] + cell[6] - cell[7] - cell[8];

    // ---- Equilibrium moments (branched) ----
    T ucx2 = ucx * ucx;
    T ucy2 = ucy * ucy;
    T uc2 = ucx2 + ucy2;
    T meq[LatSet::q];
    if constexpr (isIncompressible) {
      // He-Luo incompressible: velocity moments have NO rho/p factor
      meq[0] = zeroth;
      meq[1] = T{-2} * zeroth + T{3} * uc2;
      meq[2] = zeroth + T{9} * ucx2 * ucy2 - T{3} * uc2;
      meq[3] = ucx;
      meq[4] = ucx * (T{3} * ucy2 - T{1});
      meq[5] = ucy;
      meq[6] = ucy * (T{3} * ucx2 - T{1});
      meq[7] = ucx2 - ucy2;
      meq[8] = ucx * ucy;
    } else {
      // Compressible: velocity moments multiplied by rho = zeroth
      meq[0] = zeroth;
      meq[1] = T{-2} * zeroth + T{3} * zeroth * uc2;
      meq[2] = T{9} * zeroth * ucx2 * ucy2 - T{3} * zeroth * uc2 + zeroth;
      meq[3] = zeroth * ucx;
      meq[4] = zeroth * ucx * (T{3} * ucy2 - T{1});
      meq[5] = zeroth * ucy;
      meq[6] = zeroth * ucy * (T{3} * ucx2 - T{1});
      meq[7] = zeroth * (ucx2 - ucy2);
      meq[8] = zeroth * ucx * ucy;
    }

    // ---- Common: deviation in moment space ----
    T dm[LatSet::q];
    for (unsigned int i = 0; i < LatSet::q; ++i) dm[i] = m[i] - meq[i];

    // ---- Guo force source: incompressible divides by actual density ----
    T Fx = F[0], Fy = F[1];
    T uF = ucx * Fx + ucy * Fy;
    T invDens = T{1} / dens;
    T S[LatSet::q];
    S[0] = T{4./9.} * T{3} * (-uF) * invDens;
    S[1] = T{1./9.} * (T{3} * (Fx - uF) + T{9} * ucx * Fx) * invDens;
    S[2] = T{1./9.} * (T{3} * (-Fx - uF) + T{9} * ucx * Fx) * invDens;
    S[3] = T{1./9.} * (T{3} * (Fy - uF) + T{9} * ucy * Fy) * invDens;
    S[4] = T{1./9.} * (T{3} * (-Fy - uF) + T{9} * ucy * Fy) * invDens;
    S[5] = T{1./36.} * (T{3} * (Fx + Fy - uF) + T{9} * (ucx + ucy) * (Fx + Fy)) * invDens;
    S[6] = T{1./36.} * (T{3} * (-Fx - Fy - uF) + T{9} * (ucx + ucy) * (Fx + Fy)) * invDens;
    S[7] = T{1./36.} * (T{3} * (Fx - Fy - uF) + T{9} * (ucx - ucy) * (Fx - Fy)) * invDens;
    S[8] = T{1./36.} * (T{3} * (-Fx + Fy - uF) + T{9} * (ucx - ucy) * (Fx - Fy)) * invDens;

    // ---- Common: force in moment space (F_m = M·S, unrolled D2Q9) ----
    T Fm[LatSet::q];
    Fm[0] = S[0] + S[1] + S[2] + S[3] + S[4] + S[5] + S[6] + S[7] + S[8];
    Fm[1] = T{-4} * S[0] - S[1] - S[2] - S[3] - S[4]
          + T{2} * (S[5] + S[6] + S[7] + S[8]);
    Fm[2] = T{4} * S[0] - T{2} * (S[1] + S[2] + S[3] + S[4])
          + S[5] + S[6] + S[7] + S[8];
    Fm[3] = S[1] - S[2] + S[5] - S[6] + S[7] - S[8];
    Fm[4] = T{-2} * S[1] + T{2} * S[2] + S[5] - S[6] + S[7] - S[8];
    Fm[5] = S[3] - S[4] + S[5] - S[6] - S[7] + S[8];
    Fm[6] = T{-2} * S[3] + T{2} * S[4] + S[5] - S[6] - S[7] + S[8];
    Fm[7] = S[1] + S[2] - S[3] - S[4];
    Fm[8] = S[5] + S[6] - S[7] - S[8];

    // ---- Common: MRT collision + force source ----
    // f_new = f - InvM·S·(m - m_eq) + InvM·(I - S/2)·F_m
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll   += invM_ij * rtvec[j] * dm[j];
        source += invM_ij * (T{1} - rtvec[j] / T{2}) * Fm[j];
      }
      cell[i] = cell[i] - coll + source;
    }
  }
};

}

#endif