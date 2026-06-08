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

    // relaxation-time vector: UseCHRelaxation=false → interface sharpening (s=1 for j-moments)
    //                        UseCHRelaxation=true  → Cahn-Hilliard     (omega_phi for j-moments)
    const T rtvec[LatSet::q] {
      (UseCHRelaxation ? omega_phi : T{1}),          // 0: rho (phi)
      T{14./10.},                                     // 1: e
      T{14./10.},                                     // 2: epsilon
      (UseCHRelaxation ? omega_phi : T{1}),          // 3: jx
      T{12./10.},                                     // 4: qx
      (UseCHRelaxation ? omega_phi : T{1}),          // 5: jy
      T{12./10.},                                     // 6: qy
      omega_phi,                                      // 7: pxx (shear)
      omega_phi                                       // 8: pxy (shear)
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
// replaces BGKForce<moment::forcerhoU, equilibrium::SecondOrder, force::Force>
template <typename T, typename TypePack, typename ForceField>
struct MRTForce<CELL<T, D2Q9<T>, TypePack>, ForceField> {
  using LatSet = D2Q9<T>;

  __any__ static void apply(CELL<T, D2Q9<T>, TypePack>& cell) {
    const T omega = cell.getOmega();
    const auto& F = cell.template get<ForceField>();

    // 1. Compute rho and u from populations (unrolled D2Q9)
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
          + cell[5] + cell[6] + cell[7] + cell[8];
    T ux = (cell[1] - cell[2] + cell[5] - cell[6] + cell[7] - cell[8]) / rho;
    T uy = (cell[3] - cell[4] + cell[5] - cell[6] - cell[7] + cell[8]) / rho;

    // 2. Half-force correction: uc = u + F/(2*rho)
    T halfInvRho = T{0.5} / rho;
    T ucx = ux + halfInvRho * F[0];
    T ucy = uy + halfInvRho * F[1];

    // Write macroscopic fields
    cell.template get<RHO<T>>() = rho;
    cell.template get<VELOCITY<T, LatSet::d>>()[0] = ucx;
    cell.template get<VELOCITY<T, LatSet::d>>()[1] = ucy;

    // 3. MRT relaxation time vector
    //    Conserved moments: 0,3,5 → s=0
    //    Shear moments: 7,8 → s=omega (determines fluid viscosity)
    //    High-order moments: 1,2,4,6 → s=s_e (tunable for stability)
    const T rtvec[LatSet::q] {
      T{0},             // rho
      T{11./10.},       // e
      T{11./10.},       // epsilon
      T{0},             // jx
      T{11./10.},       // qx
      T{0},             // jy
      T{11./10.},       // qy
      omega,            // pxx
      omega             // pxy
    };

    // 4. Moments from populations (M·f, unrolled D2Q9)
    T m[LatSet::q];
    m[0] = rho;
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

    // 5. Equilibrium moments (analytical, using corrected velocity uc)
    //    Derived from M·feq for SecondOrder D2Q9 equilibrium
    T ucx2 = ucx * ucx;
    T ucy2 = ucy * ucy;
    T uc2 = ucx2 + ucy2;
    T meq[LatSet::q];
    meq[0] = rho;
    meq[1] = T{-2} * rho + T{3} * rho * uc2;
    meq[2] = T{9} * rho * ucx2 * ucy2 - T{3} * rho * uc2 + rho;
    meq[3] = rho * ucx;
    meq[4] = rho * ucx * (T{3} * ucy2 - T{1});
    meq[5] = rho * ucy;
    meq[6] = rho * ucy * (T{3} * ucx2 - T{1});
    meq[7] = rho * (ucx2 - ucy2);
    meq[8] = rho * ucx * ucy;

    // 6. Deviation in moment space
    T dm[LatSet::q];
    dm[0] = m[0] - meq[0];
    dm[1] = m[1] - meq[1];
    dm[2] = m[2] - meq[2];
    dm[3] = m[3] - meq[3];
    dm[4] = m[4] - meq[4];
    dm[5] = m[5] - meq[5];
    dm[6] = m[6] - meq[6];
    dm[7] = m[7] - meq[7];
    dm[8] = m[8] - meq[8];

    // 7. Guo force in population space (unrolled D2Q9, using uc)
    //    S_i = w_i * [InvCs2*(c_i·F - u·F) + InvCs4*(c_i·u)*(c_i·F)]
    T Fx = F[0], Fy = F[1];
    T uF = ucx * Fx + ucy * Fy;
    T S[LatSet::q];
    // i=0: c=(0,0), w=4/9
    S[0] = T{4./9.} * T{3} * (-uF);
    // i=1: c=(1,0), w=1/9
    S[1] = T{1./9.} * (T{3} * (Fx - uF) + T{9} * ucx * Fx);
    // i=2: c=(-1,0), w=1/9
    S[2] = T{1./9.} * (T{3} * (-Fx - uF) + T{9} * ucx * Fx);
    // i=3: c=(0,1), w=1/9
    S[3] = T{1./9.} * (T{3} * (Fy - uF) + T{9} * ucy * Fy);
    // i=4: c=(0,-1), w=1/9
    S[4] = T{1./9.} * (T{3} * (-Fy - uF) + T{9} * ucy * Fy);
    // i=5: c=(1,1), w=1/36
    S[5] = T{1./36.} * (T{3} * (Fx + Fy - uF) + T{9} * (ucx + ucy) * (Fx + Fy));
    // i=6: c=(-1,-1), w=1/36
    S[6] = T{1./36.} * (T{3} * (-Fx - Fy - uF) + T{9} * (ucx + ucy) * (Fx + Fy));
    // i=7: c=(1,-1), w=1/36
    S[7] = T{1./36.} * (T{3} * (Fx - Fy - uF) + T{9} * (ucx - ucy) * (Fx - Fy));
    // i=8: c=(-1,1), w=1/36
    S[8] = T{1./36.} * (T{3} * (-Fx + Fy - uF) + T{9} * (ucx - ucy) * (Fx - Fy));

    // 8. Force in moment space (F_m = M·S, unrolled D2Q9)
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

    // 9. Combined MRT collision + force source
    //    f_new = f - InvM·S·(m - m_eq) + InvM·(I - S/2)·F_m
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