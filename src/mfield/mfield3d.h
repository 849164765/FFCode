#pragma once

#include "mfield/mfield2d.h"

namespace mfield {

// ===================================================================
//  3D 磁标势场扩展 (3D magnetic field extension)
//
//  复用 mfield2d.h 中与维度无关的标量场定义:
//    PSI, OMEGA_PSI, MU_PERCELL, CHI_PERCELL, HX, HY, HMAG
//    MU_L, MU_H, CHI_L, CHI_H, H_0, PSI_K
//  仅新增 z 分量 HZ 与 3D 专属的场包/算子。
// ===================================================================

// --- Per-cell varying field ---
struct HzBase : public FieldBase<1> {};  // H_z = -∂ψ/∂z

template <typename T>
using HZ = GenericField<GenericArray<T>, HzBase>;

// ===================================================================
//  3D Field Pack
// ===================================================================

template <typename T>
using MFFIELDS3D = TypePack<
    PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>,
    HX<T>, HY<T>, HZ<T>, HMAG<T>,
    MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>, H_0<T>, PSI_K<T>>;

// ===================================================================
//  Functor Declarations (3D)
// ===================================================================

// MF1-3D: PF → MF coupling — update per-cell μ(φ), χ(φ), ω_ψ from φ
template <typename PFCELL, typename MFCELL>
struct MFUpdateCoeffs3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell);
};

// MF2-3D: H = -∇ψ via D3Q7 6-direction isotropic gradient
// NOTE: D3Q7 的权重二阶矩 Σ w_k c_ka² = 1/4 ≠ cs² = 1/3，归一化必须除以
// Q2 = Σ w_k c_ka²（D2Q5 恰好 Q2 = 1/3 = cs²，因此 2D 直接除以 cs²）。
template <typename CELL>
struct MFComputeH3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// MF3-3D: Magnetic body force F_mag = +χ·|H|·∇|H|  (Kelvin force, PF+MF → NS)
template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

// ===================================================================
//  Communication Functions (3D)
// ===================================================================

template <typename T, typename LATTICE>
void CommunicateHZ(LATTICE& lattice) {
  lattice.template getField<HZ<T>>().Communicate();
}

// Batch communicate all per-cell magnetic fields that need ghost sync
template <typename T, typename LATTICE>
void CommunicateAllMFFields3D(LATTICE& lattice) {
  lattice.template getField<PSI<T>>().Communicate();
  lattice.template getField<HX<T>>().Communicate();
  lattice.template getField<HY<T>>().Communicate();
  lattice.template getField<HZ<T>>().Communicate();
  lattice.template getField<HMAG<T>>().Communicate();
  lattice.template getField<OMEGA_PSI<T>>().Communicate();
  lattice.template getField<CHI_PERCELL<T>>().Communicate();
}

// Broadcast block-level params from rank 0
template <typename T, typename LATTICE>
void BroadcastAllMFParams3D(LATTICE& lattice,
                            T& mu_l, T& mu_h, T& chi_l, T& chi_h, T& H0,
                            T& PsiSolver_K) {
#ifdef MPI_ENABLED
  mpi().bCast(mu_l, 0);
  mpi().bCast(mu_h, 0);
  mpi().bCast(chi_l, 0);
  mpi().bCast(chi_h, 0);
  mpi().bCast(H0, 0);
  mpi().bCast(PsiSolver_K, 0);
#endif
  lattice.template getField<MU_L<T>>().InitValue(mu_l);
  lattice.template getField<MU_H<T>>().InitValue(mu_h);
  lattice.template getField<CHI_L<T>>().InitValue(chi_l);
  lattice.template getField<CHI_H<T>>().InitValue(chi_h);
  lattice.template getField<H_0<T>>().InitValue(H0);
  lattice.template getField<PSI_K<T>>().InitValue(PsiSolver_K);
}

}  // namespace mfield

#include "mfield/mfield3d.hh"
