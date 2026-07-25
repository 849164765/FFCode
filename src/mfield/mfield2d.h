#pragma once

#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"

namespace mfield {

// ===================================================================
//  Section 1: Self-Owned Field Base Definitions
//
//  GenericField<GenericArray<T>> — per-cell varying (ghost sync):
//    PSI, OMEGA_PSI, MU_PERCELL, CHI_PERCELL, HX, HY, HMAG
//  Data<T> — block-level constant:
//    MU_L, MU_H, CHI_L, CHI_H, H0
// ===================================================================

// --- Per-cell varying fields ---
struct PsiBase         : public FieldBase<1> {};  // 磁标势 ψ
struct OmegaPsiBase    : public FieldBase<1> {};  // ω_ψ = 1/τ_ψ
struct MuPerCellBase   : public FieldBase<1> {};  // μ(φ) = μ_l + φ·(μ_h - μ_l)
struct ChiPerCellBase  : public FieldBase<1> {};  // χ(φ) = χ_l + φ·(χ_h - χ_l)
struct HxBase          : public FieldBase<1> {};  // H_x = -∂ψ/∂x
struct HyBase          : public FieldBase<1> {};  // H_y = -∂ψ/∂y
struct HMagBase        : public FieldBase<1> {};  // |H| = sqrt(H_x² + H_y²)

// --- Block-level constants ---
struct MuLBase    : public FieldBase<1> {};  // 铁磁相磁导率 μ_l
struct MuHBase    : public FieldBase<1> {};  // 非磁相磁导率 μ_h
struct ChiLBase   : public FieldBase<1> {};  // 铁磁相磁化率 χ_l
struct ChiHBase   : public FieldBase<1> {};  // 非磁相磁化率 χ_h
struct H0Base     : public FieldBase<1> {};  // 外加磁场强度 H₀ (沿 y 方向)

// ===================================================================
//  Section 2: Self-Owned Field Type Aliases
// ===================================================================

template <typename T>
using PSI         = GenericField<GenericArray<T>, PsiBase>;
template <typename T>
using OMEGA_PSI   = GenericField<GenericArray<T>, OmegaPsiBase>;
template <typename T>
using MU_PERCELL  = GenericField<GenericArray<T>, MuPerCellBase>;
template <typename T>
using CHI_PERCELL = GenericField<GenericArray<T>, ChiPerCellBase>;
template <typename T>
using HX          = GenericField<GenericArray<T>, HxBase>;
template <typename T>
using HY          = GenericField<GenericArray<T>, HyBase>;
template <typename T>
using HMAG        = GenericField<GenericArray<T>, HMagBase>;

template <typename T>
using MU_L  = Data<T, MuLBase>;
template <typename T>
using MU_H  = Data<T, MuHBase>;
template <typename T>
using CHI_L = Data<T, ChiLBase>;
template <typename T>
using CHI_H = Data<T, ChiHBase>;
template <typename T>
using H_0   = Data<T, H0Base>;

// ===================================================================
//  Section 3: Field Packs
// ===================================================================

template <typename T>
using MFFIELDS = TypePack<
    PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>,
    HX<T>, HY<T>, HMAG<T>,
    MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>, H_0<T>>;

template <typename T, unsigned int D>
using MFEXTERNALFIELDS = TypePack<PHI<T>>;

template <typename T, unsigned int D>
using MFFIELDPACK = TypePack<MFFIELDS<T>, MFEXTERNALFIELDS<T, D>>;

// ===================================================================
//  Section 4: Functor Declarations
// ===================================================================

// MF1: PF → MF coupling — update per-cell μ(φ), χ(φ), ω_ψ from φ
template <typename PFCELL, typename MFCELL>
struct MFUpdateCoeffs2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell);
};

// MF2: H = -∇ψ via D2Q5 4-direction isotropic gradient
template <typename CELL>
struct MFComputeH2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// MF3a: Magnetic body force from Maxwell stress tensor divergence
//   F_m = (H·∇χ)H - (1/2)|H|²∇χ
template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

// MF3b: Interfacial force only: F_m = -(1/2)|H|²∇χ
template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForceInterfacial2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

// MF3c: Paper Eq.(8) formula: F_m = (χ/2)∇|H|²  (Kelvin body force)
// Valid when magnetic susceptibility is constant (sharp interface limit).
template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForcePaper2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

// ===================================================================
//  Section 5: Communication Functions
// ===================================================================

template <typename T, typename LATTICE>
void CommunicatePSI(LATTICE& lattice) {
  lattice.template getField<PSI<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateHX(LATTICE& lattice) {
  lattice.template getField<HX<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateHY(LATTICE& lattice) {
  lattice.template getField<HY<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateHMAG(LATTICE& lattice) {
  lattice.template getField<HMAG<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateOMEGAPSI(LATTICE& lattice) {
  lattice.template getField<OMEGA_PSI<T>>().Communicate();
}

// Batch communicate all per-cell magnetic fields that need ghost sync
template <typename T, typename LATTICE>
void CommunicateAllMFFields(LATTICE& lattice) {
  lattice.template getField<PSI<T>>().Communicate();
  lattice.template getField<HX<T>>().Communicate();
  lattice.template getField<HY<T>>().Communicate();
  lattice.template getField<HMAG<T>>().Communicate();
  lattice.template getField<OMEGA_PSI<T>>().Communicate();
  lattice.template getField<CHI_PERCELL<T>>().Communicate();
}

// Broadcast block-level params from rank 0
template <typename T, typename LATTICE>
void BroadcastAllMFParams(LATTICE& lattice,
                          T& mu_l, T& mu_h, T& chi_l, T& chi_h, T& H0) {
#ifdef MPI_ENABLED
  mpi().bCast(mu_l, 0);
  mpi().bCast(mu_h, 0);
  mpi().bCast(chi_l, 0);
  mpi().bCast(chi_h, 0);
  mpi().bCast(H0, 0);
#endif
  lattice.template getField<MU_L<T>>().InitValue(mu_l);
  lattice.template getField<MU_H<T>>().InitValue(mu_h);
  lattice.template getField<CHI_L<T>>().InitValue(chi_l);
  lattice.template getField<CHI_H<T>>().InitValue(chi_h);
  lattice.template getField<H_0<T>>().InitValue(H0);
}

}  // namespace mfield

#include "mfield/mfield2d.hh"
