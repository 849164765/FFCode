#pragma once

#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"

namespace mfield {

// ===================================================================
//  Section 1: 3D-specific Field Base Definitions
//
//  Reuse all 2D field bases from mfield2d.h (PSI, OMEGA_PSI, MU_PERCELL,
//  CHI_PERCELL, HX, HY, HMAG, MU_L, MU_H, CHI_L, CHI_H, H_0).
//  Only add the z-component of H for 3D.
// ===================================================================

// --- Per-cell varying field (3D only) ---
struct HzBase : public FieldBase<1> {};  // H_z = -∂ψ/∂z

// ===================================================================
//  Section 2: 3D-specific Field Type Aliases
// ===================================================================

template <typename T>
using HZ = GenericField<GenericArray<T>, HzBase>;

// ===================================================================
//  Section 3: 3D Field Packs
// ===================================================================

template <typename T>
using MFFIELDS3D = TypePack<
    PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>,
    HX<T>, HY<T>, HZ<T>, HMAG<T>,
    MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>, H_0<T>>;

template <typename T, unsigned int D>
using MFEXTERNALFIELDS3D = TypePack<PHI<T>>;

template <typename T, unsigned int D>
using MFFIELDPACK3D = TypePack<MFFIELDS3D<T>, MFEXTERNALFIELDS3D<T, D>>;

// ===================================================================
//  Section 4: 3D Functor Declarations
// ===================================================================

// MF1-3D: PF → MF coupling — update per-cell μ(φ), χ(φ), ω_ψ from φ
//         Same formula as MFUpdateCoeffs2D (dimension-agnostic)
template <typename PFCELL, typename MFCELL>
struct MFUpdateCoeffs3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell);
};

// MF2-3D: H = -∇ψ via D3Q7 6-direction isotropic gradient
template <typename CELL>
struct MFComputeH3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// MF3-3D: Magnetic body force F_mag = χ·|H|·∇|H| (Kelvin force)
template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

// ===================================================================
//  Section 5: Communication Functions (3D)
// ===================================================================

template <typename T, typename LATTICE>
void CommunicateHZ(LATTICE& lattice) {
  lattice.template getField<HZ<T>>().Communicate();
}

// Batch communicate all per-cell magnetic fields that need ghost sync (3D)
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

}  // namespace mfield

#include "mfield/mfield3d.hh"
