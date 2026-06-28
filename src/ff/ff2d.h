#pragma once

#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"

namespace ff {

// ===================================================================
//  Utility: quintic smooth-step (Heaviside-like) function
//  h(φ) = φ³(6φ² - 15φ + 10)
//  Satisfies: h(0)=0, h(1)=1, h'(0)=h'(1)=0, h''(0)=h''(1)=0
//  Provides smooth transition of magnetic properties across the interface.
// ===================================================================
template <typename T>
__any__ T quinticSmoothStep(T phi) {
  T phi2 = phi * phi;
  T phi3 = phi2 * phi;
  return phi3 * (T{6} * phi2 - T{15} * phi + T{10});
}

// ===================================================================
//  Section 1: Self-Owned Field Base Definitions
//  (cf. src/ca/zhu_stefanescu2d.h: STATEBase, FSBase, ...)
//
//  GenericField<GenericArray<T>> — per-cell varying (ghost sync via Communicate):
//    LaplacianBase, ChemPotentialBase
//  Data<T> — block-level constant (no ghost cells, same mechanism as CONSTFORCE):
//    GravityBase, BetaBase, KappaBase, RhoLBase, RhoHBase, EtaLBase, EtaHBase
// ===================================================================

struct LaplacianBase : public FieldBase<1> {};
struct ChemPotentialBase : public FieldBase<1> {};

struct GravityBase : public FieldBase<1> {};
struct BetaBase : public FieldBase<1> {};
struct KappaBase : public FieldBase<1> {};
struct RhoLBase : public FieldBase<1> {};
struct RhoHBase : public FieldBase<1> {};
struct EtaLBase : public FieldBase<1> {};
struct EtaHBase : public FieldBase<1> {};
struct DeltaRhoBase : public FieldBase<1> {};

// Magnetic field bases (Data<T> — block-level constants)
struct MuLBase : public FieldBase<1> {};
struct MuHBase : public FieldBase<1> {};
struct ChiLBase : public FieldBase<1> {};
struct ChiHBase : public FieldBase<1> {};

// Magnetic field bases (GenericField<GenericArray<T>> — per-cell ghost fields)
struct PsiBase : public FieldBase<1> {};
struct HFieldBase : public FieldBase<1> {};
struct HSqBase : public FieldBase<1> {};
struct MuBase : public FieldBase<1> {};         // interpolated μ on PF lattice
struct ChiBase : public FieldBase<1> {};         // interpolated χ on MF lattice
struct FmBase : public FieldBase<1> {};          // magnetic force F_m = (χ/2)·∇|H|²

// ===================================================================
//  Section 2: Self-Owned Field Type Aliases
// ===================================================================

template <typename T>
using LAPLACIAN = GenericField<GenericArray<T>, LaplacianBase>;

template <typename T>
using CHEMICALPOTENTIAL = GenericField<GenericArray<T>, ChemPotentialBase>;

template <typename T>
using GRAVITY = Data<T, GravityBase>;

template <typename T>
using BETA = Data<T, BetaBase>;

template <typename T>
using KAPPA = Data<T, KappaBase>;

template <typename T>
using RHO_L = Data<T, RhoLBase>;

template <typename T>
using RHO_H = Data<T, RhoHBase>;

template <typename T>
using ETA_L = Data<T, EtaLBase>;

template <typename T>
using ETA_H = Data<T, EtaHBase>;

template <typename T>
using DELTARHO = Data<T, DeltaRhoBase>;

// Magnetic field: block-level constants (Data<T>)
template <typename T>
using MU_L = Data<T, MuLBase>;
template <typename T>
using MU_H = Data<T, MuHBase>;
template <typename T>
using CHI_L = Data<T, ChiLBase>;
template <typename T>
using CHI_H = Data<T, ChiHBase>;

// Magnetic field: per-cell ghost fields (GenericField<GenericArray<T>>)
template <typename T>
using PSI = GenericField<GenericArray<T>, PsiBase>;
template <typename T, unsigned int D>
using H_FIELD = GenericField<GenericArray<Vector<T, D>>, HFieldBase>;
template <typename T>
using H_SQ = GenericField<GenericArray<T>, HSqBase>;
template <typename T>
using MU = GenericField<GenericArray<T>, MuBase>;
template <typename T>
using CHI = GenericField<GenericArray<T>, ChiBase>;

template <typename T, unsigned int D>
using F_MAGNETIC = GenericField<GenericArray<Vector<T, D>>, FmBase>;

// ===================================================================
//  Section 3: Field Packs (self-owned + external, cf. CAFIELDS + REFFIELDS)
// ===================================================================

template <typename T>
using FFFIELDS = TypePack<LAPLACIAN<T>, CHEMICALPOTENTIAL<T>,
                          GRAVITY<T>, BETA<T>, KAPPA<T>,
                          RHO_L<T>, RHO_H<T>, ETA_L<T>, ETA_H<T>,
                          DELTARHO<T>>;

template <typename T, unsigned int D>
using FFEXTERNALFIELDS = TypePack<PHI<T>, GRAD<T, D>, NORMAL<T, D>, INTERFACEWIDTH<T>>;

template <typename T, unsigned int D>
using FFFIELDPACK = TypePack<FFFIELDS<T>, FFEXTERNALFIELDS<T, D>>;

template <typename T, unsigned int D>
using ALLFF_FIELDS = typename ExtractFieldPack<FFFIELDPACK<T, D>>::mergedpack;

// ===================================================================
//  Magnetic field field packs
//
//  FFMAGSELFFIELDS:  per-cell ghost fields on MF lattice (PSI, H_FIELD, H_SQ)
//  FFMAGCONSTFIELDS: block-level constants on MF lattice (MU_L, MU_H, CHI_L, CHI_H)
//  FFPFFIELDS_EXT:   per-cell field on PF lattice (MU — interpolated permeability)
//                    CHI is stored on MF lattice, computed in the PF→MF copy step
// ===================================================================
template <typename T, unsigned int D>
using FFMAGSELFFIELDS = TypePack<PSI<T>, H_FIELD<T, D>, H_SQ<T>, CHI<T>>;

template <typename T>
using FFMAGCONSTFIELDS = TypePack<MU_L<T>, MU_H<T>>;

template <typename T>
using FFMFCONSTFIELDS = TypePack<CHI_L<T>, CHI_H<T>>;

template <typename T>
using FFPFFIELDS_EXT = TypePack<MU<T>>;

// ===================================================================
//  Section 4: Functor Declarations
// ===================================================================

template <typename CELL>
struct FF2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FF3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FFLaplacian2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FFChemPotential2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FFLaplacian3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FFChemPotential3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

template <typename CELL>
struct FFChemPotentialGradient2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

template <typename PFCELL, typename NSCELL>
struct FFSurfaceTension2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFGravityForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// FFViscoForce2D: F_v from non-equilibrium moments (Fortran MRT two-pass algorithm)
// Computes mgneq from NS populations using first-pass MRT relaxation,
// then F_v = -3*mu*DeltaRho/rho * (C · grad_phi)
template <typename PFCELL, typename NSCELL>
struct FFViscoForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// FFPreForce2D: F_p = -(p/3) * DeltaRho * grad_phi (pressure gradient force)
template <typename PFCELL, typename NSCELL>
struct FFPreForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFRhoOmegaUpdate2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFSurfaceTension3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFGravityForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFPreForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

template <typename PFCELL, typename NSCELL>
struct FFRhoOmegaUpdate3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// FFViscoForce3D: F_v from non-equilibrium moments (Fortran MRT two-pass algorithm)
// M version: 3D generic version using M/InvM matrix operations.
template <typename PFCELL, typename NSCELL>
struct FFViscoForce3DM {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};
template <typename PFCELL, typename NSCELL>
struct FFViscoForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// ===================================================================
//  Magnetic Field Functor Declarations
//
//  MFGradient2D:  H = -∇ψ using second-order central difference on D2Q5
//  MFHsq2D:       |H|² = H_x² + H_y²
//  MFForce2D:     F_m = (χ/2) ∇(|H|²) → add to NSCELL::FORCE
//  FFMuUpdate2D:  μ = μ_l + (μ_h - μ_l) * h(φ)  with h(φ)=quinticSmoothStep
//  Quintic smooth-step replaces linear interpolation for better stability
// ===================================================================

template <typename MFCELL>
struct MFGradient2D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell);
};

template <typename MFCELL>
struct MFHsq2D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell);
};

template <typename MFCELL, typename NSCELL>
struct MFForce2D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell, NSCELL& ns_cell);
};

template <typename PFCELL>
struct FFMuUpdate2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  using GenericRho = typename PFCELL::GenericRho;
  __any__ static void apply(PFCELL& pf_cell);
};

template <typename MFCELL>
struct MFGradient3D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell);
};

template <typename MFCELL>
struct MFHsq3D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell);
};

template <typename MFCELL, typename NSCELL>
struct MFForce3D {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell, NSCELL& ns_cell);
};

template <typename PFCELL>
struct FFMuUpdate3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  using GenericRho = typename PFCELL::GenericRho;
  __any__ static void apply(PFCELL& pf_cell);
};

// ===================================================================
//  Section 5: Communication Functions for Self-Owned Fields
//  (cf. src/ca/zhu_stefanescu2d.h: BlockZhuStefanescu2DManager::Communicate())
//
//  GenericField<GenericArray<T>> (per-cell, ghost cells):
//    CommunicateAllSelfFields — ghost sync each time step (CHEMICALPOTENTIAL)
//
//  Data<T> (block-level constant, same mechanism as CONSTFORCE/CONSTRHO):
//    BroadcastAllParams — MPI_Bcast from rank 0, called once after construction.
//    Data fields have no ghost cells (isField=false), no Communicate() needed.
//    Broadcast covers environments where non-root ranks cannot read the INI file.
// ===================================================================

template <typename T, typename LATTICE>
void CommunicateLAPLACIAN(LATTICE& lattice) {
  lattice.template getField<LAPLACIAN<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateCHEMICALPOTENTIAL(LATTICE& lattice) {
  lattice.template getField<CHEMICALPOTENTIAL<T>>().Communicate();
}

// Batch communicate Type A self-owned fields that need ghost sync.
template <typename T, typename LATTICE>
void CommunicateAllSelfFields(LATTICE& lattice) {
  lattice.template getField<CHEMICALPOTENTIAL<T>>().Communicate();
  // LAPLACIAN: not needed — FFChemPotential2D reads from same cell only
}

// Broadcast Type B params from rank 0 and set all uniform fields.
// Must be called once after lattice construction (before the main loop).
template <typename T, typename LATTICE>
void BroadcastAllParams(LATTICE& lattice,
                        T& rho_l, T& rho_h, T& eta_l, T& eta_h,
                        T& gravity, T& Beta, T& Kappa) {
#ifdef MPI_ENABLED
  mpi().bCast(rho_l, 0);
  mpi().bCast(rho_h, 0);
  mpi().bCast(eta_l, 0);
  mpi().bCast(eta_h, 0);
  mpi().bCast(gravity, 0);
  mpi().bCast(Beta, 0);
  mpi().bCast(Kappa, 0);
#endif
  lattice.template getField<RHO_L<T>>().InitValue(rho_l);
  lattice.template getField<RHO_H<T>>().InitValue(rho_h);
  lattice.template getField<ETA_L<T>>().InitValue(eta_l);
  lattice.template getField<ETA_H<T>>().InitValue(eta_h);
  lattice.template getField<GRAVITY<T>>().InitValue(gravity);
  lattice.template getField<BETA<T>>().InitValue(Beta);
  lattice.template getField<KAPPA<T>>().InitValue(Kappa);
}

// Magnetic field communication functions
template <typename T, typename LATTICE>
void CommunicatePSI(LATTICE& lattice) {
  lattice.template getField<PSI<T>>().Communicate();
}

template <typename T, typename LATTICE>
void CommunicateHSQ(LATTICE& lattice) {
  lattice.template getField<H_SQ<T>>().Communicate();
}

// Broadcast magnetic susceptibility params from rank 0 and set uniform fields.
// Must be called once after MF lattice construction.
// Note: MU_L, MU_H are set on PF lattice separately (used by FFMuUpdate2D).
template <typename T, typename LATTICE>
void BroadcastMagParams(LATTICE& lattice,
                        T& chi_l, T& chi_h) {
#ifdef MPI_ENABLED
  mpi().bCast(chi_l, 0);
  mpi().bCast(chi_h, 0);
#endif
  lattice.template getField<CHI_L<T>>().InitValue(chi_l);
  lattice.template getField<CHI_H<T>>().InitValue(chi_h);
}

}  // namespace ff

#include "ff/ff2d.hh"
#include "ff/ff3d.hh"
