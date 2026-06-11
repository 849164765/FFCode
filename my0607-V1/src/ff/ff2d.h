#pragma once

#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"

namespace ff {

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

// ===================================================================
//  Section 3: Field Packs (self-owned + external, cf. CAFIELDS + REFFIELDS)
// ===================================================================

template <typename T>
using FFFIELDS = TypePack<LAPLACIAN<T>, CHEMICALPOTENTIAL<T>,
                          GRAVITY<T>, BETA<T>, KAPPA<T>,
                          RHO_L<T>, RHO_H<T>, ETA_L<T>, ETA_H<T>>;

template <typename T, unsigned int D>
using FFEXTERNALFIELDS = TypePack<PHI<T>, GRAD<T, D>, NORMAL<T, D>, INTERFACEWIDTH<T>>;

template <typename T, unsigned int D>
using FFFIELDPACK = TypePack<FFFIELDS<T>, FFEXTERNALFIELDS<T, D>>;

template <typename T, unsigned int D>
using ALLFF_FIELDS = typename ExtractFieldPack<FFFIELDPACK<T, D>>::mergedpack;

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

template <typename PFCELL, typename NSCELL>
struct FFRhoOmegaUpdate2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// ---- Post-stream macroscopic computation ----
// Compute rho and preliminary velocity u* = Σc_i f_i / rho from distributions.
// Writes ρ→RHO and u*→VELOCITY. No force correction: u_final = u* + F/(2ρ)
// must be applied separately (via FFFinalVelocityCompute2D).
template <typename CELL>
struct FFRhoUCompute2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// Force-corrected velocity: u_final = u* + F/(2ρ)
// Reads FORCE from cell field, updates VELOCITY by adding force correction.
template <typename CELL>
struct FFFinalVelocityCompute2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// Density-gradient forces: F_p + F_v (coupled PF→NS)
// Eq.(27) F_v = ν(∇u+u∇)·∇ρ   Eq.(28) F_p = -(p/ρ)∇ρ
// Adds both terms to ns_cell FORCE field.
template <typename PFCELL, typename NSCELL>
struct FFDensityForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// Interpolate ρ(φ)=ρ_l+φ(ρ_h-ρ_l) → NS RHO and ω(φ)=1/(0.5+ν(φ)/cs²) → NS OMEGA
// Called after PF processing, before NS collision.
template <typename PFCELL, typename NSCELL>
struct FFRhoOmegaInterp2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename PFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, NSCELL& ns_cell);
};

// Compute u* = Σc_i·g_i / Σg_i from post-stream populations.
// Writes VELOCITY only; does NOT overwrite RHO or OMEGA.
template <typename CELL>
struct FFVelocityCompute2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

// Velocity + Pressure per Guo et al. Eqs.(35)-(36)
//  u = Σ c_i·g_i + F/(2·RHO)                  (Eq.35)
//  p = RHO * cs² * Σ g_i                      (Eq.36)
// Reads NS POP, RHO(=ρ(φ)), FORCE; writes VELOCITY and PRESSURE.
template <typename CELL>
struct FFVelocityPressureCompute2D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
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

}  // namespace ff

#include "ff/ff2d.hh"
#include "ff/ff3d.hh"
