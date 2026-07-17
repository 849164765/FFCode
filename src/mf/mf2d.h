#pragma once

#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"

namespace mf {

struct H2Base : public FieldBase<1> {};
struct GradH2Base : public FieldBase<2> {};
struct ChiCellBase : public FieldBase<1> {};
struct ChiFerroBase : public FieldBase<1> {};
struct Mu0Base : public FieldBase<1> {};
struct H0Base : public FieldBase<1> {};

template <typename T> using H2 = GenericField<GenericArray<T>, H2Base>;
template <typename T> using GRAD_H2 = GenericField<GenericArray<Vector<T, 2>>, GradH2Base>;
template <typename T> using CHI_CELL = GenericField<GenericArray<T>, ChiCellBase>;
template <typename T> using CHI_FERRO = Data<T, ChiFerroBase>;
template <typename T> using MU0 = Data<T, Mu0Base>;
template <typename T> using H0_FIELD = Data<T, H0Base>;

template <typename T>
using MFSELFFIELDS = TypePack<H2<T>, GRAD_H2<T>, CHI_CELL<T>, CHI_FERRO<T>, MU0<T>, H0_FIELD<T>>;

// ===================================================================
//  Functor Declarations
// ===================================================================

// MFCollision: D2Q5 MRT collision for Laplace equation (no source term)
template <typename CELL>
struct MFCollision {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static_assert(LatSet::q == 5, "MFCollision requires D2Q5 lattice");
  __any__ static void apply(CELL& cell);
};

// MFMacro: ψ = Σh_α from D2Q5 populations
template <typename CELL>
struct MFMacro {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  static_assert(LatSet::q == 5, "MFMacro requires D2Q5 lattice");
  __any__ static void apply(CELL& cell);
};

// MFHField: H = -∇ψ, |H|² = Hx²+Hy² (4-point central difference)
template <typename MFCELL>
struct MFHField {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  static_assert(LatSet::q == 5, "MFHField requires D2Q5 lattice");
  __any__ static void apply(MFCELL& cell);
};

// MFGradH2: ∇(|H|²) via 4-point central differencing
template <typename MFCELL>
struct MFGradH2 {
  using T = typename MFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  static_assert(LatSet::q == 5, "MFGradH2 requires D2Q5 lattice");
  __any__ static void apply(MFCELL& cell);
};

// MFChiUpdate: PF → MF coupling: χ_cell = χ_ferro · φ
template <typename PFCELL, typename MFCELL>
struct MFChiUpdate {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell);
};

// MFKelvinForce: MF → NS coupling: F_m = (μ₀·χ/2)·∇(|H|²)
template <typename MFCELL, typename NSCELL>
struct MFKelvinForce {
  using T = typename MFCELL::FloatType;
  using LatSet = typename NSCELL::LatticeSet;
  __any__ static void apply(MFCELL& mf_cell, NSCELL& ns_cell);
};

// MFInitPsi: initialize ψ from applied H₀ field: ψ = -H₀·y + const
// (ψ satisfies H_y = -∂ψ/∂y = H₀)
template <typename CELL>
struct MFInitPsi {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;
  __any__ static void apply(CELL& cell);
};

// ===================================================================
//  Communication Functions
// ===================================================================

template <typename T, typename LATTICE>
void CommunicateAllMFSelfFields(LATTICE& lattice) {
  lattice.template getField<H2<T>>().Communicate();
  lattice.template getField<GRAD_H2<T>>().Communicate();
  lattice.template getField<CHI_CELL<T>>().Communicate();
}

template <typename T, typename LATTICE>
void BroadcastAllMFParams(LATTICE& lattice,
                          T& chi_ferro, T& mu0, T& h0) {
#ifdef MPI_ENABLED
  mpi().bCast(chi_ferro, 0);
  mpi().bCast(mu0, 0);
  mpi().bCast(h0, 0);
#endif
  lattice.template getField<CHI_FERRO<T>>().InitValue(chi_ferro);
  lattice.template getField<MU0<T>>().InitValue(mu0);
  lattice.template getField<H0_FIELD<T>>().InitValue(h0);
}

}  // namespace mf

#include "mf/mf2d.hh"
