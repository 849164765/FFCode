#pragma once

#include "mfield/mfield2d.h"

namespace mfield {

// The scalar magnetic fields and phase-independent parameters are shared with
// the 2D implementation.  A separate 3D field pack only adds H_z so both
// headers can be included by freelb.h without redefining the common aliases.
struct HzBase : public FieldBase<1> {};

template <typename T>
using HZ = GenericField<GenericArray<T>, HzBase>;

template <typename T>
using MFFIELDS3D = TypePack<
    PSI<T>, OMEGA_PSI<T>, MU_PERCELL<T>, CHI_PERCELL<T>,
    HX<T>, HY<T>, HZ<T>, HMAG<T>,
    MU_L<T>, MU_H<T>, CHI_L<T>, CHI_H<T>, H_0<T>, PSI_K<T>>;

template <typename PFCELL, typename MFCELL>
struct MFUpdateCoeffs3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell);
};

template <typename CELL>
struct MFComputeH3D {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  __any__ static void apply(CELL& cell);
};

template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MFMagneticForce3D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;
  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell);
};

template <typename T, typename LATTICE>
void CommunicateHZ(LATTICE& lattice) {
  lattice.template getField<HZ<T>>().Communicate();
}

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
