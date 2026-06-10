#pragma once

#include "data_struct/block_lattice_base.h"
#include "data_struct/field_struct.h"
#include "lbm/lattice_set.h"
#include "lbm/unit_converter.h"

namespace ff {

// --- FF-specific field bases ---
struct GRADBase : public FieldBase<1> {};
struct NORMALBase : public FieldBase<1> {};
struct LAPLACIANBase : public FieldBase<1> {};
struct CHEMICALPOTENTIALBase : public FieldBase<1> {};
struct ACSOURCEBase : public FieldBase<1> {};
struct FORCE_SF_Base : public FieldBase<1> {};
struct FORCE_Buoy_Base : public FieldBase<1> {};
struct FORCE_Visc_Base : public FieldBase<1> {};
struct FORCE_P_Base : public FieldBase<1> {};
struct RHO_LBase : public FieldBase<1> {};
struct RHO_HBase : public FieldBase<1> {};
struct ETA_LBase : public FieldBase<1> {};
struct ETA_HBase : public FieldBase<1> {};
struct GRAVITYBase : public FieldBase<1> {};

// --- FF type aliases (2D only) ---
template <typename T>
using GRAD = GenericField<GenericArray<Vector<T, 2>>, GRADBase>;
template <typename T>
using NORMAL = GenericField<GenericArray<Vector<T, 2>>, NORMALBase>;
template <typename T>
using LAPLACIAN = GenericField<GenericArray<T>, LAPLACIANBase>;
template <typename T>
using CHEMICALPOTENTIAL = GenericField<GenericArray<T>, CHEMICALPOTENTIALBase>;
template <typename T>
using ACSOURCE = GenericField<GenericArray<Vector<T, 2>>, ACSOURCEBase>;
template <typename T>
using FORCE_SF = GenericField<GenericArray<Vector<T, 2>>, FORCE_SF_Base>;
template <typename T>
using FORCE_Buoy = GenericField<GenericArray<Vector<T, 2>>, FORCE_Buoy_Base>;
template <typename T>
using FORCE_Visc = GenericField<GenericArray<Vector<T, 2>>, FORCE_Visc_Base>;
template <typename T>
using FORCE_P = GenericField<GenericArray<Vector<T, 2>>, FORCE_P_Base>;
template <typename T>
using RHO_L = Data<T, RHO_LBase>;
template <typename T>
using RHO_H = Data<T, RHO_HBase>;
template <typename T>
using ETA_L = Data<T, ETA_LBase>;
template <typename T>
using ETA_H = Data<T, ETA_HBase>;
template <typename T>
using GRAVITY = Data<T, GRAVITYBase>;
// --- Field packs (2D) ---
// Owned by FF2DMgr (computed fields)
template <typename T>
using FFFIELDS = TypePack<NORMAL<T>, LAPLACIAN<T>, ACSOURCE<T>,
                            FORCE_SF<T>, FORCE_Buoy<T>,
                            FORCE_Visc<T>, FORCE_P<T>>;

// Referenced from PFLattice/NSLattice
template <typename T>
using REFLBMFIELDS = TypePack<PHI<T>, GRAD<T>, CHEMICALPOTENTIAL<T>,
                              RHO_L<T>, RHO_H<T>,
                              ETA_L<T>, ETA_H<T>, GRAVITY<T>,
                              RHO<T>, PRESSURE<T>,
                              FORCE<T, 2>, VELOCITY<T, 2>>;

template <typename T>
using FIELDPACK = TypePack<FFFIELDS<T>, REFLBMFIELDS<T>>;

template <typename T>
using ALLFIELDS = typename ExtractFieldPack<FIELDPACK<T>>::mergedpack;

// --- Allen-Cahn source scheme ---
// Discrete source term: S_i = w_i · (c_i · V)
//   V = (4φ(1-φ)/W)·n — pre-computed by FF2DBlock::computeACSource()
template <typename CELL>
struct PhaseFieldSource {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  static constexpr unsigned int scalardir = 0;

  __any__ static void apply(const Vector<T, LatSet::d>& u,
                            const Vector<T, LatSet::d>& V,
                            std::array<T, LatSet::q>& Si) {
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cdotV = T{0};
      for (unsigned int l = 0; l < LatSet::d; ++l)
        cdotV += latset::c<LatSet>(i)[l] * V[l];
      Si[i] = latset::w<LatSet>(i) * cdotV;
    }
  }

  __any__ static auto& getForce(CELL& cell) {
    return cell.template get<ACSOURCE<T>>();
  }
};

/// Single-block FF2D: computes per-block phase-field quantities
template <typename T, typename LatSet>
class FF2DBlock : public BlockLatticeBase<T, LatSet, ALLFIELDS<T>> {
 private:
  PhaseFieldConverter<T>& _conv;

 public:
  template <typename... FIELDPTRS>
  FF2DBlock(Block<T, 2>& geo, PhaseFieldConverter<T>& conv,
            std::tuple<FIELDPTRS...> fieldptrs);

  void computeGradient();
  void computeLaplacian();
  void computeChemicalPotential();
  void computeACSource();
  void computeForceSF();
  void computeForceBuoy();
  void computeForceVisc();
  void computeForceP();
  void computeForceTotal();
  void apply();

  AbstractConverter<T>& converter() { return _conv; }
  const AbstractConverter<T>& converter() const { return _conv; }
};

/// Multi-block FF2D manager
template <typename T, typename LatSet>
class FF2DManager : public BlockLatticeManagerBase<T, LatSet, FIELDPACK<T>> {
 private:
  PhaseFieldConverter<T>& _conv;
  std::vector<FF2DBlock<T, LatSet>> _blocks;

 public:
  template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
  FF2DManager(BlockGeometry<T, 2>& blockgeo, PhaseFieldConverter<T>& conv,
              INITVALUEPACK& initvalues, FIELDPTRTYPES*... fieldptrs);

  void Init(BlockGeometryHelper<T, 2>& GeoHelper);
  void computeGradient();
  void computeLaplacian();
  void computeChemicalPotential();
  void computeACSource();
  void computeForceSF();
  void computeForceBuoy();
  void computeForceVisc();
  void computeForceP();
  void computeForceTotal();
  void apply();
  void Communicate();
};

}  // namespace ff

#include "ff/ff2d.hh"
#include "ff/ff3d.hh"
