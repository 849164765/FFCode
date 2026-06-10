#pragma once

#include <vector>
#include "ff/ff2d.h"

namespace ff {

// ===================================================================
// FF2DBlock implementations
// ===================================================================

template <typename T, typename LatSet>
template <typename... FIELDPTRS>
FF2DBlock<T, LatSet>::FF2DBlock(Block<T, 2>& geo,
                                PhaseFieldConverter<T>& conv,
                                std::tuple<FIELDPTRS...> fieldptrs)
    : BlockLatticeBase<T, LatSet, ALLFIELDS<T>>(geo, fieldptrs),
      _conv(conv) {}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeGradient() {
  // writes to GRAD (referenced from PFLattice), NORMAL (owned)
  auto& gradField = this->template getField<GRAD<T>>();
  auto& normalField = this->template getField<NORMAL<T>>();
  auto& phiField = this->template getField<PHI<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;

      Vector<T, LatSet::d> grad;
      grad[0] = T{0};
      grad[1] = T{0};

      for (unsigned int k = 1; k < LatSet::q; ++k) {
        std::size_t nid = id + this->Delta_Index[k];
        T phi_k = phiField.get(nid);
        T wk = latset::w<LatSet>(k);
        const auto& ck = latset::c<LatSet>(k);
        grad[0] += wk * ck[0] * phi_k;
        grad[1] += wk * ck[1] * phi_k;
      }
      grad[0] /= LatSet::cs2;
      grad[1] /= LatSet::cs2;

      gradField.SetField(id, grad);

      T gradMag = grad.getnorm();
      T delta = T(0.005);
      Vector<T, LatSet::d> n;
      if (gradMag < delta) {
        n[0] = T{0};
        n[1] = T{0};
      } else {
        n[0] = grad[0] / gradMag;
        n[1] = grad[1] / gradMag;
      }
      normalField.SetField(id, n);
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeLaplacian() {
  auto& laplacianField = this->template getField<LAPLACIAN<T>>();
  auto& phiField = this->template getField<PHI<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      T phi_self = phiField.get(id);
      T laplacian = T{0};

      for (unsigned int k = 0; k < LatSet::q; ++k) {
        std::size_t nid = id + this->Delta_Index[k];
        T phi_k = phiField.get(nid);
        T wk = latset::w<LatSet>(k);
        laplacian += wk * (phi_k - phi_self);
      }
      laplacian *= T{2} / LatSet::cs2;

      laplacianField.SetField(id, laplacian);
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeChemicalPotential() {
  // λ = 4β * φ*(φ-1)*(φ-0.5) - κ * ∇²φ
  // writes to CHEMICALPOTENTIAL (referenced from PFLattice)
  auto& chemPotField = this->template getField<CHEMICALPOTENTIAL<T>>();
  auto& phiField = this->template getField<PHI<T>>();
  auto& lapField = this->template getField<LAPLACIAN<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      T phi = phiField.get(id);
      T laplacian = lapField.get(id);

      T double_well = phi * (phi - T{1}) * (phi - T{0.5});
      T chem_pot = T{4} * _conv.getLatticeBeta() * double_well
                    - _conv.getLatticeKappa() * laplacian;

      chemPotField.SetField(id, chem_pot);
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeACSource() {
  // ACSOURCE = (4φ(1-φ)/W) · NORMAL  (Allen-Cahn source vector)
  auto& acSource = this->template getField<ACSOURCE<T>>();
  auto& phiField = this->template getField<PHI<T>>();
  auto& normalField = this->template getField<NORMAL<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  T W = _conv.getLatticeW();
  std::cout << "[computeACSource] W=" << W << std::endl;

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      T phi = phiField.get(id);
      T S = T{4} * phi * (T{1} - phi) / W;
      const Vector<T, LatSet::d>& n = normalField.get(id);
      acSource.SetField(id, n * S);
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeForceSF() {
  // Surface tension force: F_s = λ · ∇φ
  auto& forceSF = this->template getField<FORCE_SF<T>>();
  auto& chemField = this->template getField<CHEMICALPOTENTIAL<T>>();
  auto& gradField = this->template getField<GRAD<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      T lambda = chemField.get(id);
      const Vector<T, LatSet::d>& grad = gradField.get(id);
      forceSF.SetField(id, grad * lambda);
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeForceP() {
  // Pressure gradient force: F_p = -∇p = -cs²·∇ρ = -cs²·(ρ_h-ρ_l)·∇φ
  auto& forceP = this->template getField<FORCE_P<T>>();
  auto& gradField = this->template getField<GRAD<T>>();
  auto& rhoL = this->template getField<RHO_L<T>>();
  auto& rhoH = this->template getField<RHO_H<T>>();
  auto& rhoField = this->template getField<RHO<T>>();
  auto& pressureField = this->template getField<PRESSURE<T>>();

  T rho_l = rhoL.get(0);
  T rho_h = rhoH.get(0);
  T coeff = (rho_h - rho_l);

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      forceP.SetField(id, -coeff * pressureField.get(id) * gradField.get(id) / rhoField.get(id));
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeForceBuoy() {
  // F_b = (ρ(φ) - ρ_h) · g
  // ρ(φ) = ρ_l + φ·(ρ_h-ρ_l),  g = (0, g_y) in lattice units
  auto& forceBuoy = this->template getField<FORCE_Buoy<T>>();
  auto& phiField = this->template getField<PHI<T>>();
  auto& rhoL = this->template getField<RHO_L<T>>();
  auto& rhoH = this->template getField<RHO_H<T>>();
  auto& gravField = this->template getField<GRAVITY<T>>();

  T rho_l = rhoL.get(0);
  T rho_h = rhoH.get(0);
  T g = gravField.get(0);
  T drho = rho_h - rho_l;

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();

  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      T phi = phiField.get(id);
      T rho = rho_l + phi * drho;
      forceBuoy.SetField(id, Vector<T, LatSet::d>{T{0}, (rho - rho_h) * g});
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeForceVisc() {
  // F_v = ∇·[η(φ)·(∇u + ∇u^T)]
  //
  // Two-pass isotropic FD:
  //   Pass 1: compute stress T = η·(∇u+∇u^T) at each interior cell
  //   Pass 2: compute F_v = ∇·T via isotropic FD divergence
  auto& forceVisc = this->template getField<FORCE_Visc<T>>();
  auto& velField = this->template getField<VELOCITY<T, 2>>();
  auto& phiField = this->template getField<PHI<T>>();
  auto& etaL = this->template getField<ETA_L<T>>();
  auto& etaH = this->template getField<ETA_H<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();
  std::size_t Ntotal = this->getN();

  T eta_l = etaL.get(0);
  T eta_h = etaH.get(0);
  T invCs2 = T{1} / LatSet::cs2;

  // Temporary storage for stress tensor components
  std::vector<T> Txx(Ntotal, T{0});
  std::vector<T> Txy(Ntotal, T{0});
  std::vector<T> Tyy(Ntotal, T{0});

  // Pass 1: T = η(φ)·(∇u + ∇u^T)  [isotropic FD for ∇u]
  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;

      T phi = phiField.get(id);
      T eta = eta_l + phi * (eta_h - eta_l);

      T duxdx = T{0}, duxdy = T{0};
      T duydx = T{0}, duydy = T{0};

      for (unsigned int k = 1; k < LatSet::q; ++k) {
        std::size_t nid = id + this->Delta_Index[k];
        const auto& un = velField.get(nid);
        T wk = latset::w<LatSet>(k);
        const auto& ck = latset::c<LatSet>(k);
        duxdx += wk * ck[0] * un[0];
        duxdy += wk * ck[1] * un[0];
        duydx += wk * ck[0] * un[1];
        duydy += wk * ck[1] * un[1];
      }

      Txx[id] = eta * T{2} * duxdx * invCs2;
      Txy[id] = eta * (duxdy + duydx) * invCs2;
      Tyy[id] = eta * T{2} * duydy * invCs2;
    }
  }

  // Pass 2: F_v = ∇·T  [isotropic FD for tensor divergence]
  for (int j = overlap + 1; j < Ny - overlap - 1; ++j) {
    for (int i = overlap + 1; i < Nx - overlap - 1; ++i) {
      std::size_t id = i + j * Nx;

      T Fx = T{0}, Fy = T{0};
      for (unsigned int k = 1; k < LatSet::q; ++k) {
        std::size_t nid = id + this->Delta_Index[k];
        T wk = latset::w<LatSet>(k);
        const auto& ck = latset::c<LatSet>(k);
        Fx += wk * (ck[0] * Txx[nid] + ck[1] * Txy[nid]);
        Fy += wk * (ck[0] * Txy[nid] + ck[1] * Tyy[nid]);
      }
      Fx *= invCs2;
      Fy *= invCs2;

      forceVisc.SetField(id, Vector<T, LatSet::d>{Fx, Fy});
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::computeForceTotal() {
  // F_total = F_s + F_g + F_b + F_v → NS FORCE field
  auto& forceTotal = this->template getField<FORCE<T, 2>>();
  auto& forceSF = this->template getField<FORCE_SF<T>>();
  auto& forceBuoy = this->template getField<FORCE_Buoy<T>>();
  auto& forceVisc = this->template getField<FORCE_Visc<T>>();
  auto& forceP = this->template getField<FORCE_P<T>>();

  int overlap = this->getOverlap();
  int Nx = this->getNx();
  int Ny = this->getNy();
  for (int j = overlap; j < Ny - overlap; ++j) {
    for (int i = overlap; i < Nx - overlap; ++i) {
      std::size_t id = i + j * Nx;
      forceTotal.SetField(id, forceSF.get(id)
               + forceBuoy.get(id) + forceVisc.get(id) + forceP.get(id));
    }
  }
}

template <typename T, typename LatSet>
void FF2DBlock<T, LatSet>::apply() {
  computeGradient();
  computeLaplacian();
  computeChemicalPotential();
  computeACSource();
  computeForceSF();
  computeForceP();
  computeForceBuoy();
  computeForceVisc();
  computeForceTotal();
}

// ===================================================================
// FF2DManager implementations
// ===================================================================

template <typename T, typename LatSet>
template <typename INITVALUEPACK, typename... FIELDPTRTYPES>
FF2DManager<T, LatSet>::FF2DManager(BlockGeometry<T, 2>& blockgeo,
                                    PhaseFieldConverter<T>& conv,
                                    INITVALUEPACK& initvalues,
                                    FIELDPTRTYPES*... fieldptrs)
    : BlockLatticeManagerBase<T, LatSet, FIELDPACK<T>>(blockgeo, initvalues,
                                                        fieldptrs...),
      _conv(conv) {
  for (int i = 0; i < this->BlockGeo.getBlockNum(); ++i) {
    _blocks.emplace_back(
        this->BlockGeo.getBlock(i), _conv,
        ExtractFieldPtrs<T, LatSet, FIELDPACK<T>>::getFieldPtrTuple(
            i, this->Fields, this->FieldPtrs));
  }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::Init(BlockGeometryHelper<T, 2>& GeoHelper) {
  // Fields already created by base constructor with init values.
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeGradient() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) {
    block.computeGradient();
  }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeLaplacian() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) {
    block.computeLaplacian();
  }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeChemicalPotential() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) {
    block.computeChemicalPotential();
  }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeACSource() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) {
    block.computeACSource();
  }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeForceSF() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) { block.computeForceSF(); }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeForceP() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) { block.computeForceP(); }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeForceBuoy() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) { block.computeForceBuoy(); }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeForceVisc() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) { block.computeForceVisc(); }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::computeForceTotal() {
#pragma omp parallel for num_threads(Thread_Num)
  for (auto& block : _blocks) { block.computeForceTotal(); }
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::apply() {
  computeGradient();
  computeLaplacian();
  computeChemicalPotential();
  computeACSource();
  computeForceSF();
  computeForceP();
  computeForceBuoy();
  computeForceVisc();
  computeForceTotal();
}

template <typename T, typename LatSet>
void FF2DManager<T, LatSet>::Communicate() {
  this->template getField<NORMAL<T>>().NormalCommunicate();
  this->template getField<LAPLACIAN<T>>().NormalCommunicate();
}

}  // namespace ff
