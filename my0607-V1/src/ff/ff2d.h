#pragma once

#include <array>
#include "data_struct/Vector.h"
#include "lbm/lattice_set.h"
#include "utils/alias.h"

namespace ff {

/// 相场源项: F_α = ω_α·(e_α·n) · 4φ(1-φ)/W
template <typename CELLTYPE>
struct PhaseFieldSource {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell, std::array<T, LatSet::q>& Fi) {
    T W = cell.template get<INTERFACEWIDTH<T>>();
    T phi = cell.template get<GenericRho>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL<T, LatSet::d>>();

    const T coeff = T{4} * phi * (T{1} - phi) / W;
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = latset::c<LatSet>(i) * n;
      Fi[i] = latset::w<LatSet>(i) * en * coeff;
    }
  }
};

/// 梯度算子: ∇φ = (1/c_s²) Σ_i w_i c_i φ(x+c_i)
template <typename CELLTYPE>
struct PhaseFieldGrad {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    Vector<T, LatSet::d> grad_phi;
    grad_phi[0] = T{0};
    grad_phi[1] = T{0};

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T phi_i = cell.getNeighbor(i).template get<GenericRho>();
      T wi = latset::w<LatSet>(i);
      const auto& ci = latset::c<LatSet>(i);
      grad_phi[0] += wi * phi_i * ci[0];
      grad_phi[1] += wi * phi_i * ci[1];
    }

    grad_phi[0] /= LatSet::cs2;
    grad_phi[1] /= LatSet::cs2;

    cell.template get<GRAD<T, LatSet::d>>() = grad_phi;
  }
};

/// 法向量算子: n = ∇φ/|∇φ|  (远场低于阈值时置零)
template <typename CELLTYPE>
struct PhaseFieldNormal {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    Vector<T, LatSet::d> n;
    n[0] = T{0};
    n[1] = T{0};
    T delta = T{0.005};
    const Vector<T, LatSet::d>& grad_phi = cell.template get<GRAD<T, LatSet::d>>();
    T grad_norm = grad_phi.getnorm();
    if (grad_norm > delta) {
      n[0] = grad_phi[0] / grad_norm;
      n[1] = grad_phi[1] / grad_norm;
    }
    cell.template get<NORMAL<T, LatSet::d>>() = n;
  }
};

template <typename CELLTYPE>
struct NSFieldSource {
  using CELL = CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

 __any__ static void apply(CELL& cell, std::array<T, LatSet::q>& Gi) {
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& F = cell.template get<FORCE<T, LatSet::d>>();
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T wi = latset::w<LatSet>(i);
      const auto& ci = latset::c<LatSet>(i);
      T udote = u * ci;
      T Fdote = F * ci;
      T udotF = u * F;
      T edotF = ci * F;
      Gi[i] = wi * (edotF * LatSet::InvCs2 + ((udote * Fdote) - udotF / LatSet::InvCs2) * LatSet::InvCs4);
    }


 };
  
};
                    
  template <typename T, typename LATTICE>       
  void BroadcastAllParams(LATTICE& lattice, T& interfaceWidth, T& rho_l, T& rho_h, T&eta_l, T& eta_h, T& gravity, T& Beta, T& Kappa) {                                   
  #ifdef MPI_ENABLED                            
    mpi().bCast(interfaceWidth, 0);             
    mpi().bCast(rho_l, 0);                      
    mpi().bCast(rho_h, 0);                      
    mpi().bCast(eta_l, 0);                      
    mpi().bCast(eta_h, 0);                      
    mpi().bCast(gravity, 0);                    
    mpi().bCast(Beta, 0);                       
    mpi().bCast(Kappa, 0);                      
  #endif                                        
    lattice.template getField<INTERFACEWIDTH<T>>().InitValue(interfaceWidth);                 
    lattice.template getField<RHO_L<T>>().InitValue(rho_l);        
    lattice.template getField<RHO_H<T>>().InitValue(rho_h);        
    lattice.template getField<ETA_L<T>>().InitValue(eta_l);        
    lattice.template getField<ETA_H<T>>().InitValue(eta_h);        
    lattice.template getField<BETA<T>>().InitValue(Beta);          
    lattice.template  getField<KAPPA<T>>().InitValue(Kappa);        
  }

}  // namespace ff
