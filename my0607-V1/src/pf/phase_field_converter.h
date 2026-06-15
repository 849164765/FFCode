// phase_field_converter.h
// Unit converter for Allen-Cahn phase field parameters
//
// Converts interface width, mobility, and surface tension to lattice
// parameters: Kappa = 3*sigma*W/2, Beta = 12*sigma/W.
//
// Usage:
//   BaseConverter<T> bc(LatSet::cs2);
//   bc.SimplifiedConvertFromViscosity(Ni, charU, VisKine);
//   PhaseFieldConverter<T> phiConv(bc);
//   phiConv.Converter(interfaceWidth, mobility, surfaceTension);

#pragma once

#include "lbm/unit_converter.h"

template <typename T>
struct PhaseFieldConverter : public AbstractConverter<T> {
  BaseConverter<T>& BaseConv;

  T Interface_Width{};
  T Mobility{};
  T SurfaceTension{};
  T Kappa{};
  T Beta{};
  T Lattice_RT_phi{};
  T OMEGA_phi{};

  explicit PhaseFieldConverter(BaseConverter<T>& bc) : BaseConv(bc) {}

  void Converter(T W, T M, T sigma) {
    Interface_Width = W;
    Mobility = M;
    SurfaceTension = sigma;
    Kappa = T(3) * sigma * Interface_Width * T(0.5);
    Beta = T(12) * sigma / Interface_Width;
    Lattice_RT_phi = T(0.5) + Mobility * BaseConv.cs2;
    OMEGA_phi = T(1) / Lattice_RT_phi;
  }

  T getLattice_RT() const override { return Lattice_RT_phi; }
  T getOMEGA() const override { return OMEGA_phi; }
  T getNSLatticeRT() const { return BaseConv.getLattice_RT(); }
  T getNSOMEGA() const { return BaseConv.getOMEGA(); }
  T getLatticeRho(T rho) const override { return BaseConv.getLatticeRho(rho); }
  T getLatRhoInit() const override { return BaseConv.getLatRhoInit(); }
  T getPhysRho(T rho) const override { return BaseConv.getPhysRho(rho); }
  T getLatticeU(T u) const override { return BaseConv.getLatticeU(u); }
  T getPhysU(T u) const override { return BaseConv.getPhysU(u); }

  T getOmegaPhi() const { return OMEGA_phi; }
  T getTauPhi() const { return Lattice_RT_phi; }
};
