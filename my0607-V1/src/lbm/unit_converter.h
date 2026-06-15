/*Lattice Boltzmann Unit Converter*/
#pragma once

#include <iostream>
#include <vector>

#include "data_struct/Vector.h"

// Lattice_deltaX = 1   ->  Conv_L = deltaX
// Lattice_deltaT = 1   ->  Conv_Time = deltaT
// Lattice_rho = 1   ->  Conv_rho = rho
template <typename T>
struct AbstractConverter {


  virtual T getLattice_RT() const = 0;
  virtual T getOMEGA() const = 0;
  virtual T getLatticeRho(T rho_phys) const = 0;
  virtual T getLatRhoInit() const = 0;
  virtual T getPhysRho(T Lattice_rho) const = 0;
  virtual T getLattice_gbeta() const { return 0; }
  virtual T getLatticeU(T U_phys) const { return 0; }
  virtual T getPhysU(T Lattice_U) const { return 0; }
};

template <typename T>
struct BaseConverter final : public AbstractConverter<T> {
  T deltaX;   // mm          // spacing between two lattice cells
  T deltaT;   // s           // time step
  T charL;    // mm          // i.e. simulation domain size
  T charU;    // mm / s      // maximal or expected velocity during simulation
  T charP;    // g / (mm s^2) = Pa
  T VisKine;  // mm^2 / s  // kinematic viscosity
  T rho;      // g / mm^3    // density
  T Re;       // Reynolds number, Re = U0*L0/nu = charU*charL/VisKine
  /*-----------*/
  T Conv_L;        // mm
  T Conv_Time;     // s
  T Conv_U;        // mm / s
  T Conv_rho;      // g / mm^3
  T Conv_Mass;     // g
  T Conv_VisKine;  // mm^2 / s
  T Conv_Force;    // g mm / s^2
  T Conv_Acc;      // mm / s^2
  T Conv_P;        // g / (mm s^2) = Pa
  /*-----------*/
  T Lattice_charL;    // lattice domain size = Ni
  T Lattice_charU;    // char lattice velocity
  T Lattice_VisKine;  // lattice kinematic viscosity
  T Lattice_Re;       // grid Reynolds number
  T Lattice_g;        // lattice gravity                 9.8 m / s^2 = 9800 mm / s^2

  T Lattice_RT;  // relaxation time
  T OMEGA;

  T cs2;

  BaseConverter(T cs2_) : cs2(cs2_) {}
  T getLattice_RT() const override { return Lattice_RT; }
  T getOMEGA() const override { return OMEGA; }
  T getLatticeRho(T rho_phys) const override { return rho_phys / Conv_rho; }
  T getLatRhoInit() const override { return T(1); }
  // = getPhysStrainRate()
  T getLatTime(T Phys_Time) const { return Phys_Time / deltaT; }
  // = getPhysTime()
  T getLatStrainRate(T Phys_StrainRate) const { return Phys_StrainRate * deltaT; }
  T getLatVisKine(T Phys_VisKine) const { return Phys_VisKine / Conv_VisKine; }

  // U method
  T getLatticeU(T U_phys) const override { return U_phys / Conv_U; }
  T getPhysU(T Lattice_U) const override { return Lattice_U * Conv_U; }
  template <unsigned int D>
  Vector<T, D> getLatticeU(const Vector<T, D> &U_phys) {
    return U_phys / Conv_U;
  }
  template <unsigned int D>
  Vector<T, D> getPhysU(const Vector<T, D> &Lattice_U) {
    return Lattice_U * Conv_U;
  }

  T getPhysRho(T Lattice_rho) const override { return Lattice_rho * Conv_rho; }
  // = getLatStrainRate()
  T getPhysTime(T Lattice_Time) { return Lattice_Time * deltaT; }
  // = getLatTime()
  T getPhysStrainRate(T Lattice_StrainRate) { return Lattice_StrainRate / deltaT; }
  T getPhysVisKine(T Lattice_VisKine) { return Lattice_VisKine * Conv_VisKine; }
  /*--------------------Basic Converters--------------------*/
  void Converter(T deltaX_, T deltaT_, T rho_, T charL_, T charU_, T VisKine_,
                 T charP_ = T(0)) {
    /*------------phys param-------------*/
    deltaX = deltaX_;
    deltaT = deltaT_;
    charL = charL_;
    charU = charU_;
    VisKine = VisKine_;
    rho = rho_;
    Re = charU_ * charL_ / VisKine_;
    /*----------conversion factors--------------*/
    Conv_L = deltaX;
    Conv_Time = deltaT;
    Conv_rho = rho;
    Conv_U = Conv_L / Conv_Time;                      // dx / dt
    Conv_Mass = Conv_rho * Conv_L * Conv_L * Conv_L;  // m = rho * dx^3
    Conv_VisKine = Conv_L * Conv_L / Conv_Time;       // nu = dx^2 / dt
    Conv_Force = Conv_Mass * Conv_L / Conv_Time /
                 Conv_Time;                     // F = m * dx / dt^2 = rho * dx^4 / dt^2
    Conv_Acc = Conv_L / Conv_Time / Conv_Time;  // a = dx / dt^2
    Conv_P = Conv_Force / Conv_L / Conv_L;      // P = F / dx^2 = rho * dx^2 / dt^2
    /*-----------lattice param------------*/
    Lattice_charU = charU / Conv_U;
    Lattice_charL = charL / Conv_L;
    Lattice_VisKine = VisKine / Conv_VisKine;
    Lattice_Re = Lattice_charU / Lattice_VisKine;
    Lattice_g = T(9810) / Conv_Acc;

    OMEGA = T(1) / Lattice_RT;
  }

  /*-------------------- Converters--------------------*/
  void ConvertFromRT(T deltaX_, T LatRT_, T rho_, T charL_, T charU_, T VisKine_) {
    Lattice_RT = LatRT_;
    Converter(deltaX_, (LatRT_ - T(0.5)) * cs2 * deltaX_ * deltaX_ / VisKine_, rho_,
              charL_, charU_, VisKine_);
  }
  void ConvertFromTimeStep(T deltaX_, T deltaT_, T rho_, T charL_, T charU_, T VisKine_) {
    Lattice_RT = T(0.5) + deltaT_ * VisKine_ / (cs2 * deltaX_ * deltaX_);
    Converter(deltaX_, deltaT_, rho_, charL_, charU_, VisKine_);
  }
  // deltaX, deltaT, rho = 1
  void SimplifiedConvertFromViscosity(int Ni_, T charU_, T VisKine_) {
    Lattice_RT = T(0.5) + VisKine_ / cs2;
    Converter(T(1), T(1), T(1), static_cast<T>(Ni_), charU_, VisKine_);
  }
  // deltaX, deltaT, rho = 1
  void SimplifiedConverterFromRT(int Ni_, T charU_, T LatRT_) {
    Lattice_RT = LatRT_;
    Converter(T(1), T(1), T(1), static_cast<T>(Ni_), charU_, (LatRT_ - T(0.5)) * cs2);
  }

  void check(int &check_status);
};



template <typename T>
class UnitConvManager {
  /*Unit converter for LB*/
  // Unit phys = Unit LB * Conversionfactor
 public:
  BaseConverter<T> *BaseConv = nullptr;
  std::vector<AbstractConverter<T> *> ConvList;

  UnitConvManager(AbstractConverter<T> *convlist) {
    if (convlist) ConvList.push_back(convlist);
  }
  UnitConvManager(BaseConverter<T> *bc) : BaseConv(bc) {
    if (bc) ConvList.push_back(bc);
  }

  void Check_and_Print();
};

// refine converter
template <typename T>
struct RefineConverter {
  // physical relaxation time(RT) reamains unchanged between different refinement levels
  // (Lat_RTC - 0.5) * deltaTC = (LatRTF - 0.5) * deltaTF
  // Lat_RTC = 0.5*Lat_RTF + 0.25
  // Lat_RTF = 2 * Lat_RTC - 0.5
  // OmegaC = 1 / (0.5 / OmegaF + 0.25)
  // OmegaF = 1 / (2 / OmegaC - 0.5)
  static inline T getOmegaF(T omegaC) { return T(1) / (T(2) / omegaC - T(0.5)); }
  static inline T getOmegaC(T omegaF) { return T(1) / (T(0.5) / omegaF + T(0.25)); }
  // recursively get OmegaF
  static T getOmegaF(T omegaC, std::uint8_t level) {
    if (level == std::uint8_t(0)) return omegaC;
    if (level == std::uint8_t(1)) return getOmegaF(omegaC);
    return getOmegaF(getOmegaF(omegaC, level - 1));
  }
  // recursively get OmegaC
  static T getOmegaC(T omegaF, std::uint8_t level) {
    if (level == std::uint8_t(0)) return omegaF;
    if (level == std::uint8_t(1)) return getOmegaC(omegaF);
    return getOmegaC(getOmegaC(omegaF, level - 1));
  }

  // continuous distribution function remains unchanged between different refinement
  // levels gC + (0.5*1/Lat_RTC)(geq - gC) = gF + (0.5*1/Lat_RTF)(geq - gF) gC = gF +
  // (0.5*1/Lat_RTF)(geq - gF) = gF + (0.5*OmegaF)(gF - geq) gF = gC -
  // (0.25*1/Lat_RTC)(geq - gC) = gC - (0.25*OmegaC)(gC - geq)
  static inline T getPopF(T popC, T popeq, T omegaC) {
    return popC - T(0.25) * omegaC * (popC - popeq);
  }
  template <unsigned int q>
  static inline void computePopF(std::array<T *, q> popC, const std::array<T, q> &feq,
                                 T omegaC) {
    for (unsigned int i = 0; i < q; ++i) {
      *(popC[i]) -= T(0.25) * omegaC * (*(popC[i]) - feq[i]);
    }
  }

  static inline T getPopC(T popF, T popeq, T omegaF) {
    return popF + T(0.5) * omegaF * (popF - popeq);
  }
  template <unsigned int q>
  static inline void computePopC(std::array<T *, q> popF, const std::array<T, q> &feq,
                                 T omegaF) {
    for (unsigned int i = 0; i < q; ++i) {
      *(popF[i]) += T(0.5) * omegaF * (*(popF[i]) - feq[i]);
    }
  }

  // get lattice g beta
  static inline T getLattice_gbetaF(T gbeta0, std::uint8_t level) {
    return gbeta0 / std::pow(T(2), level);
  }
};