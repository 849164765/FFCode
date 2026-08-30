#pragma once
// Mag2D.hh — 完整 Maxwell 磁应力形式的 2D 磁力算子
//
// 论文 Eq.(8) 只给出恒 χ 时的体积 Kelvin 力 F=(μ0χ/2)∇|H|²。
// 对两相界面（χ 随 φ 变化），完整的 Maxwell 应力散度为：
//
//   F_m = (χ/2)·∇(H²) - (1/2)·H²·∇χ
//
// 其中：
//   ∇χ = (χ_h - χ_l)·∇φ
//   ∇(H²) 用 D2Q5 各向同性梯度由邻居 H² = Hx²+Hy² 计算。
//
// 这一形式同时包含：
//   - 体积 Kelvin 项：0.5·χ·∇(H²)
//   - 界面磁压项：-0.5·H²·∇χ
// 是“另一个 Maxwell 磁应力形式”的完整实现。

#include "data_struct/cell.h"

namespace mfield {

template <typename PFCELL, typename MFCELL, typename NSCELL>
struct MagMagneticForce2D {
  using T = typename PFCELL::FloatType;
  using LatSet = typename MFCELL::LatticeSet;  // D2Q5

  __any__ static void apply(PFCELL& pf_cell, MFCELL& mf_cell, NSCELL& ns_cell) {
    T chi_l = mf_cell.template get<CHI_L<T>>();
    T chi_h = mf_cell.template get<CHI_H<T>>();
    T Hx = mf_cell.template get<HX<T>>();
    T Hy = mf_cell.template get<HY<T>>();
    T H2 = Hx*Hx + Hy*Hy;
    T Hmag = std::sqrt(H2);
    T phi = pf_cell.template get<typename PFCELL::GenericRho>();
    T chi = chi_l + phi*(chi_h - chi_l);

    // ∇(H²) via D2Q5 isotropic gradient
    Vector<T,2> gH2{0,0};
    for (unsigned int k=1;k<LatSet::q;++k) {
      auto nc = mf_cell.getNeighbor(k);
      T hx = nc.template get<HX<T>>();
      T hy = nc.template get<HY<T>>();
      T h2n = hx*hx + hy*hy;
      gH2[0] += latset::w<LatSet>(k) * latset::c<LatSet>(k)[0] * h2n;
      gH2[1] += latset::w<LatSet>(k) * latset::c<LatSet>(k)[1] * h2n;
    }
    gH2[0] /= LatSet::cs2;
    gH2[1] /= LatSet::cs2;

    // ∇χ = (χ_h - χ_l)·∇φ
    const auto& gphi = pf_cell.template get<GRAD<T,2>>();
    T dchi = chi_h - chi_l;
    Vector<T,2> gchi{dchi*gphi[0], dchi*gphi[1]};

    // 完整 Maxwell 磁应力体积力：
    //   F = 0.5·χ·∇H²  -  0.5·H²·∇χ
    T Fx = T{0.5}*chi*gH2[0] - T{0.5}*H2*gchi[0];
    T Fy = T{0.5}*chi*gH2[1] - T{0.5}*H2*gchi[1];

    auto& ns_force = ns_cell.template get<FORCE<T,2>>();
    ns_force[0] += Fx;
    ns_force[1] += Fy;
  }
};

} // namespace mfield
