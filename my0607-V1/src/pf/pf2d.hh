// pf2d.hh
// Phase field global parameters (set by user before simulation)

#pragma once

namespace pf {

// Global parameters for the Allen-Cahn phase field model.
// Must be set before simulation starts.
template <typename T>
T Kappa{};       // gradient coefficient: 3*sigma*W/2

template <typename T>
T Beta{};        // double-well coefficient: 12*sigma/W

template <typename T>
T Omega_phi{};   // phase field relaxation frequency: 1/(0.5 + Mobility/cs2)

}  // namespace pf
