// collisionMRT.h

#pragma once

#include "lbm/moment.ur.h"

namespace mrtdata {
// mrt transformation matrix: mapping from population to moment space
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> M[Q][Q] = {};

// Inverse of M
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> InvM[Q][Q] = {};

// Standard d'Humières 2002 D3Q19 MRT basis (fully orthogonal, filled for D=3,Q=19)
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> M_std[Q][Q] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> InvM_std[Q][Q] = {};

// relaxation times
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s[Q] = {};

// relaxation times
template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s2[Q] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ int shearIndexes = {};

// relevant indexes of r. t. for shear viscosity
template <unsigned int D, unsigned int Q>
__constexpr__ int shearViscIndexes[shearIndexes<D,Q>] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s_phase_source[Q] = {};

template <unsigned int D, unsigned int Q>
__constexpr__ Fraction<> s1[Q] = {};

// mrt transformation matrix
template <>
__constexpr__ Fraction<> M<2,9>[9][9] = {
{ 1,  1,  1,  1,  1,  1,  1,  1,  1},
{-4, -1, -1, -1, -1,  2,  2,  2,  2},
{ 4, -2, -2, -2, -2,  1,  1,  1,  1},
{ 0,  1, -1,  0,  0,  1, -1,  1, -1},
{ 0, -2,  2,  0,  0,  1, -1,  1, -1},
{ 0,  0,  0,  1, -1,  1, -1, -1,  1},
{ 0,  0,  0, -2,  2,  1, -1, -1,  1},
{ 0,  1,  1, -1, -1,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  1,  1, -1, -1}
};

template <>
__constexpr__ Fraction<> InvM<2,9>[9][9] = {
  {{1, 9}, {-1,  9}, { 1,  9},       0,        0,       0,        0,       0,       0},
  {{1, 9}, {-1, 36}, {-1, 18}, { 1, 6}, {-1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18}, {-1, 6}, { 1,  6},       0,        0, { 1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, { 1, 6}, {-1,  6}, {-1, 4},       0},
  {{1, 9}, {-1, 36}, {-1, 18},       0,        0, {-1, 6}, { 1,  6}, {-1, 4},       0},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, { 1, 6}, { 1, 12},       0, { 1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, {-1, 6}, {-1, 12},       0, { 1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, { 1, 6}, { 1, 12}, {-1, 6}, {-1, 12},       0, {-1, 4}},
  {{1, 9}, { 1, 18}, { 1, 36}, {-1, 6}, {-1, 12}, { 1, 6}, { 1, 12},       0, {-1, 4}}
};

template <>
__constexpr__ Fraction<> s<2,9>[9] = {
0, {11, 10}, {11, 10}, 0, {11, 10}, 0, {11, 10}, 0, 0
};

template <>
__constexpr__ Fraction<> s_phase_source<2,9>[9] = {
  0, {11, 10}, {11, 10}, 0, {11, 10}, 0, {11, 10}, 1, 1
};

template <>
__constexpr__ Fraction<> s1<2,9>[9] = {
  1, {14, 10}, {14, 10}, 1, {12, 10}, 1, {12, 10}, 1, 1
};

template <>
__constexpr__ int shearIndexes<2,9> = 2;

template <>
__constexpr__ int shearViscIndexes<2,9>[shearIndexes<2,9>] = { 7, 8};

// ===================================================================
// D2Q5 MRT — magnetic field potential solver (M_5 from Guo et al. 2025, Eq. 41)
// Moments: m = (ψ, J_x, J_y, e, p_xx)
// ψ: conserved (zeroth moment)
// J_x, J_y: flux moments (relaxed with s_h)
// e, p_xx: higher-order moments (relaxed with s=1)
// InvM_5 = true inverse of M_5 (M is orthogonal: M^{-1}[i][j] = M[j][i] / ||M_j||²)
// Row norms: ||r0||²=5, ||r1||²=2, ||r2||²=2, ||r3||²=20, ||r4||²=4
// ===================================================================

template <>
__constexpr__ Fraction<> M<2,5>[5][5] = {
  { 1,  1,  1,  1,  1},
  { 0,  1, -1,  0,  0},
  { 0,  0,  0,  1, -1},
  { 4, -1, -1, -1, -1},
  { 0,  1,  1, -1, -1}
};

template <>
__constexpr__ Fraction<> InvM<2,5>[5][5] = {
  {{1, 5},       0,       0,   {1, 5},       0},
  {{1, 5},   {1, 2},       0,  {-1, 20},  {1, 4}},
  {{1, 5},  {-1, 2},       0,  {-1, 20},  {1, 4}},
  {{1, 5},       0,   {1, 2},  {-1, 20}, {-1, 4}},
  {{1, 5},       0,  {-1, 2},  {-1, 20}, {-1, 4}}
};

template <>
__constexpr__ int shearIndexes<2,5> = 0;

template <>
__constexpr__ int shearViscIndexes<2,5>[0] = {};

template <>
__constexpr__ Fraction<> M<3,19>[19][19] = {
{ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{-1,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{ 1, -2, -2, -2, -2, -2, -2,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
{ 0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0},
{ 0, -2,  2,  0,  0,  0,  0,  1, -1,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0},
{ 0,  0,  0,  1, -1,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1},
{ 0,  0,  0, -2,  2,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1},
{ 0,  0,  0,  0,  0,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1},
{ 0,  0,  0,  0,  0, -2,  2,  0,  0,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1},
{ 0,  2,  2, -1, -1, -1, -1,  1,  1,  1,  1, -2, -2,  1,  1,  1,  1, -2, -2},
{ 0, -2, -2,  1,  1,  1,  1,  1,  1,  1,  1, -2, -2,  1,  1,  1,  1, -2, -2},
{ 0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1,  0,  0,  1,  1, -1, -1,  0,  0},
{ 0,  0,  0, -1, -1,  1,  1,  1,  1, -1, -1,  0,  0,  1,  1, -1, -1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0, -1, -1},
{ 0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0,  1, -1, -1,  1,  0,  0},
{ 0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1, -1,  1,  0,  0, -1,  1},
{ 0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  0,  0, -1,  1,  1, -1}
};

template <>
__constexpr__ Fraction<> InvM<3,19>[19][19] = {
  {{1, 3}, {-1, 2}, {1, 6},       0,        0,       0,       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18}, {1, 6}, {-1, 6},       0,       0,        0,        0,  {1,12}, {-1,12},        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18}, {-1,6},  {1, 6},       0,       0,        0,        0,  {1,12}, {-1,12},        0,        0,        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0, {1, 6}, {-1,6},        0,        0, {-1,24},  {1,24},  {1, 8}, {-1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0, {-1,6},  {1,6},        0,        0, {-1,24},  {1,24},  {1, 8}, {-1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0,      0,       0,  {1, 6}, {-1, 6}, {-1,24},  {1,24}, {-1, 8},  {1, 8},        0,        0,        0,        0,        0,        0},
  {{1,18},       0, {-1,18},      0,        0,      0,       0, {-1, 6},  {1, 6}, {-1,24},  {1,24}, {-1, 8},  {1, 8},        0,        0,        0,        0,        0,        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24}, {1,12}, {1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16},  {1, 4},        0,        0,  {1, 8},  {1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24}, {-1,12},{-1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16},  {1, 4},        0,        0, {-1, 8}, {-1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24},      0,       0,  {1,12},  {1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0,  {1, 4},        0, {-1, 8},        0,  {1, 8}},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24},      0,       0, {-1,12}, {-1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0,  {1, 4},        0,  {1, 8},        0, {-1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {1,12}, {1,24},  {1,12},  {1,24}, {-1,24}, {-1,24},        0,        0,        0,        0,  {1, 4},        0, {-1, 8}, {-1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {-1,12},{-1,24}, {-1,12}, {-1,24}, {-1,24}, {-1,24},        0,        0,        0,        0,  {1, 4},        0,  {1, 8},  {1, 8}},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24}, {-1,12},{-1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16}, {-1, 4},        0,        0,  {1, 8}, {-1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24}, {1,12}, {1,24},        0,        0,  {1,48},  {1,48},  {1,16},  {1,16}, {-1, 4},        0,        0, {-1, 8},  {1, 8},        0},
  {{1,36}, {1,24}, {1,72}, {1,12}, {1,24},      0,       0, {-1,12}, {-1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0, {-1, 4},        0, {-1, 8},        0, {-1, 8}},
  {{1,36}, {1,24}, {1,72}, {-1,12},{-1,24},      0,       0,  {1,12},  {1,24},  {1,48},  {1,48}, {-1,16}, {-1,16},        0, {-1, 4},        0,  {1, 8},        0,  {1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {1,12}, {1,24}, {-1,12}, {-1,24}, {-1,24}, {-1,24},        0,        0,        0,        0, {-1, 4},        0, {-1, 8},  {1, 8}},
  {{1,36}, {1,24}, {1,72},      0,        0, {-1,12},{-1,24}, {1,12},  {1,24}, {-1,24}, {-1,24},        0,        0,        0,        0, {-1, 4},        0,  {1, 8}, {-1, 8}}
};

template <>
__constexpr__ Fraction<> s<3,19>[19] = {
0, {119, 100}, {14, 10}, 0, {12, 10}, 0, {12, 10}, 0, {12, 10}, 0, {14, 10}, 0, {14, 10}, 0, 0, 0, {99, 50}, {99, 50}, {99, 50}
};

template <>
__constexpr__ Fraction<> s1<3,19>[19] = {
1, {119, 100}, {14, 10}, 1, {12, 10}, 1, {12, 10}, 1, {12, 10}, 0, {14, 10}, 0, {14, 10}, 0, 0, 0, {99, 50}, {99, 50}, {99, 50}
};

template <>
__constexpr__ int shearIndexes<3,19> = 5;

template <>
__constexpr__ int shearViscIndexes<3,19>[shearIndexes<3,19>] = { 9, 11, 13, 14, 15};

// ===================================================================
// Standard D3Q19 MRT basis (d'Humières et al. 2002)
// Gram-Schmidt orthogonalized polynomials on the D3Q19 lattice.
// M_std is fully orthogonal: M_std * InvM_std = I.
// InvM_std = M_std^T * diag(1/||r_i||^2).
//
// Moment ordering (rows 0-18):
//   0: rho   (density)
//   1: e     (energy)
//   2: eps   (energy squared)
//   3: jx    (x-momentum)
//   4: qx    (x-energy flux)
//   5: jy    (y-momentum)
//   6: qy    (y-energy flux)
//   7: jz    (z-momentum)
//   8: qz    (z-energy flux)
//   9: 2cx²-cy²-cz²  (diagonal stress)
//  10: cy²-cz²       (diagonal stress diff)
//  11: cx·cy         (shear xy)
//  12: cy·cz         (shear yz)
//  13: cx·cz         (shear xz)
//  14: cx·(cy²-cz²)  (3rd order cubic x)
//  15: cy·(cz²-cx²)  (3rd order cubic y)
//  16: cz·(cx²-cy²)  (3rd order cubic z)
//  17: cx²·cy²       (4th order even)
//  18: cy²·cz²       (4th order even)
// ===================================================================

template <>
__constexpr__ Fraction<> M_std<3,19>[19][19] = {
    {   1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1},
    { -30,  -11,  -11,  -11,  -11,  -11,  -11,    8,    8,    8,    8,    8,    8,    8,    8,    8,    8,    8,    8},
    {  12,   -4,   -4,   -4,   -4,   -4,   -4,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1},
    {   0,    1,   -1,    0,    0,    0,    0,    1,   -1,    1,   -1,    0,    0,    1,   -1,    1,   -1,    0,    0},
    {   0,   -4,    4,    0,    0,    0,    0,    1,   -1,    1,   -1,    0,    0,    1,   -1,    1,   -1,    0,    0},
    {   0,    0,    0,    1,   -1,    0,    0,    1,   -1,    0,    0,    1,   -1,   -1,    1,    0,    0,    1,   -1},
    {   0,    0,    0,   -4,    4,    0,    0,    1,   -1,    0,    0,    1,   -1,   -1,    1,    0,    0,    1,   -1},
    {   0,    0,    0,    0,    0,    1,   -1,    0,    0,    1,   -1,    1,   -1,    0,    0,   -1,    1,   -1,    1},
    {   0,    0,    0,    0,    0,   -4,    4,    0,    0,    1,   -1,    1,   -1,    0,    0,   -1,    1,   -1,    1},
    {   0,    2,    2,   -1,   -1,   -1,   -1,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0},
    {   0,    0,    0,    1,    1,   -1,   -1,    1,    1,   -1,   -1,    0,    0,    1,    1,   -1,   -1,    0,    0},
    {   0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,   -1,   -1,    0,    0,    0,    0},
    {   0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,   -1,   -1},
    {   0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    0,    0,    0,    0,   -1,   -1,    0,    0},
    {   0,    0,    0,    0,    0,    0,    0,    1,   -1,   -1,    1,    0,    0,    1,   -1,   -1,    1,    0,    0},
    {   0,    0,    0,    0,    0,    0,    0,   -1,    1,    0,    0,    1,   -1,    1,   -1,    0,    0,    1,   -1},
    {   0,    0,    0,    0,    0,    0,    0,    0,    0,    1,   -1,   -1,    1,    0,    0,   -1,    1,    1,   -1},
    {   0,    0,    0,   -1,   -1,    1,    1,    1,    1,    0,    0,   -1,   -1,    1,    1,    0,    0,   -1,   -1},
    {   0,    0,    0,   -1,   -1,    1,    1,    0,    0,   -1,   -1,    1,    1,    0,    0,   -1,   -1,    1,    1}
};

template <>
__constexpr__ Fraction<> InvM_std<3,19>[19][19] = {
    {{   1,   19}, {  -5,  399}, {   1,   21},            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0,            0},
    {{   1,   19}, { -11, 2394}, {  -1,   63}, {   1,   10}, {  -1,   10},            0,            0,            0,            0, {   1,    6},            0,            0,            0,            0,            0,            0,            0,            0,            0},
    {{   1,   19}, { -11, 2394}, {  -1,   63}, {  -1,   10}, {   1,   10},            0,            0,            0,            0, {   1,    6},            0,            0,            0,            0,            0,            0,            0,            0,            0},
    {{   1,   19}, { -11, 2394}, {  -1,   63},            0,            0, {   1,   10}, {  -1,   10},            0,            0, {  -1,   12}, {   1,   12},            0,            0,            0,            0,            0,            0, {  -1,    4}, {  -1,    4}},
    {{   1,   19}, { -11, 2394}, {  -1,   63},            0,            0, {  -1,   10}, {   1,   10},            0,            0, {  -1,   12}, {   1,   12},            0,            0,            0,            0,            0,            0, {  -1,    4}, {  -1,    4}},
    {{   1,   19}, { -11, 2394}, {  -1,   63},            0,            0,            0,            0, {   1,   10}, {  -1,   10}, {  -1,   12}, {  -1,   12},            0,            0,            0,            0,            0,            0, {   1,    4}, {   1,    4}},
    {{   1,   19}, { -11, 2394}, {  -1,   63},            0,            0,            0,            0, {  -1,   10}, {   1,   10}, {  -1,   12}, {  -1,   12},            0,            0,            0,            0,            0,            0, {   1,    4}, {   1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {   1,   10}, {   1,   40}, {   1,   10}, {   1,   40},            0,            0,            0, {   1,   12}, {   1,    4},            0,            0, {   1,    8}, {  -1,    8},            0, {   1,    4},            0},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {  -1,   10}, {  -1,   40}, {  -1,   10}, {  -1,   40},            0,            0,            0, {   1,   12}, {   1,    4},            0,            0, {  -1,    8}, {   1,    8},            0, {   1,    4},            0},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {   1,   10}, {   1,   40},            0,            0, {   1,   10}, {   1,   40},            0, {  -1,   12},            0,            0, {   1,    4}, {  -1,    8},            0, {   1,    8},            0, {  -1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {  -1,   10}, {  -1,   40},            0,            0, {  -1,   10}, {  -1,   40},            0, {  -1,   12},            0,            0, {   1,    4}, {   1,    8},            0, {  -1,    8},            0, {  -1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252},            0,            0, {   1,   10}, {   1,   40}, {   1,   10}, {   1,   40},            0,            0,            0, {   1,    4},            0,            0, {   1,    8}, {  -1,    8}, {  -1,    4}, {   1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252},            0,            0, {  -1,   10}, {  -1,   40}, {  -1,   10}, {  -1,   40},            0,            0,            0, {   1,    4},            0,            0, {  -1,    8}, {   1,    8}, {  -1,    4}, {   1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {   1,   10}, {   1,   40}, {  -1,   10}, {  -1,   40},            0,            0,            0, {   1,   12}, {  -1,    4},            0,            0, {   1,    8}, {   1,    8},            0, {   1,    4},            0},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {  -1,   10}, {  -1,   40}, {   1,   10}, {   1,   40},            0,            0,            0, {   1,   12}, {  -1,    4},            0,            0, {  -1,    8}, {  -1,    8},            0, {   1,    4},            0},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {   1,   10}, {   1,   40},            0,            0, {  -1,   10}, {  -1,   40},            0, {  -1,   12},            0,            0, {  -1,    4}, {  -1,    8},            0, {  -1,    8},            0, {  -1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252}, {  -1,   10}, {  -1,   40},            0,            0, {   1,   10}, {   1,   40},            0, {  -1,   12},            0,            0, {  -1,    4}, {   1,    8},            0, {   1,    8},            0, {  -1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252},            0,            0, {   1,   10}, {   1,   40}, {  -1,   10}, {  -1,   40},            0,            0,            0, {  -1,    4},            0,            0, {   1,    8}, {   1,    8}, {  -1,    4}, {   1,    4}},
    {{   1,   19}, {   4, 1197}, {   1,  252},            0,            0, {  -1,   10}, {  -1,   40}, {   1,   10}, {   1,   40},            0,            0,            0, {  -1,    4},            0,            0, {  -1,    8}, {  -1,    8}, {  -1,    4}, {   1,    4}}
};

} // namespace mrtdata

namespace mrt {

template <typename LatSet>
constexpr typename LatSet::FloatType M(unsigned int i, unsigned int j) {
  return mrtdata::M<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType InvM(unsigned int i, unsigned int j) {
  return mrtdata::InvM<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType s(unsigned int i) {
  return mrtdata::s<LatSet::d,LatSet::q>[i].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType s1(unsigned int i) {
  return mrtdata::s1<LatSet::d,LatSet::q>[i].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr int shearIndexes() {
  return mrtdata::shearIndexes<LatSet::d,LatSet::q>;
}

template <typename LatSet>
constexpr int shearViscIndexes(unsigned int i) {
  return mrtdata::shearViscIndexes<LatSet::d,LatSet::q>[i];
}

// Standard d'Humières MRT accessors — used by FFViscoForce3DM.
template <typename LatSet>
constexpr typename LatSet::FloatType M_standard(unsigned int i, unsigned int j) {
  return mrtdata::M_std<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

template <typename LatSet>
constexpr typename LatSet::FloatType InvM_standard(unsigned int i, unsigned int j) {
  return mrtdata::InvM_std<LatSet::d,LatSet::q>[i][j].template operator()<typename LatSet::FloatType>();
}

} // namespace mrt

namespace collision {

// a typical MRT collision process with:
// macroscopic variables updated
// equilibrium distribution function calculated
template <typename CELL, bool WriteToField = false>
struct MRT_Feq_RhoU {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = equilibrium::SecondOrder<CELL>;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    const T omega = cell.getOmega();
    // relaxation time vector
    T rtvec[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      rtvec[i] = mrt::s<LatSet>(i);
    }
    for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i) {
      rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega;
    }
    // relaxation time matrix
    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    // update macroscopic variables
    T rho{};
    Vector<T, LatSet::d> u{};
    moment::template rhoU<CELL, WriteToField>::apply(cell, rho, u);

    // compute Momenta
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }
    // compute Equilibrium
    std::array<T, LatSet::q> feq{};
    equilibrium::SecondOrder<CELL>::apply(feq, rho, u);
    T momentaEq[LatSet::q] {};
    for(unsigned int i = 0; i < LatSet::q; ++i) {
      momentaEq[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
      }
    }
    // MRT collision
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * (momenta[j] - momentaEq[j]);
      }
      cell[i] -= momentaEqS;
    }
  }

};

// MRT source collision for phase-field equations.
//
// Template parameters:
//   UseCHRelaxation = false (default): interface sharpening mode.
//     j-moments use s=1 (full relaxation), suitable for NORMAL=∇φ/|∇φ|
//     (e.g., simpledrop2dMRT).
//   UseCHRelaxation = true: Cahn-Hilliard diffusion mode.
//     j-moments use s=omega_phi (matches BGKSource for correct
//     effective mobility), suitable for NORMAL=∇μ
//     (e.g., rti2dMRT, bubble).
template <typename EquilibriumScheme, typename NORMAL, bool WriteToField = false, bool UseCHRelaxation = false>
struct MRTSource {
  using CELL = typename EquilibriumScheme::CELLTYPE;
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = EquilibriumScheme;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    const Vector<T, LatSet::d>& u = cell.template get<VELOCITY<T, LatSet::d>>();
    const Vector<T, LatSet::d>& n = cell.template get<NORMAL>();

    T phi = cell.template get<GenericRho>();
    // Clamp phi to [0,1] with NaN/Inf guard before Allen-Cahn source term.
    // IEEE 754 makes NaN<x false, so isfinite must be checked first.
    if (!std::isfinite(phi) || phi < T{0}) phi = T{0};
    if (phi > T{1}) phi = T{1};
    std::array<T, LatSet::q> feq;
    EquilibriumScheme::apply(feq, phi, u);

    T omega_phi = cell.getOmega();
    T interfacewidth = cell.template get<INTERFACEWIDTH<T>>();

    T rtvec[LatSet::q] {};
    if constexpr (UseCHRelaxation) {
      // Cahn-Hilliard: Fortran-aligned relaxation vector
      if constexpr (std::is_same_v<LatSet, D2Q9<T>> || std::is_same_v<LatSet, D2Q5<T>>) {
        // Fortran Sf = {1.0, 1.1, 1.1, 1/tau, 1/tau, 1/tau, 1/tau, 1.2, 1.2}
        rtvec[0] = T{1};                                    // rho/phi (conserved)
        rtvec[1] = T{11./10.};                               // e
        rtvec[2] = T{11./10.};                               // epsilon
        rtvec[3] = omega_phi;                                // jx (1/tau)
        rtvec[4] = omega_phi;                                // qx (1/tau)
        if constexpr (LatSet::q > 5) {
          rtvec[5] = omega_phi;                              // jy (1/tau)
          rtvec[6] = omega_phi;                              // qy (1/tau)
          rtvec[7] = T{12./10.};                             // pxx (1.2)
          rtvec[8] = T{12./10.};                             // pxy (1.2)
        }
      } else if constexpr (std::is_same_v<LatSet, D3Q19<T>>) {
        // D3Q19 CH relaxation vector aligned with mrtdata::M<3,19> row order:
        // 0:rho, 1:e, 2:epsilon, 3:jx, 4:qx, 5:jy, 6:qy, 7:jz, 8:qz,
        // 9-18: stress/shear and higher-order moments (9,11,13,14,15 are shear modes)
        rtvec[0] = T{1};            // rho/phi (conserved)
        rtvec[1] = T{11./10.};       // e
        rtvec[2] = T{11./10.};       // epsilon
        rtvec[3] = omega_phi;        // jx (1/tau_phi)
        rtvec[4] = omega_phi;        // qx (1/tau_phi)
        rtvec[5] = omega_phi;        // jy (1/tau_phi)
        rtvec[6] = omega_phi;        // qy (1/tau_phi)
        rtvec[7] = omega_phi;        // jz (1/tau_phi)
        rtvec[8] = omega_phi;        // qz (1/tau_phi)
        rtvec[9] = T{12./10.};       // pxx
        rtvec[10] = T{12./10.};      // pww
        rtvec[11] = T{12./10.};      // pxy
        rtvec[12] = T{12./10.};      // pyz
        rtvec[13] = T{12./10.};      // pxz
        rtvec[14] = T{12./10.};      // higher-order stress
        rtvec[15] = T{12./10.};      // higher-order stress
        rtvec[16] = T{12./10.};      // higher-order ghost
        rtvec[17] = T{12./10.};      // higher-order ghost
        rtvec[18] = T{12./10.};      // higher-order ghost
      } else {
        // Generic fallback: use omega_phi for all moments
        for (unsigned int i = 0; i < LatSet::q; ++i) {
          rtvec[i] = omega_phi;
        }
      }
    } else {
      // Interface sharpening: s1 + shear override (original behaviour)
      for (unsigned int i = 0; i < LatSet::q; ++i) {
        rtvec[i] = mrt::s1<LatSet>(i);
      }
      for (int i = 0; i < mrt::shearIndexes<LatSet>(); ++i) {
        rtvec[mrt::shearViscIndexes<LatSet>(i)] = omega_phi;
      }
    }

    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momenta[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }

    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      momentaEq[i] = T{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];
      }
    }

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T momentaEqS{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momentaEqS += InvM_S[i][j] * (momenta[j] - momentaEq[j]);
      }
      cell[i] -= momentaEqS;
    }

    T A = T{4} * phi * (T{1} - phi) / interfacewidth;

    T S_body[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T en = T{};
      Vector<T, LatSet::d> c = latset::c<LatSet>(i);
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        en += c[d] * n[d];
      }
      S_body[i] = latset::w<LatSet>(i) * en * A;
    }

    T source_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j) {
      source_m[j] = T{};
      for (unsigned int k = 0; k < LatSet::q; ++k) {
        source_m[j] += mrt::M<LatSet>(j, k) * S_body[k];
      }
    }

    T source_m_corrected[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j) {
      source_m_corrected[j] = (T{1} - rtvec[j] / T{2}) * source_m[j];
    }

    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        source += mrt::InvM<LatSet>(i, j) * source_m_corrected[j];
      }
      cell[i] += source;
    }
  }

};

// MRT collision with Guo force scheme for Navier-Stokes
// replaces BGKForce<moment::forcerhoU, equilibrium::SecondOrder, force::Force>
template <typename CELL, typename ForceField>
struct MRTForce {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // Per-cell omega (spatially varying relaxation) for variable viscosity
    T omega{};
    if constexpr (cell.template hasField<OMEGA<T>>()) {
      omega = cell.template get<OMEGA<T>>();
    } else {
      omega = cell.getOmega();
    }
    const auto& F = cell.template get<ForceField>();

    // Compile-time check: incompressible (He-Luo) vs compressible LBM
    // DENSITY field present => He-Luo incompressible (Fortran-aligned)
    constexpr bool isIncompressible = CELL::template hasField<DENSITY<T>>();

    // 1. Compute zeroth moment p=Σf (pressure) and first moment Σ(c·f)
    //    Incompressible: u = Σ(c·f) (velocity, NOT divided by rho)
    //    Compressible:   u = Σ(c·f)/rho (momentum divided by rho)
    T zeroth{};  // p = Σf (pressure)
    Vector<T, LatSet::d> u{};
    for (unsigned int k = 0; k < LatSet::q; ++k) {
      zeroth += cell[k];
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        u[d] += latset::c<LatSet>(k)[d] * cell[k];
      }
    }

    // rho: incompressible reads from DENSITY field (phi-based, set by
    // FFRhoOmegaUpdate3D); compressible falls back to rho = Σf
    T rho{};
    if constexpr (isIncompressible) {
      rho = cell.template get<DENSITY<T>>();
    } else {
      rho = zeroth;
      for (unsigned int d = 0; d < LatSet::d; ++d) u[d] /= rho;
    }

    // 2. Half-force correction: u_c = u + F/(2*rho)
    Vector<T, LatSet::d> uc = u;
    for (unsigned int d = 0; d < LatSet::d; ++d) uc[d] += F[d] / (T{2} * rho);

    // Write macroscopic fields
    // Incompressible: do NOT write GenericRho (=DENSITY), which is set from phi
    // by FFRhoOmegaUpdate3D and must not be overwritten with Σf.
    // Compressible: write rho = Σf to GenericRho.
    if constexpr (!isIncompressible) {
      cell.template get<GenericRho>() = rho;
    }
    cell.template get<VELOCITY<T, LatSet::d>>() = uc;

    // 3. MRT relaxation time vector
    //    D3Q19: s_q for q-moments, omega for others (Fortran-aligned)
    //    D2Q9: aligned with Fortran BubbleRising.f90 L773/L817
    //          Sf = (0, omega, omega, 0, s_q, 0, s_q, omega, omega)
    //    其他格子：全设为 omega，使 MRTForce 数值等价于 BGKForce
    //    稳定后可按需调大指定矩的 s_k 来利用 MRT 优势
    T rtvec[LatSet::q] {};
    if constexpr (std::is_same_v<LatSet, D3Q19<T>>) {
      // D3Q19 MRT relaxation: s_q for q-moments, omega for others
      T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
      rtvec[0] = T{0};       // rho
      rtvec[1] = omega;      // e
      rtvec[2] = omega;      // epsilon
      rtvec[3] = T{0};       // jx
      rtvec[4] = s_q;        // qx
      rtvec[5] = T{0};       // jy
      rtvec[6] = s_q;        // qy
      rtvec[7] = T{0};       // jz
      rtvec[8] = s_q;        // qz
      for (unsigned int i = 9; i < LatSet::q; ++i) rtvec[i] = omega;  // stress/higher-order
    } else if constexpr (std::is_same_v<LatSet, D2Q9<T>>) {
      // D2Q9 MRT relaxation: aligned with Fortran BubbleRising.f90 L773/L817
      // Sf = (0, omega, omega, 0, s_q, 0, s_q, omega, omega)
      T s_q = T{8} * (T{2} - omega) / (T{8} - omega);
      rtvec[0] = T{0};       // rho
      rtvec[1] = omega;      // e
      rtvec[2] = omega;      // eps
      rtvec[3] = T{0};       // jx
      rtvec[4] = s_q;        // qx
      rtvec[5] = T{0};       // jy
      rtvec[6] = s_q;        // qy
      rtvec[7] = omega;      // pxx
      rtvec[8] = omega;      // pxy
    } else {
      for (unsigned int i = 0; i < LatSet::q; ++i) rtvec[i] = omega;
    }

    // 4. Moments from population
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];

    // 5. Equilibrium in moment space
    std::array<T, LatSet::q> feq{};
    if constexpr (isIncompressible) {
      // Incompressible (He-Luo): feq = w*(p + 3*cu + 4.5*cu² - 1.5*u²), p = Σf
      equilibrium::IncompressibleSecondOrder<CELL>::apply(feq, zeroth, uc);
    } else {
      // Compressible: feq = w*rho*(1 + 3*cu + 4.5*cu² - 1.5*u²)
      equilibrium::SecondOrder<CELL>::apply(feq, rho, uc);
    }
    T momentaEq[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i)
      for (unsigned int j = 0; j < LatSet::q; ++j)
        momentaEq[i] += mrt::M<LatSet>(i, j) * feq[j];

    // 6. Force source in population space (Guo force scheme)
    T uF = T{};
    for (unsigned int d = 0; d < LatSet::d; ++d) uF += uc[d] * F[d];
    T fi_force[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T cF = T{}, cu = T{};
      for (unsigned int d = 0; d < LatSet::d; ++d) {
        cF += latset::c<LatSet>(i)[d] * F[d];
        cu += latset::c<LatSet>(i)[d] * uc[d];
      }
      // Incompressible Guo force needs /rho factor (Fortran L824: DF = w*[...] * dt/rho, dt=1)
      // Compressible Guo force has no /rho factor
      if constexpr (isIncompressible) {
        fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                       + LatSet::InvCs4 * cu * cF) / rho;
      } else {
        fi_force[i] = latset::w<LatSet>(i) * (LatSet::InvCs2 * (cF - uF)
                       + LatSet::InvCs4 * cu * cF);
      }
    }

    // 7. Force source in moment space
    T force_m[LatSet::q] {};
    for (unsigned int j = 0; j < LatSet::q; ++j)
      for (unsigned int k = 0; k < LatSet::q; ++k)
        force_m[j] += mrt::M<LatSet>(j, k) * fi_force[k];

    // 8. Combined MRT collision + force source
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      T source{};
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        T invM_ij = mrt::InvM<LatSet>(i, j);
        coll += invM_ij * rtvec[j] * (momenta[j] - momentaEq[j]);
        source += invM_ij * (T{1} - rtvec[j] / T{2}) * force_m[j];
      }
      cell[i] = cell[i] - coll + source;
    }
  }
};

// MRT collision for magnetic field potential equation (D2Q5).
//
// Solves (1/ε) ∂ψ/∂t = ∇·(μ∇ψ) via pseudo-time relaxation.
// M_5 moments: m = (ψ, J_x, J_y, e, p_xx).
// Relaxation: S = diag(1, ω_mag, ω_mag, 1, 1).
//   ψ:   conserved (s=1, equivalent to s=0 since m[0]=m_eq[0]=ψ identically)
//   J_x, J_y: relaxed by ω_mag = 1/τ_mag (diffusion rate)
//   e, p_xx:  fully relaxed to equilibrium (=0) with s=1
//
// Equilibrium: h_α^eq = w_eq * ψ with w_eq = 1/q.
// All non-conserved equilibrium moments are identically zero (M_5 rows sum to 0
// when acting on constant feq), so no explicit equilibrium computation is needed.
template <typename CELL, bool WriteToField = false>
struct MRTMag {
  using T = typename CELL::FloatType;
  using LatSet = typename CELL::LatticeSet;
  using equilibriumscheme = equilibrium::MagEquilibrium<CELL>;
  using GenericRho = typename CELL::GenericRho;

  __any__ static void apply(CELL& cell) {
    // Per-cell omega for variable-coefficient diffusion ∇·(μ∇ψ):
    // ω(x) = 1/(0.5 + μ(x)/(ε·cs²))
    // Falls back to lattice-level omega if OMEGA field is not in TypePack.
    T omega_mag;
    if constexpr (cell.template hasField<OMEGA<T>>()) {
      omega_mag = cell.template get<OMEGA<T>>();
    } else {
      omega_mag = cell.getOmega();  // ω = 1/τ, τ = 0.5 + ε·μ_avg/cs²
    }

    // Relaxation vector: S = (1, ω_mag, ω_mag, 1, 1)
    T rtvec[LatSet::q];
    rtvec[0] = T{1};
    if constexpr (LatSet::q > 1) {
      rtvec[1] = omega_mag;
      rtvec[2] = omega_mag;
      for (unsigned int i = 3; i < LatSet::q; ++i) rtvec[i] = T{1};
    }

    // Build InvM_S = InvM * diag(rtvec)
    T InvM_S[LatSet::q][LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        InvM_S[i][j] = mrt::InvM<LatSet>(i, j) * rtvec[j];
      }
    }

    // Forward transform: m = M * f
    T momenta[LatSet::q] {};
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      for (unsigned int j = 0; j < LatSet::q; ++j) {
        momenta[i] += mrt::M<LatSet>(i, j) * cell[j];
      }
    }

    // All non-conserved equilibrium moments are identically zero:
    //   m_eq[j] = 0 for j > 0  (proven for constant feq with M_5)
    //   m[0] = m_eq[0] = ψ   (conserved moment, s=1 → no relaxation)
    //
    // Collision: cell[i] -= Σ_{j>0} InvM_S[i][j] * m[j]
    // The conserved moment j=0 is skipped because m[0] = m_eq[0] = ψ
    // identically, so (m[0] - m_eq[0]) = 0 contributes nothing.
    for (unsigned int i = 0; i < LatSet::q; ++i) {
      T coll{};
      for (unsigned int j = 1; j < LatSet::q; ++j) {
        coll += InvM_S[i][j] * momenta[j];
      }
      cell[i] -= coll;
    }
  }
};

}