// mrt_matrix.h
// D2Q9 and D2Q5 MRT transformation matrices and functions
// Paper: Guo et al., Phys. Fluids 37, 022148 (2025), Eq.20 (M9), Eq.41 (M5)
//
// D2Q9 velocity order (per lattice_set.h D2Q9<T>::c[]):
//   [0]=(0,0), [1]=(1,0), [2]=(-1,0), [3]=(0,1), [4]=(0,-1),
//   [5]=(1,1), [6]=(-1,-1), [7]=(1,-1), [8]=(-1,1)
//
// D2Q5 velocity order (per lattice_set.h D2Q5<T>::c[]):
//   [0]=(0,0), [1]=(1,0), [2]=(-1,0), [3]=(0,1), [4]=(0,-1)

#pragma once

namespace mrt {

// ================================================================
// M9: 9x9 forward transformation matrix  (m = M9 * f)
// Row i: moment i, Column j: population f_j (c[] order above)
// Moments: [0]=rho, [1]=e, [2]=epsilon, [3]=jx, [4]=qx,
//           [5]=jy, [6]=qy, [7]=pxx, [8]=pxy
// ================================================================
template <typename T>
constexpr T M9[9][9] = {
  // rho
  {1,  1,  1,  1,  1,  1,  1,  1,  1},
  // e
  {4, -1, -1, -1, -1,  2,  2,  2,  2},
  // epsilon
  {4, -2, -2, -2, -2,  1,  1,  1,  1},
  // jx  (columns 2/3 swapped from paper Eq.20 to match our c[] order)
  {0,  1, -1,  0,  0,  1, -1, -1,  1},
  // qx  (columns 2/3 swapped)
  {0, -2,  2,  0,  0,  1, -1, -1,  1},
  // jy  (columns 2/3 swapped)
  {0,  0,  0,  1, -1,  1,  1, -1, -1},
  // qy  (columns 2/3 swapped)
  {0,  0,  0, -2,  2,  1,  1, -1, -1},
  // pxx (columns 2/3 swapped)
  {0,  1,  1, -1, -1,  0,  0,  0,  0},
  // pxy
  {0,  0,  0,  0,  0,  1, -1,  1, -1}
};

// ================================================================
// invM9: 9x9 inverse transformation matrix (f = invM9 * m)
// Computed via Gaussian elimination, verified M9 * invM9 = I
// ================================================================
template <typename T>
constexpr T invM9[9][9] = {
  { T(1),      -T(1),       T(1),       T(0),   T(0),
    T(0),       T(0),       T(0),       T(0)},
  { T(1)/T(3), -T(1)/T(4),  T(1)/T(6),  T(1)/T(6), -T(1)/T(6),
    T(0),       T(0),       T(1)/T(4),  T(0)},
  { T(1)/T(3), -T(1)/T(4),  T(1)/T(6), -T(1)/T(6),  T(1)/T(6),
    T(0),       T(0),       T(1)/T(4),  T(0)},
  { T(1)/T(3), -T(1)/T(4),  T(1)/T(6),  T(0),   T(0),
    T(1)/T(6), -T(1)/T(6), -T(1)/T(4),  T(0)},
  { T(1)/T(3), -T(1)/T(4),  T(1)/T(6),  T(0),   T(0),
   -T(1)/T(6),  T(1)/T(6), -T(1)/T(4),  T(0)},
  {-T(1)/T(3),  T(1)/T(2), -T(5)/T(12), T(1)/T(6),  T(1)/T(12),
    T(1)/T(6),  T(1)/T(12), T(0),       T(1)/T(4)},
  {-T(1)/T(3),  T(1)/T(2), -T(5)/T(12),-T(1)/T(6), -T(1)/T(12),
    T(1)/T(6),  T(1)/T(12), T(0),      -T(1)/T(4)},
  {-T(1)/T(3),  T(1)/T(2), -T(5)/T(12),-T(1)/T(6), -T(1)/T(12),
   -T(1)/T(6), -T(1)/T(12), T(0),       T(1)/T(4)},
  {-T(1)/T(3),  T(1)/T(2), -T(5)/T(12), T(1)/T(6),  T(1)/T(12),
   -T(1)/T(6), -T(1)/T(12), T(0),      -T(1)/T(4)}
};

// ================================================================
// M5: 5x5 forward transformation matrix for D2Q5  (m = M5 * h)
// Moments: [0]=psi, [1]=jx, [2]=jy, [3]=e, [4]=pxy
// Columns 0..4 correspond to D2Q5 c[]:
//   [0]=(0,0), [1]=(1,0), [2]=(-1,0), [3]=(0,1), [4]=(0,-1)
// ================================================================
template <typename T>
constexpr T M5[5][5] = {
  // psi
  {1,  1,  1,  1,  1},
  // jx
  {0,  1, -1,  0,  0},
  // jy
  {0,  0,  0,  1, -1},
  // e
  {4, -1, -1, -1, -1},
  // pxy
  {0,  1,  1, -1, -1}
};

// ================================================================
// invM5: 5x5 inverse transformation matrix (h = invM5 * m)
// Computed analytically, verified M5 * invM5 = I
// ================================================================
template <typename T>
constexpr T invM5[5][5] = {
  // h0
  { T(1)/T(5),  T(0),       T(0),        T(1)/T(5),   T(0)},
  // h1
  { T(1)/T(5),  T(1)/T(2),  T(0),       -T(1)/T(20),  T(1)/T(4)},
  // h2
  { T(1)/T(5), -T(1)/T(2),  T(0),       -T(1)/T(20),  T(1)/T(4)},
  // h3
  { T(1)/T(5),  T(0),       T(1)/T(2),  -T(1)/T(20), -T(1)/T(4)},
  // h4
  { T(1)/T(5),  T(0),      -T(1)/T(2),  -T(1)/T(20), -T(1)/T(4)}
};

// ================================================================
// M9Transform:  forward  m = M9 * f  (hand-unrolled)
// ================================================================
template <typename T>
inline void M9Transform(T m[9], const T f[9]) {
  m[0] = f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8];
  m[1] = T(4)*f[0] - f[1] - f[2] - f[3] - f[4]
       + T(2)*(f[5] + f[6] + f[7] + f[8]);
  m[2] = T(4)*f[0] - T(2)*(f[1] + f[2] + f[3] + f[4])
       + f[5] + f[6] + f[7] + f[8];
  m[3] = f[1] - f[2] + f[5] - f[6] - f[7] + f[8];
  m[4] = T(-2)*f[1] + T(2)*f[2] + f[5] - f[6] - f[7] + f[8];
  m[5] = f[3] - f[4] + f[5] + f[6] - f[7] - f[8];
  m[6] = T(-2)*f[3] + T(2)*f[4] + f[5] + f[6] - f[7] - f[8];
  m[7] = f[1] + f[2] - f[3] - f[4];
  m[8] = f[5] - f[6] + f[7] - f[8];
}

// ================================================================
// M9InverseTransform:  inverse  f = invM9 * m  (hand-unrolled)
// ================================================================
template <typename T>
inline void M9InverseTransform(T f[9], const T m[9]) {
  f[0] = m[0] - m[1] + m[2];
  f[1] = ( T(4)*m[0] - T(3)*m[1] + T(2)*m[2]
       +   T(2)*m[3] - T(2)*m[4] + T(3)*m[7]) / T(12);
  f[2] = ( T(4)*m[0] - T(3)*m[1] + T(2)*m[2]
       -   T(2)*m[3] + T(2)*m[4] + T(3)*m[7]) / T(12);
  f[3] = ( T(4)*m[0] - T(3)*m[1] + T(2)*m[2]
       +   T(2)*m[5] - T(2)*m[6] - T(3)*m[7]) / T(12);
  f[4] = ( T(4)*m[0] - T(3)*m[1] + T(2)*m[2]
       -   T(2)*m[5] + T(2)*m[6] - T(3)*m[7]) / T(12);
  f[5] = (-T(4)*m[0] + T(6)*m[1] - T(5)*m[2]
       +   T(2)*m[3] +     m[4] + T(2)*m[5] +     m[6] + T(3)*m[8]) / T(12);
  f[6] = (-T(4)*m[0] + T(6)*m[1] - T(5)*m[2]
       -   T(2)*m[3] -     m[4] + T(2)*m[5] +     m[6] - T(3)*m[8]) / T(12);
  f[7] = (-T(4)*m[0] + T(6)*m[1] - T(5)*m[2]
       -   T(2)*m[3] -     m[4] - T(2)*m[5] -     m[6] + T(3)*m[8]) / T(12);
  f[8] = (-T(4)*m[0] + T(6)*m[1] - T(5)*m[2]
       +   T(2)*m[3] +     m[4] - T(2)*m[5] -     m[6] - T(3)*m[8]) / T(12);
}

// ================================================================
// M5Transform:  forward  m = M5 * h  (hand-unrolled)
// ================================================================
template <typename T>
inline void M5Transform(T m[5], const T h[5]) {
  m[0] = h[0] + h[1] + h[2] + h[3] + h[4];
  m[1] = h[1] - h[2];
  m[2] = h[3] - h[4];
  m[3] = T(4)*h[0] - h[1] - h[2] - h[3] - h[4];
  m[4] = h[1] + h[2] - h[3] - h[4];
}

// ================================================================
// M5InverseTransform:  inverse  h = invM5 * m  (hand-unrolled)
// ================================================================
template <typename T>
inline void M5InverseTransform(T h[5], const T m[5]) {
  h[0] = ( m[0] + m[3]) / T(5);
  h[1] = ( T(4)*m[0] + T(10)*m[1] - m[3] + T(5)*m[4]) / T(20);
  h[2] = ( T(4)*m[0] - T(10)*m[1] - m[3] + T(5)*m[4]) / T(20);
  h[3] = ( T(4)*m[0] + T(10)*m[2] - m[3] - T(5)*m[4]) / T(20);
  h[4] = ( T(4)*m[0] - T(10)*m[2] - m[3] - T(5)*m[4]) / T(20);
}

// ================================================================
// PrecomputeInvM9S: invM9_S = invM9 * diag(S) * M9
// For D2Q9 MRT collision: f -= invM9_S * (f - f_eq)
// ================================================================
template <typename T>
inline void PrecomputeInvM9S(T invM9_S[9][9], const T S[9]) {
  for (int i = 0; i < 9; ++i) {
    for (int j = 0; j < 9; ++j) {
      T sum = T(0);
      for (int k = 0; k < 9; ++k) {
        sum += invM9<T>[i][k] * S[k] * M9<T>[k][j];
      }
      invM9_S[i][j] = sum;
    }
  }
}

// ================================================================
// PrecomputeInvM5S: invM5_S = invM5 * diag(S) * M5
// For D2Q5 MRT collision: h -= invM5_S * (h - h_eq)
// ================================================================
template <typename T>
inline void PrecomputeInvM5S(T invM5_S[5][5], const T S[5]) {
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      T sum = T(0);
      for (int k = 0; k < 5; ++k) {
        sum += invM5<T>[i][k] * S[k] * M5<T>[k][j];
      }
      invM5_S[i][j] = sum;
    }
  }
}

// ================================================================
// MRTApply: apply MRT collision in population space
// f_out = f - invM_S * (f - f_eq)
// ================================================================
template <typename T, int Q>
inline void MRTApply(T f[Q], const T f_eq[Q], const T invM_S[Q][Q]) {
  T delta[Q];
  for (int i = 0; i < Q; ++i) {
    delta[i] = f[i] - f_eq[i];
  }
  for (int i = 0; i < Q; ++i) {
    T correction = T(0);
    for (int j = 0; j < Q; ++j) {
      correction += invM_S[i][j] * delta[j];
    }
    f[i] -= correction;
  }
}

}  // namespace mrt
