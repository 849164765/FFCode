// test_mrt_matrix.cpp — verify D2Q5 MRT matrices M·InvM = I
#include <iostream>
#include <cmath>

int main() {
  constexpr int Q = 5;

  // M<2,5> matrix (from collisionMRT.h)
  constexpr double M[Q][Q] = {
    { 1,  1,  1,  1,  1},
    { 0,  1, -1,  0,  0},
    { 0,  0,  0,  1, -1},
    {-4,  1,  1,  1,  1},
    { 0,  1,  1, -1, -1}
  };

  // InvM<2,5> matrix (from collisionMRT.h)
  constexpr double InvM[Q][Q] = {
    { 1.0/5,    0.0,     0.0, -1.0/5,     0.0},
    { 1.0/5,  1.0/2,     0.0,  1.0/20,  1.0/4},
    { 1.0/5, -1.0/2,     0.0,  1.0/20,  1.0/4},
    { 1.0/5,    0.0,   1.0/2,  1.0/20, -1.0/4},
    { 1.0/5,    0.0,  -1.0/2,  1.0/20, -1.0/4}
  };

  // Compute M · InvM
  double result[Q][Q] = {{0}};
  for (int i = 0; i < Q; ++i)
    for (int j = 0; j < Q; ++j)
      for (int k = 0; k < Q; ++k)
        result[i][j] += M[i][k] * InvM[k][j];

  // Check identity
  double max_err = 0;
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j) {
      double expected = (i == j) ? 1.0 : 0.0;
      double err = std::abs(result[i][j] - expected);
      if (err > max_err) max_err = err;
    }
  }

  std::cout << "M · InvM =\n";
  for (int i = 0; i < Q; ++i) {
    for (int j = 0; j < Q; ++j)
      printf("%8.4f", result[i][j]);
    std::cout << "\n";
  }

  std::cout << "\nMax error: " << max_err << "\n";
  if (max_err < 1e-14) {
    std::cout << "PASS: M · InvM = I\n";
    return 0;
  } else {
    std::cout << "FAIL: M · InvM != I\n";
    return 1;
  }
}
