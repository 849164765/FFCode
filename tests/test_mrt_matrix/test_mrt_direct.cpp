// test_mrt_direct.cpp — verify D2Q5 MRT matrix values directly
#include <iostream>
#include <cmath>

// Reproduce min setup from head.h
#define __constant__
#define __any__
#define __constexpr__ constexpr

// Minimal Fraction
template <typename T = int>
struct Fraction {
  T numerator;
  T denominator;
  constexpr Fraction(T num, T den) : numerator(num), denominator(den) {}
  constexpr Fraction(T num) : numerator(num), denominator(1) {}
  constexpr Fraction() : numerator(0), denominator(1) {}
  template <typename U>
  constexpr U operator()() const { return U(numerator) / U(denominator); }
};

// D2Q5 MRT matrices
constexpr Fraction<> M25[5][5] = {
  { 1,  1,  1,  1,  1},
  { 0,  1, -1,  0,  0},
  { 0,  0,  0,  1, -1},
  {-4,  1,  1,  1,  1},
  { 0,  1,  1, -1, -1}
};

constexpr Fraction<> InvM25[5][5] = {
  {{1, 5},     {0},      {0}, {-1, 5},     {0}},
  {{1, 5},  {1, 2},      {0}, {1, 20},  {1, 4}},
  {{1, 5}, {-1, 2},      {0}, {1, 20},  {1, 4}},
  {{1, 5},     {0},   {1, 2}, {1, 20}, {-1, 4}},
  {{1, 5},     {0},  {-1, 2}, {1, 20}, {-1, 4}}
};

int main() {
  // Check M values
  std::cout << "M[0][0] = " << M25[0][0].template operator()<double>() << std::endl;
  std::cout << "M[3][0] = " << M25[3][0].template operator()<double>() << std::endl;

  // Check InvM values
  std::cout << "InvM[0][0] = " << InvM25[0][0].template operator()<double>() << std::endl;
  std::cout << "InvM[1][1] = " << InvM25[1][1].template operator()<double>() << std::endl;

  // Check M·InvM = I
  constexpr int Q = 5;
  double result[Q][Q] = {{0}};
  for (int i = 0; i < Q; ++i)
    for (int j = 0; j < Q; ++j)
      for (int k = 0; k < Q; ++k)
        result[i][j] += M25[i][k].template operator()<double>() * InvM25[k][j].template operator()<double>();

  double max_err = 0;
  for (int i = 0; i < Q; ++i)
    for (int j = 0; j < Q; ++j) {
      double expected = (i==j) ? 1.0 : 0.0;
      double err = std::abs(result[i][j] - expected);
      if (err > max_err) max_err = err;
    }

  std::cout << "\nM·InvM max error: " << max_err << std::endl;
  if (max_err < 1e-14) std::cout << "PASS" << std::endl;
  else std::cout << "FAIL" << std::endl;

  return (max_err < 1e-14) ? 0 : 1;
}
