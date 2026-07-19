// compile_check_mfield.cpp — verify mfield2d.h compiles
#include "freelb.h"
#include "freelb.hh"

using T = double;
using MFLatSet = D2Q5<T>;

int main() {
  // Verify types exist
  mfield::PSI<T>* p = nullptr;          (void)p;
  mfield::HX<T>* hx = nullptr;          (void)hx;
  mfield::HY<T>* hy = nullptr;          (void)hy;
  mfield::HMAG<T>* hm = nullptr;        (void)hm;
  mfield::OMEGA_PSI<T>* o = nullptr;    (void)o;
  mfield::MU_PERCELL<T>* mu = nullptr;  (void)mu;
  mfield::CHI_PERCELL<T>* chi = nullptr;(void)chi;

  // Data types
  mfield::MU_L<T>* mu_l = nullptr;      (void)mu_l;
  mfield::MU_H<T>* mu_h = nullptr;      (void)mu_h;
  mfield::CHI_L<T>* chi_l = nullptr;    (void)chi_l;
  mfield::CHI_H<T>* chi_h = nullptr;    (void)chi_h;
  mfield::H_0<T>* h0 = nullptr;         (void)h0;

  std::cout << "PASS: mfield2d.h all types compile" << std::endl;
  return 0;
}
