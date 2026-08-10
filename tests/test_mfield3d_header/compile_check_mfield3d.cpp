// compile_check_mfield3d.cpp - verify the 3D magnetic field API compiles
#include "freelb.h"
#include "freelb.hh"

using T = double;

int main() {
  mfield::HZ<T>* hz = nullptr;
  mfield::MFFIELDS3D<T>* fields = nullptr;
  (void)hz;
  (void)fields;

  std::cout << "PASS: mfield3d.h types compile" << std::endl;
  return 0;
}
