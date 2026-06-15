#pragma once

#include "lbm/unit_converter.h"

template <typename T>
void BaseConverter<T>::check(int &check_status) {
  MPI_RANK(0)
  std::cout << "[BasicConverter]:\n"
            << "deltaX:           " << deltaX << "\n"
            << "deltaT:           " << deltaT << "\n"
            << "charL:            " << charL << "\n"
            << "charU:            " << charU << "\n"
            << "VisKine:          " << VisKine << "\n"
            << "Re:               " << Re << "\n"
            << "Conv_Acc:         " << Conv_Acc << "\n"
            << "Lattice_charU:    " << Lattice_charU << "\n"
            << "Lattice_charL:    " << Lattice_charL << "\n"
            << "Lattice_VisKine:  " << Lattice_VisKine << "\n"
            << "Lattice_Re:       " << Lattice_Re << "\n";
  if (Lattice_Re > 100) {
    std::cout << "[Warning]: Lattice_Re > 100, Press any key to continue..."
              << std::endl;
    getchar();
    check_status = 1;
  }

  if (Lattice_charU > 0.4) {
    std::cout << "Error: Lattice_charU > 0.4" << std::endl;
    exit(-1);
  } else {
    if (Lattice_RT < 0.55) {
      T Minstable_Lattice_RT =
          0.5 + 0.125 * Lattice_charU;  // 0.5 + 0.125 * 0.4 = 0.55, 0.5 + 0.125
                                        // * 0.03 = 0.50375
      if (Lattice_RT < Minstable_Lattice_RT) {
        std::cout << "Error: Lattice_charU too large, Minstable_Lattice_RT = "
                  << Minstable_Lattice_RT << std::endl;
        exit(-1);
      }
    }
  }
}

template <typename T>
void UnitConvManager<T>::Check_and_Print() {
  MPI_RANK(0)
  int check_status = 0;
  Printer::Print_Banner("Convert Log");

  BaseConv->check(check_status);

  if (check_status == 0) {
    std::cout << "Simulation parameters correctly set." << std::endl;
  } else {
    std::cout << "Simulation parameters set with warnings." << std::endl;
  }

  std::cout << "[Lattice Parameters]:\n"
            << "Lattice_RT =     " << BaseConv->getLattice_RT() << "\n"
            << "OMEGA =          " << BaseConv->getOMEGA() << "\n"
            << "-------------------------------------------------\n"
            << std::endl;
}