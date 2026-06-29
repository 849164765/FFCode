// plot_writer.h

#pragma once

#include "io/basic_writer.h"

template <typename T>
class FLBplot {
 private:
  // x-axis: step
  // y-axis: data
  std::string name;
  std::string xvar;
  std::string yvar;
  std::string work_dir;

 public:
  FLBplot(std::string work_dir_, std::string name = "/FLBplot",
          std::string xvar = "step", std::string yvar = "data")
      : work_dir(work_dir_), name(name), xvar(xvar), yvar(yvar) {
    plotHeader();
  }
  void plotHeader() {
    std::ofstream write;
    write.open(work_dir + name + ".dat");
    write << xvar << "\t" << yvar << "\n";
  }
  // each time call this function, a step-data pair is written
  template <typename X = int, typename U = T>
  void write(X xvalue, U yvalue) {
    std::ofstream write;
    write.open(work_dir + name + ".dat", std::ios::app);
    write << xvalue << "\t" << yvalue << "\n";
  }
  template <typename U = T>
  void write(Counter &counter, U yvalue) {
    write<int, U>(counter(), yvalue);
  }
};