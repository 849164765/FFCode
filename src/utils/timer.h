#pragma once

// std::chrono::high_resolution_clock, std::chrono::milliseconds
#include <chrono>
// std::int64_t
#include <cstdint>
// std::cout, std::endl
#include <iostream>

// Thread_Num
#include "head.h"
// MPI_RANK(x), mpi().getSize()
#include "parallel/mpi_manager.h"

// timer and counters

struct Counter {
  std::int64_t count;
  Counter() { count = 0; }
  ~Counter() = default;
  // reset count
  void reset() { count = 0; }
  // define ++ operator, use: Counter c; ++c; or c++;
  void operator++() { count++; }
  // get count
  // define () operator, use: Counter c; c(); or c.operator()();
  std::int64_t operator()() const { return count; }
};

struct Timer : public Counter {
  std::chrono::high_resolution_clock::time_point START;
  std::chrono::high_resolution_clock::time_point END;

  Timer() { START_TIMER(); }
  ~Timer() = default;

  void START_TIMER() { START = std::chrono::high_resolution_clock::now(); }
  void END_TIMER() { END = std::chrono::high_resolution_clock::now(); }
  void reset() {
    START_TIMER();
    Counter::reset();
  }

  std::chrono::milliseconds GetDuration() {
    END_TIMER();
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START);
  }
  // get time elapsed in milliseconds(ms)
  std::int64_t GetDurationCount() {
    END_TIMER();
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START).count();
  }
  // get time elapsed in milliseconds(ms) but no END_TIMER()
  std::int64_t GetDurationCount_Only() {
    return std::chrono::duration_cast<std::chrono::milliseconds>(END - START).count();
  }
  // get time elapsed in seconds(s)
  double GetTimeElapsed() { return GetDurationCount() / double(1000); }
  // get time elapsed in seconds(s) but no END_TIMER()
  double GetTimeElapsed_Only() { return GetDurationCount_Only() / double(1000); }

  void Print_MainLoopPerformance(std::size_t n) {
    MPI_RANK(0)
    END_TIMER();
    std::int64_t count = this->operator()();  // this->operator()() or (*this)()
    double MLUPs = static_cast<double>(count) * static_cast<double>(n) /
                   GetDurationCount() / double(1000);
    std::cout << "[Main_Loop Performance]:"
              << "\n"
              << "Time Elapsed:  " << GetTimeElapsed_Only() << " s"
              << "\n"
              << "Total Step:    " << count << "\n"
              << "Average_MLUPs: " << MLUPs << std::endl;
#ifdef _OPENMP
    std::cout << "MLUPs/Thread:  " << MLUPs / Thread_Num << std::endl;
#endif
#ifdef MPI_ENABLED
    std::cout << "MLUPs/Process: " << MLUPs / mpi().getSize() << std::endl;
#endif
  }

  void Print_InnerLoopPerformance(std::size_t n, std::size_t steps) {
    MPI_RANK(0)
    std::int64_t count = this->operator()();  // this->operator()() or (*this)()
    double MLUPs = static_cast<double>(steps) * static_cast<double>(n) /
                   GetDurationCount() / double(1000);
    std::cout << "[Step: " << count << "]  "
              << "MLUPs: " << MLUPs << "  ";
    // this->reset();
    START_TIMER();
  }
};