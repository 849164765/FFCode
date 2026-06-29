// cell_dynamics.h
#pragma once

#include "lbm/collision.h"
#include "lbm/moment.h"

// a template for cell dynamics
// you could use it like:
// using MyDynamics = CellDynamics<task1, task2, task3, ...>;
template <typename... Tasks>
struct CellDynamics {
  template <typename CELL, typename... Args>
  static void Apply(CELL& cell, Args... args) {
    (applyHelper<Tasks>(cell, args...), ...);
  }

  template <typename Task, typename CELL, typename FirstArg, typename... Args>
  static auto applyHelper(CELL& cell, FirstArg firstArg, Args... args)
    -> std::enable_if_t<
      std::is_invocable_v<decltype(&Task::apply), CELL&, FirstArg, Args...>, void> {
    Task::apply(cell, firstArg, args...);
  }
  template <typename Task, typename CELL, typename FirstArg, typename... Args>
  static auto applyHelper(CELL& cell, FirstArg firstArg, Args... args)
    -> std::enable_if_t<
      !std::is_invocable_v<decltype(&Task::apply), CELL&, FirstArg, Args...>, void> {
    applyHelper<Task>(cell, args...);
  }
  template <typename Task, typename CELL>
  static auto applyHelper(CELL& cell)
    -> std::enable_if_t<std::is_invocable_v<decltype(&Task::apply), CELL&>, void> {
    Task::apply(cell);
  }
};

namespace dynamics {

}  // namespace dynamics