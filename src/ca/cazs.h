// basic cazs(cellular automata zhu-stefanescu) header file

#pragma once

#include "data_struct/field_struct.h"

namespace CA {

enum CAType : std::uint8_t { Boundary = 1, Interface = 2, Fluid = 4, Solid = 8 };

enum CAFlag : std::uint8_t { None = 1, toInterface = 2, toSolid = 4 };

}  // namespace CA