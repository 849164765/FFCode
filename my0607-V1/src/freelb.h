//freelb.h
//include header files for LBM

// field type aliases (needed by tmp.h and lattice_set.h)
#include "utils/alias.h"

// lbm core (must come before boundary/geometry which depend on it)
#include "lbm/lattice_set.h"
#include "lbm/equilibrium.h"
#include "lbm/moment.h"
#include "lbm/collision.h"

// timer
#include "utils/timer.h"
// geometry
#include "geometry/geometry.h"
// boundary
#include "boundary/boundary.h"

// io
#include "io/vtkWriter.h"
#include "io/vtm_writer.h"
#include "io/ini_reader.h"
#include "io/vtu_writer.h"