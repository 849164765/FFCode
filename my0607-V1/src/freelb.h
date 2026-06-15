//freelb.h
//include haeder files for LBM

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

#include "lbm/lattice_set.h"
// lbm dynamics
#include "lbm/lbm.h"

// fflbm — ferrofluid LBM operators
#include "fflbm/phasefield.h"
#include "fflbm/phasefield_source.h"
#include "fflbm/magnetic_field.h"
#include "fflbm/magnetic_boundary.h"
#include "fflbm/coupling.h"