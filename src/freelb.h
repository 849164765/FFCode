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
// (each example includes the specific mfield header it needs)