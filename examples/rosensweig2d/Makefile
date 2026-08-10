TARGET = rosensweig2d.exe
ROOT := ../..

FLAGS += -O2 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput
CXXC := mpic++

include $(ROOT)/make.mk
