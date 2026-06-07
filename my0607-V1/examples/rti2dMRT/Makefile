TARGET = rti2dMRT.exe
ROOT := ../..

FLAGS += -g -O0 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
# FLAGS += -O3 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs -march=native
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

# use mpic++ to compile
CXXC := mpic++

include $(ROOT)/make.mk
