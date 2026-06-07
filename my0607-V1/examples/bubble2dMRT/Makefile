TARGET = bubble2dMRT.exe
ROOT := ../..

FLAGS += -g -O0 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
# FLAGS += -O3 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs -march=native
# FLAGS += -O3 -fopenmp
# FLAGS += -g -DMPI_DEBUG
# FLAGS += -g
# FLAGS += -pg
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

# use mpic++ to compile
CXXC := mpic++

include $(ROOT)/make.mk
