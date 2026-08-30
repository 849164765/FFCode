TARGET = kh2dStd.exe
ROOT := ../..

FLAGS += -g -O2 -Wall  -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

CXXC := mpic++

include $(ROOT)/make.mk
