TARGET = rosen3d.exe
ROOT := ../..

FLAGS += -g -O2 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

CXXC := mpic++

override SRCS = rosen3d.cpp

include $(ROOT)/make.mk
