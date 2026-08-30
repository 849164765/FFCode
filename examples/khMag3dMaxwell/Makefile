TARGET = khMag3dMaxwell.exe
ROOT := ../..

FLAGS += -g -O2 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

CXXC := mpic++

# 仅编译主程序
override SRCS = khMag3dMaxwell.cpp

include $(ROOT)/make.mk
