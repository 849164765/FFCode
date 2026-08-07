TARGET = bubbleMag3d.exe
ROOT := ../..

FLAGS += -g -O2 -Wall -Wdisabled-optimization -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -DMPI_ENABLED
FLAGS += -D_UNROLLFOR

OUTPUT += vtkoutput

CXXC := mpic++

# 仅编译主程序 (field_validate/droplet_validate 为独立验证工具, 各自链接)
# override: 防止 make.mk 的 SRCS=$(wildcard *.cpp) 把验证程序一起编译进来
override SRCS = bubbleMag3d.cpp

include $(ROOT)/make.mk
