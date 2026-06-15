# tests/Makefile
# Build and run unit tests

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -g -O0
ROOT     := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))..)
INCLUDES := -I$(ROOT)/src/
DEFINES  := -DThread_Num=4

TARGETS := test_collision test_pf

.PHONY: all compile run clean

all: compile run

compile: $(TARGETS)

test_collision: test_collision.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) $< -o $@

test_pf: test_pf.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(DEFINES) $< -o $@

run: $(TARGETS)
	@echo "==> test_collision"
	./test_collision
	@echo ""
	@echo "==> test_pf"
	./test_pf

clean:
	rm -f $(TARGETS)
