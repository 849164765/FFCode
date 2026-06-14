# tests/Makefile — build individual test executables
ROOT := ..

FLAGS += -g -O0 -Wall -Wno-unknown-pragmas -Wno-unused-local-typedefs
FLAGS += -DFLOAT_TYPE=double
FLAGS += -D_UNROLLFOR
FLAGS += -std=c++17 -fno-diagnostics-show-template-tree

CXXC ?= g++
LINKFLAGS := -L$(ROOT)/lib -lrt

TESTS := test_phaseA test_phaseB test_phaseC

.PHONY: all clean run-all

all: $(addsuffix .exe, $(TESTS))

%.exe: %.o
	$(CXXC) $(FLAGS) -o $@ $^ $(LINKFLAGS)

%.o: %.cpp
	$(CXXC) $(FLAGS) -I$(ROOT)/src/ -c $< -o $@

clean:
	rm -rf *.o *.d *.exe output

run-all: $(addsuffix .exe, $(TESTS))
	@for t in $(TESTS); do \
		echo "=== Running $$t ==="; \
		./$$t.exe || exit 1; \
	done
	@echo "All tests passed."
