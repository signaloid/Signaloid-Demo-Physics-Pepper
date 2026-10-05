#
#	Top-level Makefile for native Monte Carlo builds.
#
#	Build the native Monte Carlo executable with:
#
#		make local-build
#
#	This builds the demo against the host toolchain and the UxHw
#	compatibility shim in `submodules/compat`, rather than against the
#	Signaloid cloud compiler.
#
#	Copyright (c) 2026, Signaloid.
#
#	Targets:
#	  local-build - Build for native execution with the UxHw compat layer
#	  local-run   - Build and run one distributional pass
#	  clean       - Remove build artifacts
#

CXX      = g++
CC       = gcc

#
#	Define `_GNU_SOURCE`, matching `config.mk`, so that `strtok_r` in
#	`common.c` and `random` in `uxhw.c` are declared. Without it both are
#	implicitly declared as returning `int`, which truncates a pointer.
#
CXXFLAGS = -std=c++11 -O2 -Isrc/ -D_GNU_SOURCE
CFLAGS   = -std=c11 -O2 -Isrc/ -D_GNU_SOURCE

LDFLAGS  =
LIBS     = -lgsl -lgslcblas -lm

# Platform-specific paths (macOS with MacPorts or Homebrew)
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    ifneq ($(wildcard /opt/local/include/gsl),)
        CXXFLAGS += -I/opt/local/include
        CFLAGS   += -I/opt/local/include
        LDFLAGS  += -L/opt/local/lib
    else ifneq ($(wildcard /opt/homebrew/include/gsl),)
        CXXFLAGS += -I/opt/homebrew/include
        CFLAGS   += -I/opt/homebrew/include
        LDFLAGS  += -L/opt/homebrew/lib
    else ifneq ($(wildcard /usr/local/include/gsl),)
        CXXFLAGS += -I/usr/local/include
        CFLAGS   += -I/usr/local/include
        LDFLAGS  += -L/usr/local/lib
    endif
endif

BINARY = demo-native-mc

# ──────────────────────────────────────────────
# Project sources
#
# The PDF backend is `pepper-lhapdfLite.c`, a self-contained C reader for
# the CT18NLO grid in `inputs/`. There is no dependency on LHAPDF6.
# ──────────────────────────────────────────────

PROJECT_CPP = \
	src/pepper.cpp \
	src/pepper-utilities.cpp \
	src/pepper-vec4.cpp \
	src/pepper-complexArith.cpp \
	src/pepper-berendsGiele.cpp \
	src/pepper-phaseSpace.cpp \
	src/pepper-kernel.cpp \
	src/pepper-uxhw.cpp \
	src/pepper-monteCarlo.cpp

PROJECT_C = \
	src/common.c \
	src/uxhw.c \
	src/pepper-lhapdfLite.c

# ──────────────────────────────────────────────
# Object files
# ──────────────────────────────────────────────

PROJECT_CPP_OBJ = $(PROJECT_CPP:.cpp=.o)
PROJECT_C_OBJ   = $(PROJECT_C:.c=.o)

ALL_OBJ = $(PROJECT_CPP_OBJ) $(PROJECT_C_OBJ)

# ──────────────────────────────────────────────
# Compilation rules
# ──────────────────────────────────────────────

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

# ──────────────────────────────────────────────
# Link
# ──────────────────────────────────────────────

$(BINARY): $(ALL_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(ALL_OBJ) $(LIBS)

# ──────────────────────────────────────────────
# Targets
# ──────────────────────────────────────────────

.PHONY: local-build local-run clean

local-build: $(BINARY)

local-run: local-build
	./$(BINARY)

clean:
	rm -f $(ALL_OBJ) $(BINARY)
