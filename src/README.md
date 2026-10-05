# Source Code

This directory holds the demo's source files. Build instructions live in the
[top-level README](../README.md). The entrypoint for both native and SCCE
builds is the top-level `Makefile` (`local-build`, `local-run`, `local-mc`,
`scce-build`, `scce-run`, `clean`).

## Demo sources

### `pepper.cpp`
Entry point. Parses command-line arguments, initializes LHAPDF, dispatches to
`calculateOutputMonteCarlo` (MC mode) or `calculateOutputUxHw` (UxHw® mode), and
formats the output.

### `pepper-kernel.cpp/h`
Shared LHAPDF initialization and cleanup, access to the PDF and `alpha_s`, and the
`calculateOutputUxHw` / `calculateOutputMonteCarlo` mode dispatch (contracting
to the first moment or to the sample mean respectively).

### `pepper-uxhw.cpp/h`
`generateXsPointUxHw`: the distributional cross-section kernel. PDF lookup,
`min`, `abs`, and step gates are wrapped in `PropagateFunction`. Returns the
differential cross section as a distribution, before any contraction to a
moment.

### `pepper-monteCarlo.cpp/h`
`generateXsPointMonteCarlo` / `generateXsMonteCarlo`: the same physics as plain
scalar arithmetic (`fmin`/`fmax`/`fabs`/step, direct PDF lookup), plus the Monte
Carlo sampling loop.

### `pepper-berendsGiele.cpp/h`
Berends-Giele recursion for 4-gluon color-ordered amplitudes. Currently unused,
because the kernels use the equivalent closed-form ME².

### `pepper-phaseSpace.cpp/h`
t-channel phase-space generator (arXiv:2302.10449). Currently unused, because
phase space is inline in the kernels.

### `pepper-complexArith.cpp/h`
Complex numbers, `CVec4` (complex four-vectors), `WeylSpinor`, and `VT` (the
sigma^mu contraction). The Weyl spinor uses a complex square root so that
negative-energy momenta are handled correctly.

### `pepper-vec4.cpp/h`
Real four-vectors with the Lorentz metric (+,-,-,-).

### `pepper-lhapdfLite.c/h`
Standalone C implementation of the LHAPDF grid reader, providing `xf(x, Q^2)`
and `alpha_s(Q^2)` by log-bicubic and cubic interpolation using only
`fopen`/`fgets`/`strtod` (no `iostream`, YAML, or exceptions). This is the PDF
backend for both native and SCCE builds, called directly by
`pepper-kernel.cpp`.

### `pepper-utilities.cpp/h`
Utility methods for parsing, setting, and reporting the usage of demo-specific
command-line arguments, along with the `OutputVariableIndex` enum. These methods
call similar methods from `common.c` for handling command-line arguments common
to all of our C/C++ demo applications.

## Shared sources

### `common.c/h`
These contain utility methods for parsing, setting, and reporting the usage of
command-line arguments common to all of our C/C++ demo applications, as well as
other methods that we commonly use across our C/C++ demo applications, e.g.,
standard methods for I/O handling. These source files are symlinks to the
original files contained in the repository
[Signaloid-Demo-CommonUtilityRoutines](https://github.com/signaloid/Signaloid-Demo-CommonUtilityRoutines)
which is included as a submodule in `submodules/common`.

### `uxhw.c/h`
These contain methods that implement the probabilistic versions of the methods
in the UxHw API (e.g., `UxHwDoubleGaussDist`) and use the GNU Scientific Library
(GSL) random number generators to achieve that. This allows building our C/C++
demo applications natively (i.e., on conventional architectures) and running
native Monte Carlo evaluations of our C/C++ demo applications without modifying
the source code. These source files are symlinks to the original files contained
in the repository
[Signaloid-Demo-UxHwCompatibilityForNativeExecution](https://github.com/signaloid/Signaloid-Demo-UxHwCompatibilityForNativeExecution)
which is included as a submodule in `submodules/compat`.

`uxhw.c` is excluded from `config.mk`: on Signaloid cores the UxHw API is
provided by the hardware, so the compatibility layer is only compiled into
native builds.

## config.mk
Signaloid cores use this file to identify the source codes they will use when
building the C/C++ demo application, along with the `CFLAGS` and `CXXFLAGS`
used for those builds.
