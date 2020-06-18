# Vector Particle-In-Cell (VPIC) Project

VPIC is a general purpose particle-in-cell simulation code for modeling
kinetic plasmas in one, two, or three spatial dimensions. It employs a
second-order, explicit, leapfrog algorithm to update charged particle
positions and velocities in order to solve the relativistic kinetic
equation for each species in the plasma, along with a full Maxwell
description for the electric and magnetic fields evolved via a second-
order finite-difference-time-domain (FDTD) solve. The VPIC code has been
optimized for modern computing architectures and uses Message Passing
Interface (MPI) calls for multi-node application as well as data
parallelism using threads. VPIC employs a variety of short-vector,
single-instruction-multiple-data (SIMD) intrinsics for high performance
and has been designed so that the data structures align with cache
boundaries. The current feature set for VPIC includes a flexible input
deck format capable of treating a wide variety of problems. These
include: the ability to treat electromagnetic materials (scalar and
tensor dielectric, conductivity, and diamagnetic material properties);
multiple emission models, including user-configurable models; arbitrary,
user-configurable boundary conditions for particles and fields; user-
definable simulation units; a suite of "standard" diagnostics, as well
as user-configurable diagnostics; a Monte-Carlo treatment of collisional
processes capable of treating binary and unary collisions and secondary
particle generation; and, flexible checkpoint-restart semantics enabling
VPIC checkpoint files to be read as input for subsequent simulations.
VPIC has a native I/O format that interfaces with the high-performance
visualization software Ensight and Paraview. While the common use cases
for VPIC employ low-order particles on rectilinear meshes, a framework
exists to treat higher-order particles and curvilinear meshes, as well
as more advanced field solvers.

# Attribution

Researchers who use the VPIC code for scientific research are asked to cite
the papers by Kevin Bowers listed below.

1. Bowers, K. J., B. J. Albright, B. Bergen, L. Yin, K. J. Barker and
D. J. Kerbyson, "0.374 Pflop/s Trillion-Particle Kinetic Modeling of
Laser Plasma Interaction on Road-runner," Proc. 2008 ACM/IEEE Conf.
Supercomputing (Gordon Bell Prize Finalist Paper).
http://dl.acm.org/citation.cfm?id=1413435

2. K.J. Bowers, B.J. Albright, B. Bergen and T.J.T. Kwan, Ultrahigh
performance three-dimensional electromagnetic relativistic kinetic
plasma simulation, Phys. Plasmas 15, 055703 (2008);
http://dx.doi.org/10.1063/1.2840133

3. K.J. Bowers, B.J. Albright, L. Yin, W. Daughton, V. Roytershteyn,
B. Bergen and T.J.T Kwan, Advances in petascale kinetic simulations
with VPIC and Roadrunner, Journal of Physics: Conference Series 180,
012055, 2009

# Getting the Code

To checkout the VPIC source, do the following:

```bash
    git clone https://github.com/lanl/vpic.git
```

## Branches

The stable release of vpic exists on `master`, the default branch.

For more cutting edge features, consider using the `devel` branch.

User contributions should target the `devel` branch.

# Requirements

The primary requirement to build VPIC is a C++11 capable compiler and
an up-to-date version of MPI.

# Build Instructions

```bash
    cd vpic 
```

VPIC uses the CMake build system. To configure a build, do the following from
the top-level source directory:
  
```bash
    mkdir build
    cd build
```

The `./arch` directory also contains various cmake scripts (including specific build options) which can help with building, but the user is left to select which compiler they wish to use.  The scripts are largely organized into folders by compiler, with specific flags and options set to match the target compiler.

Any of the arch scripts can be invoked specifying the file name from inside a build directory:

```bash
    ../arch/reference-Debug
```

After configuration, simply type: 

```bash
    make
```

Three scripts in the `./arch` directory are of particular note: lanl-ats1-hsw, lanl-ats1-knl and lanl-cts1. These scripts
provide a default way to build VPIC on LANL ATS-1 clusters such as Trinity and Trinitite and LANL CTS-1 clusters. The LANL
ATS-1 clusters are the first generation of DOE Advanced Technology Systems and consist of a partition of dual socket Intel
Haswell nodes and a partition of single socket Intel Knights Landing nodes. The LANL CTS-1 clusters are the first generation
of DOE Commodity Technology Systems and consist of dual socket Intel Broadwell nodes running the TOSS 3.3 operating system.
The lanl-ats1-hsw, lanl-ats1-knl and lanl-cts1 scripts are heavily documented and can be configured to provide a large
variety of custom builds for their respective platform types. These scripts could also serve as a good starting point for
development of a build script for other platform types. Because these scripts also configure the users build environment
via the use of module commands, the scripts run both the cmake and make commands.

From the user created build directory, these scripts can be invoked as follows:

```bash
    ../arch/lanl-ats1-hsw
```

or

```bash
    ../arch/lanl-ats1-knl
```

or

```bash
    ../arch/lanl-cts1
```

Advanced users may choose to instead invoke `cmake` directly and hand select options. Documentation on valid ways
to select these options may be found in the lanl-ats1 and lanl-cts1 build scripts mentioned above.

GCC users should ensure the `-fno-strict-aliasing` compiler flag is set (as shown in `./arch/generic-gcc-sse`).


# Building an example input deck

After you have successfully built VPIC, you should have an executable in
the `bin` directory called `vpic` (`./bin/vpic`).  To build an executable from one of
the sample input decks (found in `./sample`), simply run:

```bash
    ./bin/vpic input_deck
```

where *input_deck* is the name of your sample deck.  For example, to build
the *harris* input deck in the *sample* subdirectory
*(assuming that your build directory is located in the top-level
source directory)*:

```bash
    ./bin/vpic ../sample/harris
```

Beginners are advised to read the harris deck thoroughly, as it provides many examples of common uses cases.

# Command Line Arguments

Note: Historic VPIC users should note that the format of command line arguments was changed in the first open source release. The equals symbol is no longer accepted, and two dashes are mandatory. 

In general, command line arguments take the form `--command value`, in which two dashes are followed by a keyword, with a space delimiting the command and the value.

The following specific syntax is available to the users:

## Threading

Threading (per MPI rank) can be enabled using the following syntax: 

```bash
    ./binary.Linux --tpp n
```

Where n specifies the number of threads

### Example:

```bash
    mpirun -n 2 ./binary.Linux --tpp 2
```


To run with VPIC with two threads per MPI rank.

## Checkpoint Restart

VPIC can restart from a checkpoint dump file, using the following syntax:

```bash
    ./binary.Linux --restore <path to file>
```

### Example:

```bash
    ./binary.Linux --restore ./restart/restart0 
```

To restart VPIC using the restart file `./restart/restart0`

# Compile Time Arguments

Currently, the following options are exposed at compile time for the users consideration:

## Particle Array Resizing

- `DISABLE_DYNAMIC_RESIZING` (default `OFF`): Enable to disable the use of dynamic particle resizing
- `SET_MIN_NUM_PARTICLES` (default 128 [4kb]): Set the minimum number of particles allowable when dynamically resizing

## Threading Model

 - `USE_PTHREADS`: Use Pthreads for threading model, (default `ON`)
 - `USE_OPENMP`:   Use OpenMP for threading model

## Vectorization

The following CMake variables are used to control the vector implementation that
VPIC uses for each SIMD width.  Currently, there is support for 128 bit, 256 bit
and 512 bit SIMD widths.  The default is for each of these CMake variables to be
disabled which means that an unvectorized reference implementation of functions
will be used.

 - `USE_V4_SSE`:       Enable 4 wide (128-bit) SSE
 - `USE_V4_AVX`:       Enable 4 wide (128-bit) AVX
 - `USE_V4_AVX2`:      Enable 4 wide (128-bit) AVX2
 - `USE_V4_ALTIVEC`:   Enable 4 wide (128-bit) Altivec
 - `USE_V4_PORTABLE`:  Enable 4 wide (128-bit) portable implementation

 - `USE_V8_AVX`:       Enable 8 wide (256-bit) AVX
 - `USE_V8_AVX2`:      Enable 8 wide (256-bit) AVX2
 - `USE_V8_PORTABLE`:  Enable 8 wide (256-bit) portable implementation

 - `USE_V16_AVX512`:   Enable 16 wide (512-bit) AVX512
 - `USE_V16_PORTABLE`: Enable 16 wide (512-bit) portable implementation

Several functions in VPIC have vector implementations for each of the three SIMD
widths.  Some only have a single implementation.  An example of the latter is
move_p which only has a reference implementation and a V4 implementation.

It is possible to have a single CMake vector variable configured as ON for each
of the three supported SIMD vector widths.  It is recommended to always have a
CMake variable configured as ON for the 128 bit SIMD vector width so that move_p
will be vectorized.  In addition, it is recommended to configure as ON the CMake
variable that is associated with the native SIMD vector width of the processor
that VPIC is targeting.  If a CMake variable is configured as ON for each of the
three available SIMD vector widths, then for a given function in VPIC, the
implementation which supports the largest SIMD vector length will be chosen.  If
a V16 implementation exists, it will be chosen.  If a V16 implementation does not
exist but V8 and V4 implementations exist, the V8 implementation will be chosen.
If V16 and V8 implementations do not exist but a V4 implementation does, it will
be chosen.  If no SIMD vector implementation exists, the unvectorized reference
implementation will be chosen.

In summary, when using vector versions on a machine with 256 bit SIMD, the
V4 and V8 implementations should be configured as ON. When using a machine
with 512 bit SIMD, V4 and V16 implementations should be configured as ON.
When choosing a vector implementation for a given SIMD vector length, the
implementation that is closest to the SIMD instruction set for the targeted
processor should be chosen.  The portable versions are most commonly used for
debugging the implementation of new intrinsics versions.  However, the portable
versions are generally more performant than the unvectorized reference
implemenation.  So, one might consider using the V4_PORTABLE version on ARM
processors until a V4_NEON implementation becomes available.

## Output 

 - `VPIC_PRINT_MORE_DIGITS`: Enable more digits in timing output of status reports

## Particle sorting implementation

The CMake variable below allows building VPIC to use the legacy, thread serial
implementation of the particle sort algorithm.

 - `USE_LEGACY_SORT`: Use legacy thread serial particle sort, (default `OFF`)

The legacy particle sort implementation is the thread serial particle sort
implementation from the legacy v407 version of VPIC. This implementation
supports both in-place and out-of-place sorting of the particles. It is very
competitive with the thread parallel sort implementation for a small number
of threads per MPI rank, i.e. 4 or less, especially on KNL because sorting
the particles in-place allows the fraction of particles stored in High
Bandwidth Memory (HBM) to remain stored in HBM. Also, the memory footprint
of VPIC is reduced by the memory of a particle array which can be significant
for particle dominated problems.

The default particle sort implementation is a thread parallel implementation.
Currently, it can only perform out-of-place sorting of the particles. It will
be more performant than the legacy implementation when using many threads per
MPI rank but uses more memory because of the out-of-place sort.

# Workflow

Contributors are asked to be aware of the following workflow:

1) Pull requests are accepted into `devel` upon tests passing
2) `master` should reflect the *stable* state of the code
3) Periodic releases will be made from `devel` into `master`

# Feedback

Feedback, comments, or issues can be raised through [GitHub issues](https://github.com/lanl/vpic/issues).

A mailing list for open collaboration can also be found [here](https://groups.google.com/forum/#!forum/vpic-users)

# Versioning

Version release summary: 

## V1.1 (March 2019)

- Added V8 and V16 functionality
- Improved documentation and build processes
- Significantly improved testing and correctness capabilities

## V1.0

Initial release


# Release

This software has been approved for open source release and has been assigned **LA-CC-15-109**.

# Copyright

Copyright Triad National Security, LLC. All rights reserved.

This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC (TRIAD) for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of Triad National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
