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
parallelism using pthreads. VPIC employs a variety of short-vector,
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

VPIC uses nested submodules.  This requires the addition of the *--recursive*
flag when cloning the repository:

    % git clone --recursive git@github.com:lanl/vpic.git

This command will check out the VPIC source code, including the Cinch
build system.  Cinch is documented
[here](https://github.com/laristra/cinch).

# Requirements

The primary requirement to build VPIC is a C++11 capable compiler and
an up-to-date version of MPI.

# Build Instructions

VPIC uses the Cinch build system.  From a user-perspective, this is
equivalent to CMake.  To configure a build, do the following from
the top-level source directory:
  
    % mkdir build
    % cd build

Then call the curses version of CMake:

    % ccmake ..

**or, optionally:**

    % cmake -DENABLE_MPI ..

The `./arch` directory also contains various cmake scripts (including specific build options) which can help with building

They can be invoked using something like:

    % ../arch/generic-Release

GCC users should ensure the `-fno-strict-aliasing` compiler flag is set (as shown in `./arch/generic-gcc-sse`)

After configuration, simply type 'make'.

# Building an example input deck

After you have successfully built VPIC, you should have an executable in
the *bin* directory called *vpic*.  To build an executable from one of
the sample input decks, simply run:

    % bin/vpic input_deck

where *input_deck* is the name of your sample deck.  For example, to build
the *harris* input deck in the *sample* subdirectory
*(assuming that your build directory is located in the top-level
source directory)*:

    % bin/vpic ../sample/harris

# Release

This software has been approved for open source release and has been assigned **LA-CC-15-109**.

# Copyright

Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
