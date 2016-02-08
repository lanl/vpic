# Vector Particle-In-Cell (VPIC) Project

VPIC is a fully relativistic plasma simulation code...

# Getting the Code

VPIC uses nested submodules.  This requires the addition of the *--recursive*
flag when cloning the repository:

    % git clone --recursive git@github.com:vpic/vpic.git

This command will check out the VPIC source code, including the Cinch
build system.  Cinch is documented
[here](https://github.com/flecsi/cinch).

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
