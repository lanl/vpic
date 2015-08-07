# Vector Particle-In-Cell (VPIC) Project

VPIC is a fully relativistic plasma simulation code...

# Build Instructions

VPIC uses the CMake build system.  To configure a build, do the following from
the top-level source directory:
  
   % mkdir build
	% cd build

Then call the curses version of CMake:

	% ccmake ..

or, optionally:

	% cmake -DENABLE_MPI ..

After configuration, simply type 'make'.
