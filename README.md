# Vector Particle-In-Cell (VPIC) Project

VPIC is a fully relativistic plasma simulation code...

# Getting the Code

VPIC uses nested submodules.  This requires the addition of the *--recursive*
flag when cloning the repository:

```
% git clone --recursive git@gitlab.lanl.gov:plasma/vpic.git
```

This command will check out the VPIC source code, including the Cinch
build system.  Cinch is documented
[here](http://gitlab.lanl.gov/ngc-utils/cinch).

# Build Instructions

VPIC uses the Cinch build system.  From a user-perspective, this is
equivalent to CMake.  To configure a build, do the following from
the top-level source directory:
  
```
% mkdir build
% cd build
```

Then call the curses version of CMake:

```
% ccmake ..
```

**or, optionally:**

```
% cmake -DENABLE_MPI ..
```

After configuration, simply type 'make'.
