BUILDDIR=basic

# Basic compile settings

CC=mpicc -Werror -Wall -pedantic -std=c99   -D_XOPEN_SOURCE=600 -Wstrict-aliasing=2 -g -O2 -ffast-math -fno-unsafe-math-optimizations              

# -Wno-long-long is to accomodate some broken mpi headers under C++
CXX=mpiCC -Werror -Wall -pedantic -std=c++98 -D_XOPEN_SOURCE=600 -Wstrict-aliasing=2 -g -O2 -ffast-math -fno-unsafe-math-optimizations -Wno-long-long

LD=mpiCC -rdynamic

