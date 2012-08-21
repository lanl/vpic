include config/basic.mk

BUILDDIR=core2

CC+=-fomit-frame-pointer -march=core2 -mfpmath=sse

# V4's SSE implementation breaks some aliasing rules
CXX+=-fomit-frame-pointer -march=core2 -mfpmath=sse -fno-strict-aliasing -DUSE_V4_SSE

