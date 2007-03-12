################################################################################
# bin Makefile.am
################################################################################

bin_PROGRAMS = advance accel cyclo inbndj interp_e outbndj pcomm \
	rng_conv rng_seq v4

INCLUDES = -I@top_srcdir@/src/include \
	@EXTRA_CPPFLAGS@

advance_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
advance_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/bench/advance_p
advance_LDFLAGS = @EXTRA_LDFLAGS@
advance_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

accel_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
accel_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/accel
accel_LDFLAGS = @EXTRA_LDFLAGS@
accel_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

cyclo_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
cyclo_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/cyclo
cyclo_LDFLAGS = @EXTRA_LDFLAGS@
cyclo_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

inbndj_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
inbndj_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/inbndj
inbndj_LDFLAGS = @EXTRA_LDFLAGS@
inbndj_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

interp_e_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
interp_e_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/interp_e
interp_e_LDFLAGS = @EXTRA_LDFLAGS@
interp_e_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

outbndj_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
outbndj_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/outbndj
outbndj_LDFLAGS = @EXTRA_LDFLAGS@
outbndj_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

pcomm_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
pcomm_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/pcomm
pcomm_LDFLAGS = @EXTRA_LDFLAGS@
pcomm_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

rng_conv_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
rng_conv_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/rng_conv
rng_conv_LDFLAGS = @EXTRA_LDFLAGS@
rng_conv_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

rng_seq_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
rng_seq_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/rng_seq
rng_seq_LDFLAGS = @EXTRA_LDFLAGS@
rng_seq_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@

v4_SOURCES = @top_srcdir@/src/main.cpp \
	@top_srcdir@/src/deck_wrapper.cpp
v4_CPPFLAGS = -DINPUT_DECK=@top_srcdir@/sample/test/v4
v4_LDFLAGS = @EXTRA_LDFLAGS@
v4_LDADD = @top_builddir@/lib/libvpic.la @EXTRA_LIBS@
