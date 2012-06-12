AC_DEFUN([CCS_SIMD_ISA], [
    AC_REQUIRE([AC_PROG_CC])
    AC_REQUIRE([AC_PROG_CXX])
    AC_REQUIRE([CCS_ENABLE_SCALAR])
    AC_REQUIRE([CCS_WITH_OPT_LEVEL])


    if test "$ENABLE_SIMD" = "1" ; then
        AC_MSG_CHECKING(isa for vectorization)

        cc=`echo $CC | sed 's%^.*[[/$]]%%g'`
        cxx=`echo $CXX | sed 's%^.*[[/$]]%%g'`
        
        case "$cc" in
            *gcc*)
                machine=`$CC -dumpmachine | sed 's/-.*//g'`
                case "$machine" in
                    i?86)
                        CFLAGS="$CFLAGS -msse -msse2 -mfpmath=sse"
                        AC_SUBST(ISA, [SSE])
                        ISA_MSG="SSE"
                    ;;
                    x86_64)
                        CFLAGS="$CFLAGS -msse -msse2 -mfpmath=sse"
                        AC_SUBST(ISA, [SSE])
                        ISA_MSG="SSE"
                    ;;
                    powerpc | ppu)
                        CFLAGS="$CFLAGS -mabi=altivec -maltivec -include altivec.h"
                        AC_SUBST(ISA, [ALTIVEC])
                        ISA_MSG="altivec"
                    ;;
                    spu)
                        AC_SUBST(ISA, [SPU])
                        ISA_MSG="SPU"
                    ;;
                esac
            ;;
            icc)
                CFLAGS="$CFLAGS -msse3"
            ;;
            *)
                AC_MSG_ERROR(compiler not supported)
            ;;
        esac

        case "$cxx" in
            *g++*)
                machine=`$CXX -dumpmachine | sed 's/-.*//g'`
                case "$machine" in
                    i?86)
                        CXXFLAGS="$CXXFLAGS -msse -msse2 -mfpmath=sse"
                    ;;
                    x86_64)
                        CXXFLAGS="$CXXFLAGS -msse -msse2 -mfpmath=sse"
                    ;;
                    powerpc | ppu)
                        CXXFLAGS="$CXXFLAGS -mabi=altivec -maltivec -include altivec.h"
                    ;;
                esac
            ;;
            icpc)
                CXXFLAGS="$CXXFLAGS -msse3"
            ;;
            *)
                AC_MSG_ERROR(compiler not supported)
            ;;
        esac

        AC_MSG_RESULT($ISA_MSG)
    else
        AC_SUBST(ISA, [SCALAR])
    fi
])
