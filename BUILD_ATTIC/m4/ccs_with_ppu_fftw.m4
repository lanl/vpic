AC_DEFUN([CCS_WITH_PPU_FFTW], [

    AC_ARG_VAR([PPU_FFTW],
        [User specified path to FFTW library top-level directory])
    AC_ARG_VAR([PPU_FFTW_CPPFLAGS],
		[User specified path to FFTW library includes])
    AC_ARG_VAR([PPU_FFTW_LDFLAGS], [User specified path to FFTW libraries])
    AC_ARG_VAR([PPU_FFTW_LIBS],
        [User specified link library names, e.g., -lfftw3])

    AC_ARG_WITH([ppu-fftw],
        [AC_HELP_STRING([--with-ppu-fftw],
        [User specified path to SPE libraries])],
        [
            if test -n "$PPU_FFTW" ; then
                if test -z "$PPU_FFTW_CPPFLAGS" ; then
                    ppu_fftw_CPPFLAGS="-I$PPU_FFTW/include"
                else
                    ppu_fftw_CPPFLAGS="$PPU_FFTW_CPPFLAGS"
                fi

                if test -z "$PPU_FFTW_LDFLAGS" ; then
                    ppu_fftw_LDFLAGS="-L$PPU_FFTW/lib"
                else
                    ppu_fftw_LDFLAGS="$PPU_FFTW_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$PPU_FFTW_CPPFLAGS" ; then
                    ppu_fftw_CPPFLAGS="-I$withval/include"
                else
                    ppu_fftw_CPPFLAGS="$PPU_FFTW_CPPFLAGS"
                fi

                if test -z "$PPU_FFTW_LDFLAGS" ; then
                    ppu_fftw_LDFLAGS="-L$withval/lib"
                else
                    ppu_fftw_LDFLAGS="$PPU_FFTW_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$PPU_FFTW" ; then
                if test -z "$PPU_FFTW_CPPFLAGS" ; then
                    ppu_fftw_CPPFLAGS="-I$PPU_FFTW/include"
                else
                    ppu_fftw_CPPFLAGS="$PPU_FFTW_CPPFLAGS"
                fi

                if test -z "$PPU_FFTW_LDFLAGS" ; then
                    ppu_fftw_LDFLAGS="-L$PPU_FFTW/lib"
                else
                    ppu_fftw_LDFLAGS="$PPU_FFTW_LDFLAGS"
                fi
            else
                if test -z "$PPU_FFTW_CPPFLAGS" ; then
                    ppu_fftw_CPPFLAGS="-I/usr/include"
                else
                    ppu_fftw_CPPFLAGS="$PPU_FFTW_CPPFLAGS"
                fi

                if test -z "$PPU_FFTW_LDFLAGS" ; then
                    ppu_fftw_LDFLAGS="-L/usr/lib"
                else
                    ppu_fftw_LDFLAGS="$PPU_FFTW_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $ppu_fftw_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $ppu_fftw_LDFLAGS"

        AC_CHECK_HEADER(fftw3.h, [ppu_fftw_h=yes], [ppu_fftw_h=no], /* check */)
        AC_CHECK_LIB(fftw3, fftw_transpose, [ppu_fftw_a=yes],
			[ppu_fftw_a=no], [])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$ppu_fftw_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid fftw3.h)
        else
            AC_SUBST(PPU_FFTW_CPPFLAGS, [$ppu_fftw_CPPFLAGS])
        fi

        if test "$ppu_fftw_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libfftw3.a)
        else
            AC_SUBST(PPU_FFTW_LDFLAGS, [$ppu_fftw_LDFLAGS])
        fi

        if test -z "$PPU_FFTW_LIBS" ; then
            AC_SUBST(PPU_FFTW_LIBS, ["-lfftw3"])
        fi
    ])
