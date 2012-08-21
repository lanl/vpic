AC_DEFUN([CCS_WITH_PPU_DACS], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([PPU_DACS],
        [User specified path to PPU DaCS library top-level directory])
    AC_ARG_VAR([PPU_DACS_CPPFLAGS], [User specified path to DaCS library includes])
    AC_ARG_VAR([PPU_DACS_LDFLAGS], [User specified path to DaCS libraries])
    AC_ARG_VAR([PPU_DACS_LIBS],
        [User specified link library names, e.g., -ldacs_hybrid])

    AC_ARG_WITH([ppu_dacs],
        [AC_HELP_STRING([--with-ppu_dacs],
        [User specified path to SPE libraries])],
        [
            if test -n "$PPU_DACS" ; then
                if test -z "$PPU_DACS_CPPFLAGS" ; then
                    ppu_dacs_CPPFLAGS="-I$PPU_DACS/include"
                else
                    ppu_dacs_CPPFLAGS="$PPU_DACS_CPPFLAGS"
                fi

                if test -z "$PPU_DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_dacs_LDFLAGS="-L$PPU_DACS/lib"
                    else
                        ppu_dacs_LDFLAGS="-L$PPU_DACS/lib64"
                    fi
                else
                    ppu_dacs_LDFLAGS="$PPU_DACS_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$PPU_DACS_CPPFLAGS" ; then
                    ppu_dacs_CPPFLAGS="-I$withval/include"
                else
                    ppu_dacs_CPPFLAGS="$PPU_DACS_CPPFLAGS"
                fi

                if test -z "$PPU_DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_dacs_LDFLAGS="-L$withval/lib"
                    else
                        ppu_dacs_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    ppu_dacs_LDFLAGS="$PPU_DACS_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$PPU_DACS" ; then
                if test -z "$PPU_DACS_CPPFLAGS" ; then
                    ppu_dacs_CPPFLAGS="-I$PPU_DACS/include"
                else
                    ppu_dacs_CPPFLAGS="$PPU_DACS_CPPFLAGS"
                fi

                if test -z "$PPU_DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_dacs_LDFLAGS="-L$PPU_DACS/lib"
                    else
                        ppu_dacs_LDFLAGS="-L$PPU_DACS/lib64"
                    fi
                else
                    ppu_dacs_LDFLAGS="$PPU_DACS_LDFLAGS"
                fi
            else

                if test -z "$PPU_DACS_CPPFLAGS" ; then
                    ppu_dacs_CPPFLAGS="-I/usr/include"
                else
                    ppu_dacs_CPPFLAGS="$PPU_DACS_CPPFLAGS"
                fi

                if test -z "$PPU_DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_dacs_LDFLAGS="-L/usr/lib"
                    else
                        ppu_dacs_LDFLAGS="-L/usr/lib64"
                    fi
                else
                    ppu_dacs_LDFLAGS="$PPU_DACS_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $ppu_dacs_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $ppu_dacs_LDFLAGS"

        AC_CHECK_HEADER(dacs.h, [ppu_dacs_h=yes], [ppu_dacs_h=no], /* check */)
        AC_CHECK_LIB(dacs$1, dacs_send, [ppu_dacs_a=yes],
            [ppu_dacs_a=no], [-lstdc++ -lrt -lpthread])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$ppu_dacs_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid dacs.h)
        else
            AC_SUBST(PPU_DACS_CPPFLAGS, [$ppu_dacs_CPPFLAGS])
        fi

        if test "$ppu_dacs_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libdacs$1.a)
        else
            AC_SUBST(PPU_DACS_LDFLAGS, [$ppu_dacs_LDFLAGS])
        fi

        if test -z "$PPU_DACS_LIBS" ; then
            AC_SUBST(PPU_DACS_LIBS, ["-ldacs$1 -lstdc++ -lrt -lpthread"])
        fi
    ])
