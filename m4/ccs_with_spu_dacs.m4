AC_DEFUN([CCS_WITH_SPU_DACS], [

    AC_ARG_VAR([SPU_DACS],
        [User specified path to SPU DaCS library top-level directory])
    AC_ARG_VAR([SPU_DACS_CPPFLAGS], [User specified path to DaCS library includes])
    AC_ARG_VAR([SPU_DACS_LDFLAGS], [User specified path to DaCS libraries])
    AC_ARG_VAR([SPU_DACS_LIBS],
        [User specified link library names, e.g., -ldacs])

    AC_ARG_WITH([spu_dacs],
        [AC_HELP_STRING([--with-spu-dacs],
        [User specified path to SPE libraries])],
        [
            if test -n "$SPU_DACS" ; then
                if test -z "$SPU_DACS_CPPFLAGS" ; then
                    spu_dacs_CPPFLAGS="-I$SPU_DACS/include"
                else
                    spu_dacs_CPPFLAGS="$SPU_DACS_CPPFLAGS"
                fi

                if test -z "$SPU_DACS_LDFLAGS" ; then
                	spu_dacs_LDFLAGS="-L$SPU_DACS/lib"
                else
                    spu_dacs_LDFLAGS="$SPU_DACS_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$SPU_DACS_CPPFLAGS" ; then
                    spu_dacs_CPPFLAGS="-I$withval/include"
                else
                    spu_dacs_CPPFLAGS="$SPU_DACS_CPPFLAGS"
                fi

                if test -z "$SPU_DACS_LDFLAGS" ; then
                    spu_dacs_LDFLAGS="-L$withval/lib"
                else
                    spu_dacs_LDFLAGS="$SPU_DACS_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$SPU_DACS" ; then
                if test -z "$SPU_DACS_CPPFLAGS" ; then
                    spu_dacs_CPPFLAGS="-I$SPU_DACS/include"
                else
                    spu_dacs_CPPFLAGS="$SPU_DACS_CPPFLAGS"
                fi

                if test -z "$SPU_DACS_LDFLAGS" ; then
                    spu_dacs_LDFLAGS="-L$SPU_DACS/lib"
                else
                    spu_dacs_LDFLAGS="$SPU_DACS_LDFLAGS"
                fi
            else

                if test -z "$SPU_DACS_CPPFLAGS" ; then
                    spu_dacs_CPPFLAGS="-I/usr/include"
                else
                    spu_dacs_CPPFLAGS="$SPU_DACS_CPPFLAGS"
                fi

                if test -z "$SPU_DACS_LDFLAGS" ; then
                    spu_dacs_LDFLAGS="-L/usr/lib"
                else
                    spu_dacs_LDFLAGS="$SPU_DACS_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $spu_dacs_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $spu_dacs_LDFLAGS"

        AC_CHECK_HEADER(dacs.h, [spu_dacs_h=yes], [spu_dacs_h=no], /* check */)
        AC_CHECK_LIB(dacs, dacs_send, [spu_dacs_a=yes],
            [spu_dacs_a=no], [])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$spu_dacs_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid dacs.h)
        else
            AC_SUBST(SPU_DACS_CPPFLAGS, [$spu_dacs_CPPFLAGS])
        fi

        if test "$spu_dacs_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libdacs.a)
        else
            AC_SUBST(SPU_DACS_LDFLAGS, [$spu_dacs_LDFLAGS])
        fi

        if test -z "$SPU_DACS_LIBS" ; then
            AC_SUBST(SPU_DACS_LIBS, [-ldacs])
        fi
    ])
