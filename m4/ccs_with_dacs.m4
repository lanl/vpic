AC_DEFUN([CCS_WITH_DACS], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([DACS],
        [User specified path to DaCS library top-level directory])
    AC_ARG_VAR([DACS_CPPFLAGS], [User specified path to DaCS library includes])
    AC_ARG_VAR([DACS_LDFLAGS], [User specified path to DaCS libraries])
    AC_ARG_VAR([DACS_LIBS],
        [User specified link library names, e.g., -ldacs_hybrid])

    AC_ARG_WITH([dacs],
        [AC_HELP_STRING([--with-dacs],
        [User specified path to SPE libraries])],
        [
            if test -n "$DACS" ; then
                if test -z "$DACS_CPPFLAGS" ; then
                    dacs_CPPFLAGS="-I$DACS/include"
                else
                    dacs_CPPFLAGS="$DACS_CPPFLAGS"
                fi

                if test -z "$DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        dacs_LDFLAGS="-L$DACS/lib"
                    else
                        dacs_LDFLAGS="-L$DACS/lib64"
                    fi
                else
                    dacs_LDFLAGS="$DACS_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$DACS_CPPFLAGS" ; then
                    dacs_CPPFLAGS="-I$withval/include"
                else
                    dacs_CPPFLAGS="$DACS_CPPFLAGS"
                fi

                if test -z "$DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        dacs_LDFLAGS="-L$withval/lib"
                    else
                        dacs_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    dacs_LDFLAGS="$DACS_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$DACS" ; then
                if test -z "$DACS_CPPFLAGS" ; then
                    dacs_CPPFLAGS="-I$DACS/include"
                else
                    dacs_CPPFLAGS="$DACS_CPPFLAGS"
                fi

                if test -z "$DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        dacs_LDFLAGS="-L$DACS/lib"
                    else
                        dacs_LDFLAGS="-L$DACS/lib64"
                    fi
                else
                    dacs_LDFLAGS="$DACS_LDFLAGS"
                fi
            else
                if test -z "$DACS_CPPFLAGS" ; then
                    dacs_CPPFLAGS="-I/usr/include"
                else
                    dacs_CPPFLAGS="$DACS_CPPFLAGS"
                fi

                if test -z "$DACS_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        dacs_LDFLAGS="-L/usr/lib"
                    else
                        dacs_LDFLAGS="-L/usr/lib64"
                    fi
                else
                    dacs_LDFLAGS="$DACS_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $dacs_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $dacs_LDFLAGS"

        AC_CHECK_HEADER(dacs.h, [dacs_h=yes], [dacs_h=no], /* check */)
        AC_CHECK_LIB(dacs_$1, dacs_init, [dacs_a=yes],
            [dacs_a=no], [-lstdc++ -lrt -lpthread])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$dacs_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid dacs.h)
        else
            AC_SUBST(DACS_CPPFLAGS, [$dacs_CPPFLAGS])
        fi

        if test "$dacs_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libdacs_$1.a)
        else
            AC_SUBST(DACS_LDFLAGS, [$dacs_LDFLAGS])
        fi

        if test -z "$DACS_LIBS" ; then
            AC_SUBST(DACS_LIBS, ["-ldacs_$1 -lstdc++ -lrt -lpthread"])
        fi
    ])
