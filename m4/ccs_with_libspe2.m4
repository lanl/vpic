AC_DEFUN([CCS_WITH_LIBSPE2], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([LIBSPE2],
        [User specified path to SPE library top-level directory])
    AC_ARG_VAR([LIBSPE2_CPPFLAGS], [User specified path to SPE library includes])
    AC_ARG_VAR([LIBSPE2_LDFLAGS], [User specified path to SPE libraries])
    AC_ARG_VAR([LIBSPE2_LIBS],
        [User specified link library names, e.g., -lspe -lpthreads])

    AC_ARG_WITH([libspe2],
        [AC_HELP_STRING([--with-libspe2],
        [User specified path to SPE libraries])],
        [
            if test -n "$LIBSPE2" ; then
                if test -z "$LIBSPE2_CPPFLAGS" ; then
                    libspe2_CPPFLAGS="-I$LIBSPE2/include"
                else
                    libspe2_CPPFLAGS="$LIBSPE2_CPPFLAGS"
                fi

                if test -z "$LIBSPE2_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe2_LDFLAGS="-L$LIBSPE2/lib"
                    else
                        libspe2_LDFLAGS="-L$LIBSPE2/lib64"
                    fi
                else
                    libspe2_LDFLAGS="$LIBSPE2_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$LIBSPE2_CPPFLAGS" ; then
                    libspe2_CPPFLAGS="-I$withval/include"
                else
                    libspe2_CPPFLAGS="$LIBSPE2_CPPFLAGS"
                fi

                if test -z "$LIBSPE2_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe2_LDFLAGS="-L$withval/lib"
                    else
                        libspe2_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    libspe2_LDFLAGS="$LIBSPE2_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$LIBSPE2" ; then
                if test -z "$LIBSPE2_CPPFLAGS" ; then
                    libspe2_CPPFLAGS="-I$LIBSPE2/include"
                else
                    libspe2_CPPFLAGS="$LIBSPE2_CPPFLAGS"
                fi

                if test -z "$LIBSPE2_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe2_LDFLAGS="-L$LIBSPE2/lib"
                    else
                        libspe2_LDFLAGS="-L$LIBSPE2/lib64"
                    fi
                else
                    libspe2_LDFLAGS="$LIBSPE2_LDFLAGS"
                fi
            else
                libspe2_CPPFLAGS="-I/opt/cell/sysroot/usr/include"
                if test "$BITS" = "32" ; then
                    libspe2_LDFLAGS="-L/opt/cell/sysroot/usr/lib"
                else
                    libspe2_LDFLAGS="-L/opt/cell/sysroot/usr/lib64"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $libspe2_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $libspe2_LDFLAGS"

        AC_CHECK_HEADER(libspe2.h, [libspe2_h=yes], [libspe2_h=no], /* check */)
        AC_CHECK_LIB(spe2, spe_context_create, [libspe2_a=yes],
            [libspe2_a=no], -lpthread)

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$libspe2_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid libspe2.h)
        else
            AC_SUBST(LIBSPE2_CPPFLAGS, [$libspe2_CPPFLAGS])
        fi

        if test "$libspe2_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libspe2.a)
        else
            AC_SUBST(LIBSPE2_LDFLAGS, [$libspe2_LDFLAGS])
        fi

        if test -z "$LIBSPE2_LIBS" ; then
            AC_SUBST(LIBSPE2_LIBS, ["-lspe2 -lpthread"])
        fi
    ])
