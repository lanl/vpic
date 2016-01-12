AC_DEFUN([CCS_WITH_LIBSPE], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([LIBSPE],
        [User specified path to SPE library top-level directory])
    AC_ARG_VAR([LIBSPE_CPPFLAGS], [User specified path to SPE library includes])
    AC_ARG_VAR([LIBSPE_LDFLAGS], [User specified path to SPE libraries])
    AC_ARG_VAR([LIBSPE_LIBS],
        [User specified link library names, e.g., -lspe -lpthreads])

    AC_ARG_WITH([libspe],
        [AC_HELP_STRING([--with-libspe],
        [User specified path to SPE libraries])],
        [
            if test -n "$LIBSPE" ; then
                if test -z "$LIBSPE_CPPFLAGS" ; then
                    libspe_CPPFLAGS="-I$LIBSPE/include"
                else
                    libspe_CPPFLAGS="$LIBSPE_CPPFLAGS"
                fi

                if test -z "$LIBSPE_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe_LDFLAGS="-L$LIBSPE/lib"
                    else
                        libspe_LDFLAGS="-L$LIBSPE/lib64"
                    fi
                else
                    libspe_LDFLAGS="$LIBSPE_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$LIBSPE_CPPFLAGS" ; then
                    libspe_CPPFLAGS="-I$withval/include"
                else
                    libspe_CPPFLAGS="$LIBSPE_CPPFLAGS"
                fi

                if test -z "$LIBSPE_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe_LDFLAGS="-L$withval/lib"
                    else
                        libspe_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    libspe_LDFLAGS="$LIBSPE_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$LIBSPE" ; then
                if test -z "$LIBSPE_CPPFLAGS" ; then
                    libspe_CPPFLAGS="-I$LIBSPE/include"
                else
                    libspe_CPPFLAGS="$LIBSPE_CPPFLAGS"
                fi

                if test -z "$LIBSPE_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        libspe_LDFLAGS="-L$LIBSPE/lib"
                    else
                        libspe_LDFLAGS="-L$LIBSPE/lib64"
                    fi
                else
                    libspe_LDFLAGS="$LIBSPE_LDFLAGS"
                fi
            else
                libspe_CPPFLAGS="-I/opt/IBM/cell-sdk-1.1/sysroot/usr/include"
                if test "$BITS" = "32" ; then
                    libspe_LDFLAGS="-L/opt/IBM/cell-sdk-1.1/sysroot/usr/lib"
                else
                    libspe_LDFLAGS="-L/opt/IBM/cell-sdk-1.1/sysroot/usr/lib64"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $libspe_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $libspe_LDFLAGS"

        AC_CHECK_HEADER(libspe.h, [libspe_h=yes], [libspe_h=no], /* check */)
        AC_CHECK_LIB(spe, spe_create_thread, [libspe_a=yes],
            [libspe_a=no], -lpthread)

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$libspe_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid libspe.h)
        else
            AC_SUBST(LIBSPE_CPPFLAGS, [$libspe_CPPFLAGS])
        fi

        if test "$libspe_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libspe.a)
        else
            AC_SUBST(LIBSPE_LDFLAGS, [$libspe_LDFLAGS])
        fi

        if test -z "$LIBSPE_LIBS" ; then
            AC_SUBST(LIBSPE_LIBS, ["-lspe -lpthread"])
        fi
    ])
