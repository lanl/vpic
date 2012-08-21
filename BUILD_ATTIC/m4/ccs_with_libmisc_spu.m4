AC_DEFUN([CCS_WITH_LIBMISC_SPU], [
    AC_ARG_VAR([LIBMISC_SPU],
        [User specified path to SPE library top-level directory])
    AC_ARG_VAR([LIBMISC_SPU_CPPFLAGS],
        [User specified path to SPE library includes])
    AC_ARG_VAR([LIBMISC_SPU_LDFLAGS], [User specified path to SPE libraries])
    AC_ARG_VAR([LIBMISC_SPU_LIBS],
        [User specified link library names, e.g., -lspe -lpthreads])

    AC_ARG_WITH([libmisc_spu],
        [AC_HELP_STRING([--with-libmisc-spu],
        [User specified path to SPE libraries])],
        [
            if test -n "$LIBMISC_SPU" ; then
                if test -z "$LIBMISC_SPU_CPPFLAGS" ; then
                    libmisc_spu_CPPFLAGS="-I$LIBMISC_SPU/include"
                else
                    libmisc_spu_CPPFLAGS="$LIBMISC_SPU_CPPFLAGS"
                fi

                if test -z "$LIBMISC_SPU_LDFLAGS" ; then
                    libmisc_spu_LDFLAGS="-L$LIBMISC_SPU/lib"
                else
                    libmisc_spu_LDFLAGS="$LIBMISC_SPU_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$LIBMISC_SPU_CPPFLAGS" ; then
                    libmisc_spu_CPPFLAGS="-I$withval/include"
                else
                    libmisc_spu_CPPFLAGS="$LIBMISC_SPU_CPPFLAGS"
                fi

                if test -z "$LIBMISC_SPU_LDFLAGS" ; then
                    libmisc_spu_LDFLAGS="-L$withval/lib"
                else
                    libmisc_spu_LDFLAGS="$LIBMISC_SPU_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$LIBMISC_SPU" ; then
                if test -z "$LIBMISC_SPU_CPPFLAGS" ; then
                    libmisc_spu_CPPFLAGS="-I$LIBMISC_SPU/include"
                else
                    libmisc_spu_CPPFLAGS="$LIBMISC_SPU_CPPFLAGS"
                fi

                if test -z "$LIBMISC_SPU_LDFLAGS" ; then
                    libmisc_spu_LDFLAGS="-L$LIBMISC_SPU/lib"
                else
                    libmisc_spu_LDFLAGS="$LIBMISC_SPU_LDFLAGS"
                fi
            else
                libmisc_spu_CPPFLAGS="-I/opt/cell/sysroot/opt/cell/sdk/usr/spu/include"
                libmisc_spu_LDFLAGS="-L/opt/cell/sysroot/opt/cell/sdk/usr/spu/lib"
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $libmisc_spu_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $libmisc_spu_LDFLAGS"

        AC_CHECK_HEADER(libmisc.h, [libmisc_spu_h=yes],
            [libmisc_spu_h=no], /* check */)
        AC_CHECK_LIB(misc, malloc_align, [libmisc_spu_a=yes],
            [libmisc_spu_a=no])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$libmisc_spu_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmisc_spu.h)
        else
            AC_SUBST(LIBMISC_SPU_CPPFLAGS, [$libmisc_spu_CPPFLAGS])
        fi

        if test "$libmisc_spu_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmisc_spu.a)
        else
            AC_SUBST(LIBMISC_SPU_LDFLAGS, [$libmisc_spu_LDFLAGS])
        fi

        if test -z "$LIBMISC_SPU_LIBS" ; then
            AC_SUBST(LIBMISC_SPU_LIBS, ["-lmisc"])
        fi
    ])
