AC_DEFUN([CCS_WITH_SPULIBS], [
    AC_ARG_VAR([SPULIBS],
        [User specified path to SPE library top-level directory])
    AC_ARG_VAR([SPULIBS_CPPFLAGS], [User specified path to SPE library includes])
    AC_ARG_VAR([SPULIBS_LDFLAGS], [User specified path to SPE libraries])
    AC_ARG_VAR([SPULIBS_LIBS],
        [User specified link library names, e.g., -llibc])

    AC_ARG_WITH([spulibs],
        [AC_HELP_STRING([--with-spulibs],
        [User specified path to SPE libraries])],
        [
            if test -n "$SPULIBS" ; then
                if test -z "$SPULIBS_CPPFLAGS" ; then
                    spulibs_CPPFLAGS="-I$SPULIBS/include"
                else
                    spulibs_CPPFLAGS="$SPULIBS_CPPFLAGS"
                fi

                if test -z "$SPULIBS_LDFLAGS" ; then
                    spulibs_LDFLAGS="-L$SPULIBS/lib"
                else
                    spulibs_LDFLAGS="$SPULIBS_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$SPULIBS_CPPFLAGS" ; then
                    spulibs_CPPFLAGS="-I$withval/include"
                else
                    spulibs_CPPFLAGS="$SPULIBS_CPPFLAGS"
                fi

                if test -z "$SPULIBS_LDFLAGS" ; then
                    spulibs_LDFLAGS="-L$withval/lib"
                else
                    spulibs_LDFLAGS="$SPULIBS_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$SPULIBS" ; then
                if test -z "$SPULIBS_CPPFLAGS" ; then
                    spulibs_CPPFLAGS="-I$SPULIBS/include"
                else
                    spulibs_CPPFLAGS="$SPULIBS_CPPFLAGS"
                fi

                if test -z "$SPULIBS_LDFLAGS" ; then
                    spulibs_LDFLAGS="-L$SPULIBS/lib"
                else
                    spulibs_LDFLAGS="$SPULIBS_LDFLAGS"
                fi
            else
                spulibs_CPPFLAGS="-I/opt/cell/sysroot/opt/cell/sdk/usr/spu/include"
                spulibs_LDFLAGS="-I/opt/cell/sysroot/opt/cell/sdk/usr/spu/lib"
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $spulibs_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $spulibs_LDFLAGS"

dnl Note:
dnl spu_ext.h seems to be the only file that tests cleanly
        AC_CHECK_HEADER(spu_ext.h,
            [spulibs_h=yes], [spulibs_h=no], /* check */)

dnl FIXME: this needs to be updated with the next sdk drop
dnl     AC_CHECK_LIB(c, printf, [spulibs_a=yes], [spulibs_a=no])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$spulibs_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid spu_ext.h)
        else
            AC_SUBST(SPULIBS_CPPFLAGS, [$spulibs_CPPFLAGS])
        fi

dnl     if test "$spulibs_a" != yes ; then
dnl         AC_MSG_ERROR(Failed to find valid spulibs.a)
dnl     else
dnl         AC_SUBST(SPULIBS_LDFLAGS, [$spulibs_LDFLAGS])
dnl     fi
    ])
