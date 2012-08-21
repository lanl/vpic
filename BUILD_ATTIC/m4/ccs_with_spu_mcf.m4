AC_DEFUN([CCS_WITH_SPU_MCF], [
    AC_REQUIRE([CCS_WITH_MCF])

    AC_ARG_VAR([MCF_SPU_CPPFLAGS],
        [User specified path to MCF SPU library includes])
    AC_ARG_VAR([MCF_SPU_LDFLAGS], [User specified path to MCF SPU libraries])
    AC_ARG_VAR([MCF_SPU_LIBS],
        [User specified MCF SPU link library names, e.g., -lmcf_w -lcsal])
    AC_ARG_VAR([MCF_SPU_MAIN],
        [User specified MCF SPU mcf_w_main location])

    AC_ARG_WITH([ppu-mcf],
        [AC_HELP_STRING([--with-ppu-mcf],
        [User specified path to MCF SPU libraries])],
        [
            if test -n "$MCF" ; then
                if test -z "$MCF_SPU_CPPFLAGS" ; then
                    mcf_SPU_CPPFLAGS="-I$MCF/include"
                else
                    mcf_SPU_CPPFLAGS="$MCF_SPU_CPPFLAGS"
                fi

                if test -z "$MCF_SPU_LDFLAGS" ; then
                    mcf_SPU_LDFLAGS="-L$MCF/lib/spe"
                else
                    mcf_SPU_LDFLAGS="$MCF_SPU_LDFLAGS"
                fi

                if test -z "$MCF_SPU_MAIN" ; then
                    mcf_SPU_MAIN="$MCF/lib/spe/mcf_w_main.o"
                else
                    mcf_SPU_MAIN="$MCF_SPU_MAIN"
                fi
            elif test "$withval" != no ; then
                if test -z "$MCF_SPU_CPPFLAGS" ; then
                    mcf_SPU_CPPFLAGS="-I$withval/include"
                else
                    mcf_SPU_CPPFLAGS="$MCF_SPU_CPPFLAGS"
                fi

                if test -z "$MCF_SPU_LDFLAGS" ; then
                    mcf_SPU_LDFLAGS="-L$withval/lib/spe"
                else
                    mcf_SPU_LDFLAGS="$MCF_SPU_LDFLAGS"
                fi

                if test -z "$MCF_SPU_MAIN" ; then
                    mcf_SPU_MAIN="$withval/lib/spe/mcf_w_main.o"
                else
                    mcf_SPU_MAIN="$MCF_SPU_MAIN"
                fi
            fi
        ],
        [
            if test -n "$MCF" ; then
                if test -z "$MCF_SPU_CPPFLAGS" ; then
                    mcf_SPU_CPPFLAGS="-I$MCF/include"
                else
                    mcf_SPU_CPPFLAGS="$MCF_SPU_CPPFLAGS"
                fi

                if test -z "$MCF_SPU_LDFLAGS" ; then
                    mcf_SPU_LDFLAGS="-L$MCF/lib/spe"
                else
                    mcf_SPU_LDFLAGS="$MCF_SPU_LDFLAGS"
                fi

                if test -z "$MCF_SPU_MAIN" ; then
                    mcf_SPU_MAIN="$MCF/lib/spe/mcf_w_main.o"
                else
                    mcf_SPU_MAIN="$MCF_SPU_MAIN"
                fi
            else
                mcf_SPU_CPPFLAGS="-I/opt/MultiCorePlus/include"
                mcf_SPU_LDFLAGS="-L/opt/MultiCorePlus/lib/spe"
                mcf_SPU_MAIN="/opt/MultiCorePlus/lib/spe/mcf_w_main.o"
            fi
        ])

        if test -z "$MCF_SPU_LIBS" ; then
            mcf_SPU_LIBS="-lsal -lcsal -ltatl -lmcf_w"
        else
            mcf_SPU_LIBS="$MCF_SPU_LIBS"
        fi

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $mcf_SPU_CPPFLAGS -nostartfiles"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $mcf_SPU_LDFLAGS $mcf_SPU_MAIN"
        old_LIBS=$LIBS
        LIBS="$LIBS $mcf_SPU_LIBS"

        AC_CHECK_HEADER(mcf_w.h, [mcf_h=yes], [mcf_h=no], /* check */)

        AC_MSG_CHECKING(for libmcf_w.a)
        AC_LINK_IFELSE([[
#if defined __cplusplus
extern "C"
#endif

#include <mcf_w.h>

int mcf_w_main(int n_bytes, void * p_arg_ls) {
    int rc, net_rank;
    rc = mcf_w_net_get_rank(&net_rank);
    return 0;
} // main
        ]], [mcf_a=yes], [mcf_a=no])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"
        LIBS="$old_LIBS"

        if test "$mcf_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid mcf_m.h)
        else
            AC_SUBST(MCF_SPU_CPPFLAGS, [$mcf_SPU_CPPFLAGS])
        fi

        if test "$mcf_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmcf_w.a)
        else
            AC_SUBST(MCF_SPU_LDFLAGS, [$mcf_SPU_LDFLAGS])
            AC_MSG_RESULT(yes)
        fi

        if test -z "$MCF_SPU_LIBS" ; then
            AC_SUBST(MCF_SPU_LIBS, [$mcf_SPU_LIBS])
        fi
    ])
