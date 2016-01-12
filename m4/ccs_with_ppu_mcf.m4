AC_DEFUN([CCS_WITH_PPU_MCF], [
    AC_REQUIRE([CCS_WITH_MCF])
    AC_REQUIRE([CCS_WITH_LIBSPE])
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([MCF_PPU_CPPFLAGS],
        [User specified path to MCF PPU library includes])
    AC_ARG_VAR([MCF_PPU_LDFLAGS], [User specified path to MCF PPU libraries])
    AC_ARG_VAR([MCF_PPU_LIBS],
        [User specified MCF PPU link library names, e.g., -lmcf_m -lcsal])

    AC_ARG_WITH([ppu-mcf],
        [AC_HELP_STRING([--with-ppu-mcf],
        [User specified path to MCF PPU libraries])],
        [
            if test -n "$MCF" ; then
                if test -z "$MCF_PPU_CPPFLAGS" ; then
                    mcf_PPU_CPPFLAGS="-I$MCF/include"
                else
                    mcf_PPU_CPPFLAGS="$MCF_PPU_CPPFLAGS"
                fi

                if test -z "$MCF_PPU_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mcf_PPU_LDFLAGS="-L$MCF/lib"
                    else
                        mcf_PPU_LDFLAGS="-L$MCF/lib64"
                    fi
                else
                    mcf_PPU_LDFLAGS="$MCF_PPU_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$MCF_PPU_CPPFLAGS" ; then
                    mcf_PPU_CPPFLAGS="-I$withval/include"
                else
                    mcf_PPU_CPPFLAGS="$MCF_PPU_CPPFLAGS"
                fi

                if test -z "$MCF_PPU_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mcf_PPU_LDFLAGS="-L$withval/lib"
                    else
                        mcf_PPU_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    mcf_PPU_LDFLAGS="$MCF_PPU_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$MCF" ; then
                if test -z "$MCF_PPU_CPPFLAGS" ; then
                    mcf_PPU_CPPFLAGS="-I$MCF/include"
                else
                    mcf_PPU_CPPFLAGS="$MCF_PPU_CPPFLAGS"
                fi

                if test -z "$MCF_PPU_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mcf_PPU_LDFLAGS="-L$MCF/lib"
                    else
                        mcf_PPU_LDFLAGS="-L$MCF/lib64"
                    fi
                else
                    mcf_PPU_LDFLAGS="$MCF_PPU_LDFLAGS"
                fi
            else
                mcf_PPU_CPPFLAGS="-I/opt/MultiCorePlus/include"
                if test "$BITS" = "32" ; then
                    mcf_PPU_LDFLAGS="-L/opt/MultiCorePlus/lib"
                else
                    mcf_PPU_LDFLAGS="-L/opt/MultiCorePlus/lib64"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $mcf_PPU_CPPFLAGS $LIBSPE_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $mcf_PPU_LDFLAGS $LIBSPE_LDFLAGS"

        AC_CHECK_HEADER(mcf_m.h, [mcf_h=yes], [mcf_h=no], /* check */)
        AC_CHECK_LIB(mcf_m, mcf_m_net_create, [mcf_a=yes],
            [mcf_a=no], [$LIBSPE_LIBS])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$mcf_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid mcf_m.h)
        else
            AC_SUBST(MCF_PPU_CPPFLAGS, [$mcf_PPU_CPPFLAGS])
        fi

        if test "$mcf_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmcf_m.a)
        else
            AC_SUBST(MCF_PPU_LDFLAGS, [$mcf_PPU_LDFLAGS])
        fi

        if test -z "$MCF_PPU_LIBS" ; then
            AC_SUBST(MCF_PPU_LIBS, ["-lsal -lcsal -ltatl -lmcf_m"])
        fi
    ])
