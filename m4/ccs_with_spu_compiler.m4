AC_DEFUN([CCS_WITH_SPU_COMPILER], [
    AC_REQUIRE([CCS_ENABLE_STRIPPING])

    AC_ARG_VAR([SPU_COMPILER],
        [User specified SPU C/C++ compiler suite (gnu|xlc)])
    AC_ARG_VAR([SCC], [User specified SPU C compiler command])
    AC_ARG_VAR([SCFLAGS], [User specified SPU C compiler flags])
    AC_ARG_VAR([SCXX], [User specified SPU C++ compiler command])
    AC_ARG_VAR([SCXXFLAGS], [User specified SPU C++ compiler flags])
    AC_ARG_VAR([SLDFLAGS], [User specified SPU linker flags])
    AC_ARG_VAR([SCCAS], [User specified SPU assembler command])
    AC_ARG_VAR([SCCASFLAGS], [User specified SPU assembler flags])

    AC_ARG_WITH([spu-compiler],
        [AC_HELP_STRING([--with-spu-compiler],
        [User specified SPU C/C++ compiler suite (gnu|xlc)])],
        [

            if test -n "$SPU_COMPILER" ; then
                withval=$SPU_COMPILER
            fi

            if test "$withval" = "gnu" ; then
                if -n "$SCC" ; then
                    if test "$SCC" = "spuxlc" ; then
                        CCS_SPU_CC_COMPILER_XLC
                    else
                        CCS_SPU_CC_COMPILER_GNU
                    fi
                else
                    CCS_SPU_CC_COMPILER_GNU
                fi

                CCS_SPU_CXX_COMPILER_GNU
            elif test "$withval" = "xlc" ; then
                if -n "$SCC" ; then
                    if test "$SCC" = "spuxlc" ; then
                        CCS_SPU_CC_COMPILER_XLC
                    else
                        CCS_SPU_CC_COMPILER_GNU
                    fi
                else
                CCS_SPU_CC_COMPILER_XLC
                fi

                CCS_SPU_CXX_COMPILER_GNU
            fi

        ],
        [

            if test -n "$SCC" ; then
                if test "$SCC" = "spuxlc" ; then
                    CCS_SPU_CC_COMPILER_XLC
                else
                    CCS_SPU_CC_COMPILER_GNU
                fi
                CCS_SPU_CXX_COMPILER_GNU
            else
                CCS_SPU_CC_COMPILER_GNU
                CCS_SPU_CXX_COMPILER_GNU
            fi

        ])

    if test -z "$SLDFLAGS" ; then
        sldflags="-Wl,-N"
        if test "$ENABLE_STRIPPING" = "1" ; then
            pldflags="$pldflags -s"
        fi

        AC_SUBST(SLDFLAGS, [$pldflags])
    fi

    AC_SUBST(STIMING, ["spu_timing"])

        CC=${SCC}
        CFLAGS=${SCFLAGS}
        CXX=${SCXX}
        CXXFLAGS=${SCXXFLAGS}
        LDFLAGS=${SLDFLAGS}
        CCAS=${SCCAS}
        CCASFLAGS=${SCCASFLAGS}

    ])

AC_DEFUN([CCS_SPU_CC_COMPILER_GNU], [

    if test -z "$SCC" ; then
        AC_SUBST(SCC, ["spu-gcc"])
    fi

    if test -z "$SCFLAGS" ; then
        AC_SUBST(SCFLAGS, ["-Wno-main"])
    fi

    if test -z "$SCCAS" ; then
        AC_SUBST(SCCAS, ["spu-gcc"])
    fi

    ])

AC_DEFUN([CCS_SPU_CXX_COMPILER_GNU], [

    if test -z "$SCXX" ; then
        AC_SUBST(SCXX, ["spu-g++"])
    fi

    if test -z "$SCXXFLAGS" ; then
        AC_SUBST(SCXXFLAGS, ["-fno-exceptions -fno-rtti -include spu_intrinsics.h"])
    fi

    if test -z "$SCCAS" ; then
        AC_SUBST(SCCAS, ["spu-gcc"])
    fi

    ])

AC_DEFUN([CCS_SPU_CC_COMPILER_XLC], [

    if test -z "$SCC" ; then
        AC_SUBST(SCC, ["spuxlc"])
    fi

    if test -z "$SCFLAGS" ; then
        AC_SUBST(SCFLAGS, ["-D __ALTIVEC_LITERAL_STYLE__ -qcpluscmt"])
    fi

    if test -z "$SCCAS" ; then
        AC_SUBST(SCCAS, ["spuxlc"])
    fi

    ])
