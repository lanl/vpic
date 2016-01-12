AC_DEFUN([CCS_WITH_PPU_COMPILER], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])
    AC_REQUIRE([CCS_ENABLE_STRIPPING])

    AC_ARG_VAR([PPU_COMPILER],
        [User specified PPU C/C++ compiler suite (gnu|xlc)])
    AC_ARG_VAR([PCC], [User specified PPU C compiler command])
    AC_ARG_VAR([PCFLAGS], [User specified PPU C compiler flags])
    AC_ARG_VAR([PCXX], [User specified PPU C++ compiler command])
    AC_ARG_VAR([PCXXFLAGS], [User specified PPU C++ compiler flags])
    AC_ARG_VAR([PLDFLAGS], [User specified PPU linker flags])
    AC_ARG_VAR([PEMB], [User specified PPU embedder])
    AC_ARG_VAR([PEMBFLAGS], [User specified PPU embedder flags])
    AC_ARG_VAR([PAR], [User specified PPU archiver])
    AC_ARG_VAR([PRANLIB], [User specified PPU ranlib])

    AC_ARG_WITH([ppu-compiler],
        [AC_HELP_STRING([--with-ppu-compiler],
        [User specified PPU C/C++ compiler suite (gnu|xlc)])],
        [
            if test -n "$PPU_COMPILER" ; then
                withval=$PPU_COMPILER
            fi

            if test "$withval" = "gnu" ; then
                if test -n "$PCC" ; then
                    if test "$PCC" = "ppuxlc" ; then
                        CCS_PPU_CC_COMPILER_XLC
                    else
                        CCS_PPU_CC_COMPILER_GNU
                    fi
                else
                    CCS_PPU_CC_COMPILER_GNU
                fi

                if test -n "$PCXX" ; then
                    if test "$PCXX" = "ppuxlc++" ; then
                        CCS_PPU_CXX_COMPILER_XLC
                    else
                        CCS_PPU_CXX_COMPILER_GNU
                    fi
                else
                    CCS_PPU_CXX_COMPILER_GNU
                fi
            elif test "$withval" = "xlc" ; then
                if test -n "$PCC" ; then
                    if test "$PCC" = "ppu-gcc" ; then
                        CCS_PPU_CC_COMPILER_GNU
                    else
                        CCS_PPU_CC_COMPILER_XLC
                    fi
                else
                    CCS_PPU_CC_COMPILER_XLC
                fi

                if test -n "$PCXX" ; then
                    if test "$PCXX" = "ppu-g++" ; then
                        CCS_PPU_CXX_COMPILER_GNU
                    else
                        CCS_PPU_CXX_COMPILER_XLC
                    fi
                else
                    CCS_PPU_CXX_COMPILER_XLC
                fi
            fi

            if test -z "$PEMB" ; then
                if test "$BITS" = "32" ; then
                    AC_SUBST(PEMB, ["ppu32-embedspu"])
                else 
                    AC_SUBST(PEMB, ["ppu-embedspu"])
                fi
            fi

            if test -z "$PEMBFLAGS" ; then
                if test "$BITS" = "32" ; then
                    AC_SUBST(PEMBFLAGS, ["-m32"])
                else 
                    AC_SUBST(PEMBFLAGS, ["-m64"])
                fi
            fi

            if test -z "$PAR" ; then
                AC_SUBST(PAR, ["ppu-ar"])
            fi

            if test -z "$PRANLIB" ; then
                AC_SUBST(PRANLIB, ["ppu-ranlib"])
            fi

        ],
        [

            if test -n "$PCC" ; then
                if test "$PCC" = "ppuxlc" ; then
                    echo "using ppuclx"
                    CCS_PPU_CC_COMPILER_XLC
                else
                    CCS_PPU_CC_COMPILER_GNU
                fi
            else
                CCS_PPU_CC_COMPILER_GNU
            fi

            if test -n "$PCXX" ; then
                if test "$PCXX" = "ppuxlc++" ; then
                    CCS_PPU_CXX_COMPILER_XLC
                else
                    CCS_PPU_CXX_COMPILER_GNU
                fi
            else
                CCS_PPU_CXX_COMPILER_GNU
            fi

            if test -z "$PEMB" ; then
                if test "$BITS" = "32" ; then
                    AC_SUBST(PEMB, ["ppu32-embedspu"])
                else 
                    AC_SUBST(PEMB, ["ppu-embedspu"])
                fi
            fi

            if test -z "$PEMBFLAGS" ; then
                if test "$BITS" = "32" ; then
                    AC_SUBST(PEMBFLAGS, ["-m32"])
                else 
                    AC_SUBST(PEMBFLAGS, ["-m64"])
                fi
            fi

            if test -z "$PAR" ; then
                AC_SUBST(PAR, ["ppu-ar"])
            fi

            if test -z "$PRANLIB" ; then
                AC_SUBST(PRANLIB, ["ppu-ranlib"])
            fi

        ])

    if test -n "$PCC" ; then
        CC=$PCC
    fi

    if test "${PCFLAGS=unset}" != "unset" ; then
        CFLAGS=$PCFLAGS
    fi

    if test "${PCXX=unset}" != "unset" ; then
        CXX=$PCXX
    fi

    if test "${PCXXFLAGS=unset}" != "unset" ; then
        CXXFLAGS=$PCXXFLAGS
    fi

    if test -z "$PLDFLAGS" ; then
        pldflags=""
        if test "$BITS" = "32" ; then
            pldflags="-Wl,-m,elf32ppc"
        fi

        if test "$ENABLE_STRIPPING" = "1" ; then
            pldflags="$pldflags -s"
        fi

        AC_SUBST(PLDFLAGS, [$pldflags])
    fi

    if test "{$PLDFLAGS=unset}" != "unset" ; then
        LDFLAGS=$PLDFLAGS
    fi

    if test -n "$PAR" ; then
        AR=$PAR
    fi

    if test -n "$PRANLIB" ; then
        RANLIB=$PRANLIB
    fi

    ])

AC_DEFUN([CCS_PPU_CC_COMPILER_GNU], [
    
    if test -z "$PCC" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCC, ["ppu32-gcc"])
        else
            AC_SUBST(PCC, ["ppu-gcc"])
        fi
    fi

    if test -z "$PCFLAGS" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCFLAGS, ["-maltivec -include altivec.h"])
        else
            AC_SUBST(PCFLAGS, ["-maltivec -include altivec.h"])
        fi
    fi

    if test "$BITS" = "32" ; then
        AC_SUBST(PCFLAGS, ["$PCFLAGS -m32"])
    else
        AC_SUBST(PCFLAGS, ["$PCFLAGS -m64"])
    fi

    ])

AC_DEFUN([CCS_PPU_CXX_COMPILER_GNU], [
    
    if test -z "$PCXX" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCXX, ["ppu32-g++"])
        else
            AC_SUBST(PCXX, ["ppu-g++"])
        fi
    fi

    if test -z "$PCXXFLAGS" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCXXFLAGS, ["-mabi=altivec -maltivec -include altivec.h"])
        else
            AC_SUBST(PCXXFLAGS, ["-mabi=altivec -maltivec -include altivec.h"])
        fi
    fi

    if test "$BITS" = "32" ; then
        AC_SUBST(PCXXFLAGS, ["$PCXXFLAGS -m32"])
    else
        AC_SUBST(PCXXFLAGS, ["$PCXXFLAGS -m64"])
    fi

    ])

AC_DEFUN([CCS_PPU_CC_COMPILER_XLC], [

    if test -z "$PCC" ; then
        AC_SUBST(PCC, ["ppuxlc"])
    fi

    if test -z "$PCFLAGS" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCFLAGS, ["-D __ALTIVEC_LITERAL_STYLE__ -qcpluscmt"])
        else
            AC_SUBST(PCFLAGS, ["-D __ALTIVEC_LITERAL_STYLE__ -qcpluscmt"])
        fi
    fi

    if test "$BITS" = "32" ; then
        AC_SUBST(PCFLAGS, ["$PCFLAGS -q32"])
    else
        AC_SUBST(PCFLAGS, ["$PCFLAGS -q64"])
    fi

    ])

AC_DEFUN([CCS_PPU_CXX_COMPILER_XLC], [

    if test -z "$PCXX" ; then
        AC_SUBST(PCXX, ["ppuxlc++"])
    fi

    if test -z "$PCXXFLAGS" ; then
        if test "$BITS" = "32" ; then
            AC_SUBST(PCXXFLAGS, ["-D __ALTIVEC_LITERAL_STYLE__"])
        else
            AC_SUBST(PCXXFLAGS, ["-D __ALTIVEC_LITERAL_STYLE__"])
        fi
    fi

    if test "$BITS" = "32" ; then
        AC_SUBST(PCXXFLAGS, ["$PCXXFLAGS -q32"])
    else
        AC_SUBST(PCXXFLAGS, ["$PCXXFLAGS -q64"])
    fi

    ])
