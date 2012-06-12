AC_DEFUN([CCS_ENABLE_SPU_OVERLAYS], [
    #
    # User hints...
    #
    AC_ARG_VAR([SPU_OVERLAYS], [Enable SPU overlays])
    AC_ARG_VAR([SPU_LD_SCRIPT], [SPU linker script])

    AC_ARG_ENABLE([spu_overlays],
        [AC_HELP_STRING([--enable-spu_overlays],
        [enable SPU overlays])],
        [
            if test -n "$SPU_OVERLAYS" ; then
                enable_spu_overlays=1
            elif test "$enableval" = "yes" ; then
                enable_spu_overlays=1
            else
                enable_spu_overlays=0
            fi

            if test "$enable_spu_overlays" = "1" ; then
                AC_DEFINE(SPU_OVERLAYS, [1],
                    [define is you want IBM CBEA build enabled])

                if test -z "$SPU_LD_SCRIPT" ; then
                    AC_MSG_ERROR(Set SPU_LD_SCRIPT to valid linker script)
                else
                    if test -f "$SPU_LD_SCRIPT" ; then
                        AC_MSG_ERROR($SPU_LD_SCRIPT does not exist)
                    fi
                fi
            fi

            AC_SUBST(ENABLE_SPU_OVERLAYS, [$enable_spu_overlays])
        ],
        [
            if test -n "$SPU_OVERLAYS" ; then
                enable_spu_overlays=1
            else
                enable_spu_overlays=0
            fi

            if test "$enable_spu_overlays" = "1" ; then
                AC_DEFINE(SPU_OVERLAYS, [1],
                    [define is you want IBM CBEA build enabled])

                if test -z "$SPU_LD_SCRIPT" ; then
                    AC_MSG_ERROR(Set SPU_LD_SCRIPT to valid linker script)
                else
                    if test -f "$SPU_LD_SCRIPT" ; then
                        AC_MSG_ERROR($SPU_LD_SCRIPT does not exist)
                    fi
                fi
            fi

            AC_SUBST(ENABLE_SPU_OVERLAYS, [$enable_spu_overlays])
    ])
])
