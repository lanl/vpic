AC_DEFUN([CCS_WITH_BUILDSTYLE], [
    AC_ARG_VAR([BUILDSTYLE], [Specify build style])
    AC_ARG_WITH([buildstyle],
        [AC_HELP_STRING([--with-buildstyle],
        [Specify build style])],
        [
            AC_MSG_CHECKING(Build Style)
            if test -n "$BUILDSTYLE" ; then
                CCS_SET_BUILDSTYLE($BUILDSTYLE)
            elif test "$withval" != no ; then
                CCS_SET_BUILDSTYLE($withval)
            else
                CCS_SET_BUILDSTYLE(standard)
            fi
        ],
        [
            if test -n "$BUILDSTYLE" ; then
                CCS_SET_BUILDSTYLE($BUILDSTYLE)
            else
                CCS_SET_BUILDSTYLE(standard)
            fi
    ])
])

AC_DEFUN([CCS_SET_BUILDSTYLE], [
    AC_MSG_CHECKING(Build Style)
    case "$1" in

        standard)
            AC_DEFINE([BUILDSTYLE], [standard],
                [Standard build])
            AC_SUBST(buildstyle, [standard])
            AC_MSG_RESULT([standard])
        ;;

        standard_cell)
            AC_DEFINE([BUILDSTYLE], [standard_cell],
                [Standard build - IBM Cell PP/SPEs])
            AC_SUBST(buildstyle, [standard_cell])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([standard_cell])
        ;;

        standard_cell_ppe)
            AC_DEFINE([BUILDSTYLE], [standard_cell_ppe],
                [Standard build - IBM Cell PP/SPEs])
            AC_SUBST(buildstyle, [standard_cell_ppe])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([standard_cell_ppe])
        ;;

        ompi_relay)
            AC_DEFINE([BUILDSTYLE], [ompi_relay],
                [OpenMPI relay build])
            AC_SUBST(buildstyle, [ompi_relay])
            AC_MSG_RESULT([ompi_relay])
        ;;

        ompi_relay_cell)
            AC_DEFINE([BUILDSTYLE], [ompi_relay_cell],
                [OpenMPI relay build - IBM Cell PPE/SPEs])
            AC_SUBST(buildstyle, [ompi_relay_cell])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([ompi_relay_cell])
        ;;

        ompi_relay_cell_ppe)
            AC_DEFINE([BUILDSTYLE], [ompi_relay_cell_ppe],
                [OpenMPI relay build - IBM Cell PPE])
            AC_SUBST(buildstyle, [ompi_relay_cell_ppe])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([ompi_relay_cell_ppe])
        ;;

        ompi_relay_hybrid_cell)
            AC_DEFINE([BUILDSTYLE], [ompi_relay_hybrid_cell],
                [Hybrid OpenMPI relay build with IBM Cell PPE/SPEs])
            AC_SUBST(buildstyle, [ompi_relay_hybrid_cell])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([ompi_relay_hybrid_cell])
        ;;

        ompi_relay_hybrid_cell_ppe)
            AC_DEFINE([BUILDSTYLE], [ompi_relay_hybrid_cell_ppe],
                [Hybrid OpenMPI relay build - Cell PPE only])
            AC_SUBST(buildstyle, [ompi_relay_hybrid_cell_ppe])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([ompi_relay_hybrid_cell_ppe])
        ;;

        dacs_relay_hybrid_cell)
            AC_DEFINE([BUILDSTYLE], [dacs_relay_hybrid_cell],
                [Hybrid DaCS relay build with IBM Cell PPE/SPEs])
            AC_SUBST(buildstyle, [dacs_relay_hybrid_cell])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([dacs_relay_hybrid_cell])
        ;;

        dacs_relay_hybrid_cell_ppe)
            AC_DEFINE([BUILDSTYLE], [dacs_relay_hybrid_cell_ppe],
                [Hybrid DaCS relay build - Cell PPE only])
            AC_SUBST(buildstyle, [dacs_relay_hybrid_cell_ppe])
            AC_SUBST(enable_cell, [1])
            AC_MSG_RESULT([dacs_relay_hybrid_cell_ppe])
        ;;

        *)
            AC_MSG_ERROR(Invalid Build Style Option)
        ;;

    esac
])
