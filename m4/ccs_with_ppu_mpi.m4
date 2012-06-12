AC_DEFUN([CCS_WITH_PPU_MPI], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([PPU_MPIDIR],
        [User specified path to PPU MPI library top-level directory])
    AC_ARG_VAR([PPU_MPI_CPPFLAGS], [User specified path to PPU MPI library includes])
    AC_ARG_VAR([PPU_MPI_LDFLAGS], [User specified path to PPU MPI libraries])
    AC_ARG_VAR([PPU_MPI_LIBS],
        [User specified link library names, e.g., -lmpi])

    AC_ARG_WITH([ppu-mpi],
        [AC_HELP_STRING([--with-ppu-mpi],
        [User specified path to MPI libraries])],
        [
            if test -n "$PPU_MPIDIR" ; then
                if test -z "$PPU_MPI_CPPFLAGS" ; then
                    ppu_mpi_CPPFLAGS="-I$PPU_MPIDIR/include"
                else
                    ppu_mpi_CPPFLAGS="$PPU_MPI_CPPFLAGS"
                fi

                if test -z "$PPU_MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_mpi_LDFLAGS="-L$PPU_MPIDIR/lib"
                    else
                        ppu_mpi_LDFLAGS="-L$PPU_MPIDIR/lib64"
                    fi
                else
                    ppu_mpi_LDFLAGS="$PPU_MPI_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$PPU_MPI_CPPFLAGS" ; then
                    ppu_mpi_CPPFLAGS="-I$withval/include"
                else
                    ppu_mpi_CPPFLAGS="$PPU_MPI_CPPFLAGS"
                fi

                if test -z "$PPU_MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_mpi_LDFLAGS="-L$withval/lib"
                    else
                        ppu_mpi_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    ppu_mpi_LDFLAGS="$PPU_MPI_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$PPU_MPIDIR" ; then
                if test -z "$PPU_MPI_CPPFLAGS" ; then
                    ppu_mpi_CPPFLAGS="-I$PPU_MPIDIR/include"
                else
                    ppu_mpi_CPPFLAGS="$PPU_MPI_CPPFLAGS"
                fi

                if test -z "$PPU_MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_mpi_LDFLAGS="-L$PPU_MPIDIR/lib"
                    else
                        ppu_mpi_LDFLAGS="-L$PPU_MPI/lib64"
                    fi
                else
                    ppu_mpi_LDFLAGS="$PPU_MPI_LDFLAGS"
                fi
            else
                if test -z "$PPU_MPI_CPPFLAGS" ; then
                    ppu_mpi_CPPFLAGS="-I/usr/include"
                else
                    ppu_mpi_CPPFLAGS="$PPU_MPI_CPPFLAGS"
                fi

                if test -z "$PPU_MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_mpi_LDFLAGS="-L/usr/lib"
                    else
                        ppu_mpi_LDFLAGS="-L/usr/lib64"
                    fi
                else
                    ppu_mpi_LDFLAGS="$PPU_MPI_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $ppu_mpi_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $ppu_mpi_LDFLAGS"

        if test -z "$PPU_MPI_LIBS" ; then
            AC_SUBST(PPU_MPI_LIBS, ["-lmpi -lorte -lopal"])
        else
            AC_SUBST(PPU_MPI_LIBS, [$PPU_MPI_LIBS])
        fi

        supportLIBS=`echo $PPU_MPI_LIBS | sed s,-lmpi,,g`

        AC_CHECK_HEADER(mpi.h, [ppu_mpi_h=yes], [ppu_mpi_h=no], /* check */)
        AC_CHECK_LIB(mpi, MPI_Init, [ppu_mpi_a=yes],
            [ppu_mpi_a=no], [$supportLIBS])

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$ppu_mpi_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid mpi.h)
        else
            AC_SUBST(PPU_MPI_CPPFLAGS, [$ppu_mpi_CPPFLAGS])
        fi

        if test "$ppu_mpi_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmpi.a)
        else
            AC_SUBST(PPU_MPI_LDFLAGS, [$ppu_mpi_LDFLAGS])
        fi
    ])
