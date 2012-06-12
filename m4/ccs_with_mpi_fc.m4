dnl
dnl Discover MPI library attributes for Fortran
dnl
dnl Usage: CCS_WITH_MPI_FC
dnl
AC_DEFUN([CCS_WITH_MPI_FC], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([MPIFCDIR],
        [User specified path to Fortran MPI library top-level directory])
    AC_ARG_VAR([MPIFC_CPPFLAGS], [User specified path to Fortran MPI library includes])
    AC_ARG_VAR([MPIFC_LDFLAGS], [User specified path to Fortran MPI libraries])
    AC_ARG_VAR([MPIFC_LIBS],
        [User specified Fortran link library names, e.g., -lmpi])

    AC_ARG_WITH([mpifc],
        [AC_HELP_STRING([--with-mpifc],
        [User specified path to MPI FC libraries])],
        [
            if test -n "$MPIFCDIR" ; then
                if test -z "$MPIFC_CPPFLAGS" ; then
                    mpifc_CPPFLAGS="-I$MPIFCDIR/include"
                else
                    mpifc_CPPFLAGS="$MPIFC_CPPFLAGS"
                fi

                if test -z "$MPIFC_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpifc_LDFLAGS="-L$MPIFCDIR/lib"
                    else
                        mpifc_LDFLAGS="-L$MPIFCDIR/lib64"
                    fi
                else
                    mpifc_LDFLAGS="$MPIFC_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$MPIFC_CPPFLAGS" ; then
                    mpifc_CPPFLAGS="-I$withval/include"
                else
                    mpifc_CPPFLAGS="$MPIFC_CPPFLAGS"
                fi

                if test -z "$MPIFC_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpifc_LDFLAGS="-L$withval/lib"
                    else
                        mpifc_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    mpifc_LDFLAGS="$MPIFC_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$MPIFCDIR" ; then
                if test -z "$MPIFC_CPPFLAGS" ; then
                    mpifc_CPPFLAGS="-I$MPIFCDIR/include"
                else
                    mpifc_CPPFLAGS="$MPIFC_CPPFLAGS"
                fi

                if test -z "$MPIFC_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpifc_LDFLAGS="-L$MPIFCDIR/lib"
                    else
                        mpifc_LDFLAGS="-L$MPIFCDIR/lib64"
                    fi
                else
                    mpifc_LDFLAGS="$MPIFC_LDFLAGS"
                fi
            else
                AC_PATH_PROG(MPIFCCMP, mpif90, no)

                if test "$MPIFCCMP" = "no" ; then
                    if test -z "$MPIFC_CPPFLAGS" ; then
                        mpifc_CPPFLAGS="-I/usr/include"
                    else
                        mpifc_CPPFLAGS="$MPIFC_CPPFLAGS"
                    fi

                    if test -z "$MPIFC_LDFLAGS" ; then
                        if test "$BITS" = "32" ; then
                            mpifc_LDFLAGS="-L/usr/lib"
                        else
                            mpifc_LDFLAGS="-L/usr/lib64"
                        fi
                    else
                        mpifc_LDFLAGS="$MPIFC_LDFLAGS"
                    fi
                else
		    if test -z "$MPIFC_CPPFLAGS" ; then
                        mpifc_CPPFLAGS=`$MPIFCCMP --showme:compile`
                    else
                        mpifc_CPPFLAGS="$MPIFC_CPPFLAGS"
                    fi

		    if test -z "$MPIFC_LDFLAGS" ; then
                        mpifc_LDFLAGS=`$MPIFCCMP --showme:link | \
                            sed 's,-[[lW]][[^ ]]*,,g'`
                    else
                        mpifc_LDFLAGS="$MPIFC_LDFLAGS"
                    fi

                    mpifc_LIBS=`$MPIFCCMP --showme:link | \
                        sed 's, -[[^lW]][[^ ]]*,,g;s,-pthread,,g' | \
                        sed 's, *,,' | sed 's, $,,'`
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $mpifc_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $mpifc_LDFLAGS"

        if test -z "$MPIFC_LIBS" ; then
            if test -z "$mpifc_LIBS" ; then
                AC_SUBST(MPIFC_LIBS, ["-lmpi -lorte -lopal"])
            else
                AC_SUBST(MPIFC_LIBS, [$mpifc_LIBS])
            fi
        else
            AC_SUBST(MPIFC_LIBS, [$MPIFC_LIBS])
        fi

        AC_LANG_PUSH(Fortran)
        mpifcLIB=mpi_f90
        supportLIBS=`echo $MPIFC_LIBS | sed s,-lmpi_f90,,g`
 
        #AC_CHECK_HEADER(mpif.h, [mpif_h=yes], [mpif_h=no], /* check */)
        #AC_CHECK_FILE($mpifc_CPPFLAGS/mpif.h, [mpif_h=yes], [mpif_h=no])
        #AC_CHECK_LIB($mpifcLIB, mpi_init, [mpi_a=yes],
        #   [mpi_a=no], [$supportLIBS])
        mpif_h=yes
        mpi_a=yes

        AC_LANG_POP()

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$mpif_h" != yes ; then
          AC_MSG_ERROR(Failed to find valid mpif.h)
        else
          AC_SUBST(MPIFC_CPPFLAGS, [$mpifc_CPPFLAGS])
        fi

        if test "$mpi_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmpi.a)
        else
            AC_SUBST(MPIFC_LDFLAGS, [$mpifc_LDFLAGS])
        fi

    ])
