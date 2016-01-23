dnl
dnl Discover MPI library attributes
dnl
dnl Usage: CCS_WITH_MPI(C|C++)
dnl
AC_DEFUN([CCS_WITH_MPI], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([MPIDIR],
        [User specified path to MPI library top-level directory])
    AC_ARG_VAR([MPI_CPPFLAGS], [User specified path to MPI library includes])
    AC_ARG_VAR([MPI_LDFLAGS], [User specified path to MPI libraries])
    AC_ARG_VAR([MPI_LIBS],
        [User specified link library names, e.g., -lmpi])

    AC_ARG_WITH([mpi],
        [AC_HELP_STRING([--with-mpi],
        [User specified path to MPI libraries])],
        [
            if test -n "$MPIDIR" ; then
                if test -z "$MPI_CPPFLAGS" ; then
                    mpi_CPPFLAGS="-I$MPIDIR/include"
                else
                    mpi_CPPFLAGS="$MPI_CPPFLAGS"
                fi

                if test -z "$MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpi_LDFLAGS="-L$MPIDIR/lib"
                    else
                        mpi_LDFLAGS="-L$MPIDIR/lib64"
                    fi
                else
                    mpi_LDFLAGS="$MPI_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$MPI_CPPFLAGS" ; then
                    mpi_CPPFLAGS="-I$withval/include"
                else
                    mpi_CPPFLAGS="$MPI_CPPFLAGS"
                fi

                if test -z "$MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpi_LDFLAGS="-L$withval/lib"
                    else
                        mpi_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    mpi_LDFLAGS="$MPI_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$MPIDIR" ; then
                if test -z "$MPI_CPPFLAGS" ; then
                    mpi_CPPFLAGS="-I$MPIDIR/include"
                else
                    mpi_CPPFLAGS="$MPI_CPPFLAGS"
                fi

                if test -z "$MPI_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        mpi_LDFLAGS="-L$MPIDIR/lib"
                    else
                        mpi_LDFLAGS="-L$MPIDIR/lib64"
                    fi
                else
                    mpi_LDFLAGS="$MPI_LDFLAGS"
                fi
            else
                case "$1" in
                    C++)
                        AC_PATH_PROG(MPICMP, mpicxx, no)
                    ;;
                    *) dnl default to C
                        AC_PATH_PROG(MPICMP, mpicc, no)
                    ;;
                esac

                if test "$MPICMP" = "no" ; then
                    if test -z "$MPI_CPPFLAGS" ; then
                        mpi_CPPFLAGS="-I/usr/include"
                    else
                        mpi_CPPFLAGS="$MPI_CPPFLAGS"
                    fi

                    if test -z "$MPI_LDFLAGS" ; then
                        if test "$BITS" = "32" ; then
                            mpi_LDFLAGS="-L/usr/lib"
                        else
                            mpi_LDFLAGS="-L/usr/lib64"
                        fi
                    else
                        mpi_LDFLAGS="$MPI_LDFLAGS"
                    fi
                else
                    if test -z "$MPI_CPPFLAGS" ; then
                        mpi_CPPFLAGS=`$MPICMP --showme:compile`
                    else
                        mpi_CPPFLAGS="$MPI_CPPFLAGS"
                    fi

                    if test -z "$MPI_LDFLAGS" ; then
                        mpi_LDFLAGS=`$MPICMP --showme:link | \
                            sed 's,-[[lW]][[^ ]]*,,g'`
                    else
                        mpi_LDFLAGS="$MPI_LDFLAGS"
                    fi

                    mpi_LIBS=`$MPICMP --showme:link | \
                        sed 's, -[[^lW]][[^ ]]*,,g;s,-pthread,,g' | \
                        sed 's, *,,' | sed 's, $,,'`
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $mpi_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $mpi_LDFLAGS"

        if test -z "$MPI_LIBS" ; then
            if test -z "$mpi_LIBS" ; then
                AC_SUBST(MPI_LIBS, ["-lmpi -lorte -lopal"])
            else
                AC_SUBST(MPI_LIBS, [$mpi_LIBS])
            fi
        else
            AC_SUBST(MPI_LIBS, [$MPI_LIBS])
        fi

        case "$1" in
            C++)
                AC_LANG_PUSH(C++)
                mpiLIB=mpi_cxx
                supportLIBS=`echo $MPI_LIBS | sed s,-lmpi_cxx,,g`
            ;;
            *)
                AC_LANG_PUSH(C)
                mpiLIB=mpi
                supportLIBS=`echo $MPI_LIBS | sed s,-lmpi,,g`
            ;;
        esac

        AC_CHECK_HEADER(mpi.h, [mpi_h=yes], [mpi_h=no], /* check */)
        AC_CHECK_LIB($mpiLIB, MPI_Init, [mpi_a=yes],
            [mpi_a=no], [$supportLIBS])

		dnl try mpich	
		if test "$mpi_a" != yes ; then
			case "$1" in
				C++)
					AC_LANG_PUSH(C++)
					mpiLIB=mpichcxx
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpichcxx,,g`
				;;
				*)
					AC_LANG_PUSH(C)
					mpiLIB=mpich
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpich,,g`
				;;
			esac

			AC_CHECK_HEADER(mpi.h, [mpi_h=yes], [mpi_h=no], /* check */)
			AC_CHECK_LIB($mpiLIB, MPI_Init, [mpi_a=yes],
				[mpi_a=no], [$supportLIBS])
		fi

		dnl try mpich-gcc
		if test "$mpi_a" != yes ; then
			case "$1" in
				C++)
					AC_LANG_PUSH(C++)
					mpiLIB=mpichcxx
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpichcxx,,g`
				;;
				*)
					AC_LANG_PUSH(C)
					mpiLIB=mpich-gcc
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpich-gcc,,g`
				;;
			esac

			AC_CHECK_HEADER(mpi.h, [mpi_h=yes], [mpi_h=no], /* check */)
			AC_CHECK_LIB($mpiLIB, MPI_Init, [mpi_a=yes],
				[mpi_a=no], [$supportLIBS])
		fi

		dnl try mpich-xl
		if test "$mpi_a" != yes ; then
			case "$1" in
				C++)
					AC_LANG_PUSH(C++)
					mpiLIB=mpichcxx
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpichcxx,,g`
				;;
				*)
					AC_LANG_PUSH(C)
					mpiLIB=mpich-xl
					supportLIBS=`echo $MPI_LIBS | sed s,-lmpich-xl,,g`
				;;
			esac

			AC_CHECK_HEADER(mpi.h, [mpi_h=yes], [mpi_h=no], /* check */)
			AC_CHECK_LIB($mpiLIB, MPI_Init, [mpi_a=yes],
				[mpi_a=no], [$supportLIBS])
		fi

        AC_LANG_POP()

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$mpi_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid mpi.h)
        else
            AC_SUBST(MPI_CPPFLAGS, [$mpi_CPPFLAGS])
        fi

        if test "$mpi_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libmpi.a)
        else
            AC_SUBST(MPI_LDFLAGS, [$mpi_LDFLAGS])
        fi
    ])
