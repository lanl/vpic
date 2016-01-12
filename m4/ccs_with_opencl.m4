AC_DEFUN([CCS_WITH_OPENCL], [

    AC_ARG_VAR([OPENCL],
        [User specified path to OPENCL library top-level directory])
    AC_ARG_VAR([OPENCL_CPPFLAGS],
        [User specified path to OPENCL library includes])
    AC_ARG_VAR([OPENCL_LDFLAGS], [User specified path to OPENCL libraries])
    AC_ARG_VAR([OPENCL_LIBS],
        [User specified link library names, e.g., -lopencl])

    AC_ARG_WITH([opencl],
        [AC_HELP_STRING([--with-opencl],
        [User specified path to SPE libraries])],
        [
            if test -n "$OPENCL" ; then
                if test -z "$OPENCL_CPPFLAGS" ; then
                    opencl_CPPFLAGS="-I$OPENCL/include"
                else
                    opencl_CPPFLAGS="$OPENCL_CPPFLAGS"
                fi

                if test -z "$OPENCL_LDFLAGS" ; then
                    opencl_LDFLAGS="-L$OPENCL/lib"
                else
                    opencl_LDFLAGS="$OPENCL_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$OPENCL_CPPFLAGS" ; then
                    opencl_CPPFLAGS="-I$withval/include"
                else
                    opencl_CPPFLAGS="$OPENCL_CPPFLAGS"
                fi

                if test -z "$OPENCL_LDFLAGS" ; then
                    opencl_LDFLAGS="-L$withval/lib"
                else
                    opencl_LDFLAGS="$OPENCL_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$OPENCL" ; then
                if test -z "$OPENCL_CPPFLAGS" ; then
                    opencl_CPPFLAGS="-I$OPENCL/include"
                else
                    opencl_CPPFLAGS="$OPENCL_CPPFLAGS"
                fi

                if test -z "$OPENCL_LDFLAGS" ; then
                    opencl_LDFLAGS="-L$OPENCL/lib"
                else
                    opencl_LDFLAGS="$OPENCL_LDFLAGS"
                fi
            else
                if test -z "$OPENCL_CPPFLAGS" ; then
                    opencl_CPPFLAGS="-I/usr/include"
                else
                    opencl_CPPFLAGS="$OPENCL_CPPFLAGS"
                fi

                if test -z "$OPENCL_LDFLAGS" ; then
                    opencl_LDFLAGS="-L/usr/lib"
                else
                    opencl_LDFLAGS="$OPENCL_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $opencl_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $opencl_LDFLAGS"

        AC_CHECK_HEADER(opencl.h, [opencl_h=yes], [opencl_h=no], /* check */)

        dnl Kluge for Apple frameworks
        if test "$opencl_h" != yes ; then
		  		opencl_CPPFLAGS="-framework OpenCL"
				CPPFLAGS="$old_CPPFLAGS $opencl_CPPFLAGS"
            AC_CHECK_HEADER(OpenCL/opencl.h, [opencl_h=yes],
					[opencl_h=no], /* check */)
            if test $opencl_h = no; then
              opencl_CPPFLAGS=""
            fi 
        fi

        AC_CHECK_LIB(opencl, clCreateBuffer, [opencl_a=yes], [opencl_a=no], [])

        dnl Kluge for Apple frameworks
        if test "$opencl_a" != yes ; then
		  		opencl_LDFLAGS=""
				LDFLAGS="$old_LDFLAGS"
				OPENCL_LIBS="-framework OpenCL"
            AC_CHECK_LIB(c, clCreateBuffer, [opencl_a=yes], [opencl_a=no], [])
            if test $opencl_a = no; then
              OPENCL_LIBS=""
            fi 
        fi

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$opencl_h" = yes && test "$opencl_a" = yes; then
          if test -z "$OPENCL_LIBS" ; then
            OPENCL_LIBS=-lopencl
          fi
          AC_SUBST(OPENCL_CPPFLAGS, [$opencl_CPPFLAGS])
          AC_SUBST(OPENCL_LDFLAGS, [$opencl_LDFLAGS])
          AC_SUBST(OPENCL_LIBS, [$OPENCL_LIBS])
          AC_DEFINE(HAVE_OPENCL,1)
        else
          AC_DEFINE(HAVE_OPENCL,0)
        fi
    ])
