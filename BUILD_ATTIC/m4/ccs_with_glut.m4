AC_DEFUN([CCS_WITH_GLUT], [

    AC_ARG_VAR([GLUT],
        [User specified path to GLUT library top-level directory])
    AC_ARG_VAR([GLUT_CPPFLAGS], [User specified path to GLUT library includes])
    AC_ARG_VAR([GLUT_LDFLAGS], [User specified path to GLUT libraries])
    AC_ARG_VAR([GLUT_LIBS],
        [User specified link library names, e.g., -lglut])

    AC_ARG_WITH([glut],
        [AC_HELP_STRING([--with-glut],
        [User specified path to SPE libraries])],
        [
            if test -n "$GLUT" ; then
                if test -z "$GLUT_CPPFLAGS" ; then
                    glut_CPPFLAGS="-I$GLUT/include"
                else
                    glut_CPPFLAGS="$GLUT_CPPFLAGS"
                fi

                if test -z "$GLUT_LDFLAGS" ; then
                    glut_LDFLAGS="-L$GLUT/lib"
                else
                    glut_LDFLAGS="$GLUT_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$GLUT_CPPFLAGS" ; then
                    glut_CPPFLAGS="-I$withval/include"
                else
                    glut_CPPFLAGS="$GLUT_CPPFLAGS"
                fi

                if test -z "$GLUT_LDFLAGS" ; then
                    glut_LDFLAGS="-L$withval/lib"
                else
                    glut_LDFLAGS="$GLUT_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$GLUT" ; then
                if test -z "$GLUT_CPPFLAGS" ; then
                    glut_CPPFLAGS="-I$GLUT/include"
                else
                    glut_CPPFLAGS="$GLUT_CPPFLAGS"
                fi

                if test -z "$GLUT_LDFLAGS" ; then
                    glut_LDFLAGS="-L$GLUT/lib"
                else
                    glut_LDFLAGS="$GLUT_LDFLAGS"
                fi
            else
                if test -z "$GLUT_CPPFLAGS" ; then
                    glut_CPPFLAGS="-I/usr/include"
                else
                    glut_CPPFLAGS="$GLUT_CPPFLAGS"
                fi

                if test -z "$GLUT_LDFLAGS" ; then
                    glut_LDFLAGS="-L/usr/lib"
                else
                    glut_LDFLAGS="$GLUT_LDFLAGS"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $glut_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $glut_LDFLAGS"

        AC_CHECK_HEADER(glut.h, [glut_h=yes], [glut_h=no], /* check */)

        if test $glut_h = no ; then
            AC_CHECK_HEADER(GL/glut.h, [glut_h=yes], [glut_h=no], /* check */)
        fi

        dnl Kluge for Apple frameworks
        if test "$glut_h" != yes ; then
		  		glut_CPPFLAGS="-framework GLUT -framework OpenGL"
        		CPPFLAGS="$old_CPPFLAGS $glut_CPPFLAGS"
            AC_CHECK_HEADER(GLUT/glut.h, [glut_h=yes], [glut_h=no], /* check */)
            if test $glut_h = no; then
              GLUT_CPPFLAGS=""
            fi
        fi

        AC_CHECK_LIB(glut, glFinish, [glut_a=yes], [glut_a=no], [])

        dnl Kluge for Apple frameworks
        if test "$glut_a" != yes ; then
        		glut_LDFLAGS=""
        		LDFLAGS="$old_LDFLAGS"
        		GLUT_LIBS="-framework GLUT -framework OpenGL"
            AC_CHECK_LIB(c, glFinish, [glut_a=yes], [glut_a=no], [$GLUT_LIBS])
            if test $glut_a = no; then
              GLUT_LIBS=""
            fi
        fi

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$glut_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid glut.h)
        else
            AC_SUBST(GLUT_CPPFLAGS, [$glut_CPPFLAGS])
        fi

        if test "$glut_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libglut.a)
        else
            AC_SUBST(GLUT_LDFLAGS, [$glut_LDFLAGS])
        fi

        if test -z "$GLUT_LIBS" ; then
            AC_SUBST(GLUT_LIBS, ["-lglut"])
        fi
    ])
