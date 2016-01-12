AC_DEFUN([CCS_WITH_DACS_LAUNCH], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([DACS_LAUNCH],
        [Specify DaCS launch style <local,remote>])

    AC_ARG_WITH([dacs-launch],
        [AC_HELP_STRING([--with-dacs-launch],
        [Specify DaCS launch style <local,remote>])],
        [
            AC_MSG_CHECKING(DaCS launch style)
            if test -n "$DACS_LAUNCH" ; then
		dacs_launch="$DACS_LAUNCH"
            elif test "$withval" != no ; then
		dacs_launch="$withval"
            else		
		dacs_launch="local"
            fi
        ],
        [
            AC_MSG_CHECKING(DaCS launch style)
            if test -n "$DACS_LAUNCH" ; then
		dacs_launch="$DACS_LAUNCH"
            else
		dacs_launch="local"
            fi
        ])

	AC_SUBST(DACS_LAUNCH,$dacs_launch)
	AC_MSG_RESULT($dacs_launch)
    ])
