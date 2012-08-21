AC_DEFUN([CCS_ENABLE_OPENSSL], [
    #
    # User hints...
    #
    AC_ARG_VAR([ENABLE_OPENSSL], [Enable OpenSSL linkage])
    AC_ARG_ENABLE([openssl],
        [AC_HELP_STRING([--enable-openssl],
        [enable OpenSSL linkage])],
        [
            if test -n "$ENABLE_OPENSSL" ; then
                enable_openssl=1
            elif test "$enableval" = "yes" ; then
                enable_openssl=1
            else
                enable_openssl=0
            fi
        ],
        [
            if test -n "$ENABLE_OPENSSL" ; then
                enable_openssl=1
            else
                enable_openssl=0
            fi
        ])

        AC_SUBST(ENABLE_OPENSSL, [$enable_openssl])

		if test "$enable_openssl" = "1" ; then
			AC_DEFINE(ENABLE_OPENSSL, [1],
				[define if you want to use openssl])
		fi
])
