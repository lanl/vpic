AC_DEFUN([CCS_WITH_OPENSSL], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([OPENSSL],
        [User specified path to SPE library top-level directory])
    AC_ARG_VAR([OPENSSL_CPPFLAGS],
		[User specified path to SPE library includes])
    AC_ARG_VAR([OPENSSL_LDFLAGS], [User specified path to SPE libraries])
    AC_ARG_VAR([OPENSSL_LIBS],
        [User specified link library names, e.g., -lspe -lpthreads])

    AC_ARG_WITH([openssl],
        [AC_HELP_STRING([--with-openssl],
        [User specified path to OpenSSL library])],
        [
            if test -n "$OPENSSL" ; then
                if test -z "$OPENSSL_CPPFLAGS" ; then
                    openssl_CPPFLAGS="-I$OPENSSL/include"
                else
                    openssl_CPPFLAGS="$OPENSSL_CPPFLAGS"
                fi

                if test -z "$OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        openssl_LDFLAGS="-L$OPENSSL/lib"
                    else
                        openssl_LDFLAGS="-L$OPENSSL/lib64"
                    fi
                else
                    openssl_LDFLAGS="$OPENSSL_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$OPENSSL_CPPFLAGS" ; then
                    openssl_CPPFLAGS="-I$withval/include"
                else
                    openssl_CPPFLAGS="$OPENSSL_CPPFLAGS"
                fi

                if test -z "$OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        openssl_LDFLAGS="-L$withval/lib"
                    else
                        openssl_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    openssl_LDFLAGS="$OPENSSL_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$OPENSSL" ; then
                if test -z "$OPENSSL_CPPFLAGS" ; then
                    openssl_CPPFLAGS="-I$OPENSSL/include"
                else
                    openssl_CPPFLAGS="$OPENSSL_CPPFLAGS"
                fi

                if test -z "$OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        openssl_LDFLAGS="-L$OPENSSL/lib"
                    else
                        openssl_LDFLAGS="-L$OPENSSL/lib64"
                    fi
                else
                    openssl_LDFLAGS="$OPENSSL_LDFLAGS"
                fi
            else
                openssl_CPPFLAGS="-I/usr/include"
                if test "$BITS" = "32" ; then
                    openssl_LDFLAGS="-L/usr/lib"
                else
                    openssl_LDFLAGS="-L/usr/lib64"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $openssl_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $openssl_LDFLAGS"

        AC_CHECK_HEADER(openssl/evp.h,
			[openssl_h=yes], [openssl_h=no], /* check */)
        AC_CHECK_LIB(ssl, EVP_MD_CTX_init, [openssl_a=yes],
            [openssl_a=no], -lcrypto)

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$openssl_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid openssl/evp.h)
        else
            AC_SUBST(OPENSSL_CPPFLAGS, [$openssl_CPPFLAGS])
        fi

        if test "$openssl_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libssl.a)
        else
            AC_SUBST(OPENSSL_LDFLAGS, [$openssl_LDFLAGS])
        fi

		if test -z "$OPENSSL_LIBS" ; then
			AC_SUBST(OPENSSL_LIBS, ["-lssl -lcrypto"])
		fi
    ])
