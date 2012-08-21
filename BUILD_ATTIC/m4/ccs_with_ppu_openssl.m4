AC_DEFUN([CCS_WITH_PPU_OPENSSL], [
    AC_REQUIRE([CCS_WITH_ADDRESSING])

    AC_ARG_VAR([PPU_OPENSSL],
        [User specified path to OpenSSL library top-level directory])
    AC_ARG_VAR([PPU_OPENSSL_CPPFLAGS],
		[User specified path to OpenSSL library includes])
    AC_ARG_VAR([PPU_OPENSSL_LDFLAGS],
		[User specified path to OpenSSL libraries])
    AC_ARG_VAR([PPU_OPENSSL_LIBS],
        [User specified link library names, e.g., -lssl])

    AC_ARG_WITH([ppu-openssl],
        [AC_HELP_STRING([--with-ppu-openssl],
        [User specified path to OpenSSL library])],
        [
            if test -n "$PPU_OPENSSL" ; then
                if test -z "$PPU_OPENSSL_CPPFLAGS" ; then
                    ppu_openssl_CPPFLAGS="-I$PPU_OPENSSL/include"
                else
                    ppu_openssl_CPPFLAGS="$PPU_OPENSSL_CPPFLAGS"
                fi

                if test -z "$PPU_OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_openssl_LDFLAGS="-L$PPU_OPENSSL/lib"
                    else
                        ppu_openssl_LDFLAGS="-L$PPU_OPENSSL/lib64"
                    fi
                else
                    ppu_openssl_LDFLAGS="$PPU_OPENSSL_LDFLAGS"
                fi
            elif test "$withval" != no ; then
                if test -z "$PPU_OPENSSL_CPPFLAGS" ; then
                    ppu_openssl_CPPFLAGS="-I$withval/include"
                else
                    ppu_openssl_CPPFLAGS="$PPU_OPENSSL_CPPFLAGS"
                fi

                if test -z "$PPU_OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_openssl_LDFLAGS="-L$withval/lib"
                    else
                        ppu_openssl_LDFLAGS="-L$withval/lib64"
                    fi
                else
                    ppu_openssl_LDFLAGS="$PPU_OPENSSL_LDFLAGS"
                fi
            fi
        ],
        [
            if test -n "$PPU_OPENSSL" ; then
                if test -z "$PPU_OPENSSL_CPPFLAGS" ; then
                    ppu_openssl_CPPFLAGS="-I$PPU_OPENSSL/include"
                else
                    ppu_openssl_CPPFLAGS="$PPU_OPENSSL_CPPFLAGS"
                fi

                if test -z "$PPU_OPENSSL_LDFLAGS" ; then
                    if test "$BITS" = "32" ; then
                        ppu_openssl_LDFLAGS="-L$PPU_OPENSSL/lib"
                    else
                        ppu_openssl_LDFLAGS="-L$PPU_OPENSSL/lib64"
                    fi
                else
                    ppu_openssl_LDFLAGS="$PPU_OPENSSL_LDFLAGS"
                fi
            else
                ppu_openssl_CPPFLAGS="-I/usr/include"
                if test "$BITS" = "32" ; then
                    ppu_openssl_LDFLAGS="-L/usr/lib"
                else
                    ppu_openssl_LDFLAGS="-L/usr/lib64"
                fi
            fi
        ])

        old_CPPFLAGS=$CPPFLAGS
        CPPFLAGS="$CPPFLAGS $ppu_openssl_CPPFLAGS"
        old_LDFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS $ppu_openssl_LDFLAGS"

        AC_CHECK_HEADER(openssl/evp.h,
			[ppu_openssl_h=yes], [ppu_openssl_h=no], /* check */)
        AC_CHECK_LIB(ssl, EVP_MD_CTX_init, [ppu_openssl_a=yes],
            [ppu_openssl_a=no], -lcrypto)

        CPPFLAGS="$old_CPPFLAGS"
        LDFLAGS="$old_LDFLAGS"

        if test "$ppu_openssl_h" != yes ; then
            AC_MSG_ERROR(Failed to find valid openssl/evp.h)
        else
            AC_SUBST(PPU_OPENSSL_CPPFLAGS, [$ppu_openssl_CPPFLAGS])
        fi

        if test "$ppu_openssl_a" != yes ; then
            AC_MSG_ERROR(Failed to find valid libssl.a)
        else
            AC_SUBST(PPU_OPENSSL_LDFLAGS, [$ppu_openssl_LDFLAGS])
        fi

		if test -z "$PPU_OPENSSL_LIBS" ; then
			AC_SUBST(PPU_OPENSSL_LIBS, ["-lssl -lcrypto"])
		fi
    ])
