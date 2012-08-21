dnl ----------------------------------------------------------------------------
dnl CCS_WITH_ARCH_CPUS
dnl ----------------------------------------------------------------------------
AC_DEFUN([CCS_WITH_ARCH_CPUS], [
    AC_ARG_VAR([CCS_ARCH_CPUS], [Specify number of CPUs])

    AC_ARG_WITH([arch-cpus],
        [AC_HELP_STRING([--with-arch-cpus],
        [Specify number of CPUs])],
        [
            if test -n "$CCS_ARCH_CPUS" ; then
                AC_SUBST(CCS_ARCH_CPUS, $CCS_ARCH_CPUS)
            elif test "$withval" != no ; then
                AC_SUBST(CCS_ARCH_CPUS, $withval)
            fi
        ],
        [
            if test -n "$CCS_ARCH_CPUS" ; then
                AC_SUBST(CCS_ARCH_CPUS, $CCS_ARCH_CPUS)
            else
                if test -f "/proc/cpuinfo" ; then
                    cpus=`cat /proc/cpuinfo | grep "cpu cores" | wc -l`
                    AC_SUBST(CCS_ARCH_CPUS, $cpus)
                else
                    AC_SUBST(CCS_ARCH_CPUS, 1)
                fi
            fi
    ])

    AC_MSG_CHECKING(number of CPUs)
    AC_MSG_RESULT($CCS_ARCH_CPUS)
])
