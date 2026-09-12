# Process this file with autoconf to produce a configure script.
#AC_INIT
AC_CONFIG_SRCDIR([thor/src/field.cc])
AC_INIT([thor],[1.0])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(thor/src/field.cc)
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS([config.h])
# Checks for programs.
AC_PROG_AWK
# Let Autoconf discover the compiler, or honour CC, CXX and F77 from the
# environment.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_CPP
AC_PROG_RANLIB
# Use GNU vs. Clang/LLVM.
AC_PROG_CXX
AC_PROG_F77

# GSL supplies gsl-config, which reports the correct include and library paths
# on Linux, Intel macOS and Apple Silicon. --with-gsl is useful when it is not
# on PATH.
AC_ARG_WITH([gsl],
  [AS_HELP_STRING([--with-gsl=DIR],
                  [prefix of a GSL install, if gsl-config is not on PATH])])
AS_IF([test -n "$with_gsl"],
      [GSL_CONFIG="$with_gsl/bin/gsl-config"
       AS_IF([test ! -x "$GSL_CONFIG"],
             [AC_MSG_ERROR([no gsl-config at $GSL_CONFIG; check --with-gsl])])],
      [AC_PATH_PROG([GSL_CONFIG], [gsl-config], [no])])
AS_IF([test "x$GSL_CONFIG" = xno],
      [AC_MSG_ERROR([GSL not found. Install it or pass --with-gsl=DIR.])])
GSL_CFLAGS=`$GSL_CONFIG --cflags`
GSL_LIBS=`$GSL_CONFIG --libs`
GSL_VERSION=`$GSL_CONFIG --version`
AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile
          TPSA/Makefile
          thor/Makefile thor/src/Makefile])
AC_OUTPUT

AC_MSG_NOTICE([

  Thor $PACKAGE_VERSION -- library build configured

    C compiler:      $CC
    C++ compiler:    $CXX
    FORTRAN 77:      $F77
    GSL:             $GSL_VERSION
      gsl-config:    $GSL_CONFIG
    install prefix:  $prefix

  Next: make install
])
