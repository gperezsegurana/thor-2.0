# Process this file with autoconf to produce a configure script.
#AC_INIT
AC_CONFIG_SRCDIR([thor/src/field.cc])
AC_INIT([thor],[1.0])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(src/thor.cc)
AM_INIT_AUTOMAKE()
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
# Use GNU gcc; vs. LLVM compilers.
AC_PROG_CXX
AC_PROG_F77
# C++ links the project programs, so add the Fortran runtime explicitly for
# the static Thor and TPSA libraries.
AC_F77_LIBRARY_LDFLAGS

# Checks for libraries.

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

# Boost is header-only for the facilities used by these programs. It is often
# on the compiler's default include path on Linux. Homebrew can report its
# prefix on macOS; elsewhere pass the package prefix with --with-boost=DIR.
AC_ARG_WITH([boost],
  [AS_HELP_STRING([--with-boost=DIR],
                  [prefix of a Boost install, if its headers are not in the compiler's default path])])
BOOST_CFLAGS=""
AS_IF([test -n "$with_boost"],
      [BOOST_CFLAGS="-I$with_boost/include"],
      [AC_PATH_PROG([BREW], [brew], [no])
       AS_IF([test "x$BREW" != xno],
             [thor_brew_boost=`$BREW --prefix boost 2>/dev/null`
              AS_IF([test -d "$thor_brew_boost/include"],
                    [BOOST_CFLAGS="-I$thor_brew_boost/include"])])])
thor_save_CPPFLAGS=$CPPFLAGS
CPPFLAGS="$CPPFLAGS $BOOST_CFLAGS"
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/any.hpp], [],
  [AC_MSG_ERROR([Boost headers not found. Install Boost or pass --with-boost=DIR.])])
AC_LANG_POP([C++])
CPPFLAGS=$thor_save_CPPFLAGS
AC_SUBST([BOOST_CFLAGS])

AS_IF([test ! -d "$THOR_LIB/thor/lib"],
      [AC_MSG_ERROR([THOR_LIB is unset or wrong: '$THOR_LIB'.
       Set it to the Thor tree and build its libraries first.])])
AS_IF([test ! -d "$NUM_REC/lib"],
      [AC_MSG_ERROR([NUM_REC is unset or wrong: '$NUM_REC'.
       Set it to the Numerical Recipes tree.])])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CONFIG_FILES([Makefile src/Makefile src/Powell/Makefile src/Powell/src/Makefile])
AC_OUTPUT

AC_MSG_NOTICE([

  Thor $PACKAGE_VERSION -- projects build configured

    C compiler:      $CC
    C++ compiler:    $CXX
    FORTRAN 77:      $F77
    GSL:             $GSL_VERSION
      gsl-config:    $GSL_CONFIG
    Boost cflags:    ${BOOST_CFLAGS:-(compiler default path)}
    THOR_LIB:        $THOR_LIB
    NUM_REC:         $NUM_REC
    install prefix:  $prefix

  Next: make
])
