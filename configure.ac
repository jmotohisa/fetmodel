#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.56)
AC_INIT(fetmodel, 0.2a, motohisa@ist.hokudai.ac.jp)
AC_CONFIG_SRCDIR([src/ccm.c])

SHARED_VERSION_INFO="0:0:0" # CURRENT:REVISION:AGE

AC_CONFIG_HEADER([src/config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_SUBST(SHARED_VERSION_INFO)

AM_INIT_AUTOMAKE([1.11])
dnl LT_INIT

# need shared libraries for python interfaces
AC_ENABLE_SHARED

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL

##############################################################################
# Compiler flags

AC_ARG_ENABLE(debug,
              [AC_HELP_STRING([--enable-debug],[compile for debugging])],
	      enable_debug=$enableval, enable_debug=no)
if test "$enable_debug" = "yes"; then
	CFLAGS="-g -Wall -O0"
	CXXFLAGS="-g -Wall -O0"
	FFLAGS="-g -Wall -O0"
	AC_DEFINE(DEBUG,1,[define to enable debugging code])
fi

##############################################################################
# More checks
AC_PROG_LIBTOOL


# Checks for libraries.

AC_CHECK_LIB(m,main) 
##############################################################################
# GNU Scientific Library

AC_CHECK_FUNC(cblas_cgemm, [], [AC_CHECK_LIB(gslcblas, cblas_cgemm)])

AC_CHECK_LIB(gsl, gsl_sf_bessel_Jn, [],
	[AC_MSG_ERROR([Missing GNU GSL library])])

##############################################################################
# Check for pkg-config
AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
if test $PKG_CONFIG = no; then
  AC_MSG_ERROR([Cannot find pkg-config, make sure it is installed and in your PATH])
fi

##############################################################################

AC_ARG_WITH(libctl, [AC_HELP_STRING([--with-libctl=<dir>],[specify libctl directory])], with_libctl=$withval, with_libctl=yes)
if test "x$with_libctl" = "xno"; then
        :
elif test "x$with_libctl" != "xyes"; then
	LIBCTL_DIR="$with_libctl"
fi

##############################################################################

# subsequent libraries are only for libctl front end
pcsmatrix_save_LIBS=$LIBS

AC_ARG_WITH(libctl, [AC_HELP_STRING([--with-libctl=<dir>],[specify libctl directory])], with_libctl=$withval, with_libctl=yes)
if test "x$with_libctl" = "xno"; then
        :
elif test "x$with_libctl" != "xyes"; then
	LIBCTL_DIR="$with_libctl"
fi

AM_CONDITIONAL(WITH_LIBCTL, test "x$with_libctl" != xno)


##############################################################################
# Check for Guile library and its behavior (for libctl front end):

if test "x$with_libctl" != xno; then

AC_CHECK_PROG(GUILE_CONFIG, guile-config, guile-config, unknown)
if test "x$GUILE_CONFIG" = "xunknown"; then
        AC_CHECK_LIB(readline, readline)
        AC_CHECK_LIB(dl, dlopen)
        AC_CHECK_LIB(guile-ltdl, scm_lt_dlinit, [], [AC_CHECK_LIB(ltdl, dlopen)])
        AC_CHECK_LIB(guile, gh_eval_str, [],
                     [AC_CHECK_LIB(guile, scm_eval_string, [],
		         [AC_MSG_ERROR([linking to guile failed])])])
else
        CPPFLAGS="$CPPFLAGS `$GUILE_CONFIG compile`"
        LIBS="`$GUILE_CONFIG link` $LIBS"
	AC_MSG_CHECKING([if linking to guile works])
	AC_TRY_LINK_FUNC(gh_enter, AC_MSG_RESULT(yes),
 	  [AC_TRY_LINK_FUNC(scm_boot_guile, AC_MSG_RESULT(yes),
            [AC_MSG_RESULT(no)
             AC_MSG_ERROR([guile-config is broken])])])
fi
AC_CHECK_HEADERS([libguile.h guile/gh.h])

# Check how smob types work in this Guile version:

AC_CHECK_FUNCS(scm_make_smob_type)
AC_MSG_CHECKING([for SCM_SMOB_PREDICATE])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_PREDICATE(1,x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_PREDICATE, 1, [define if we have SCM_SMOB_PREDICATE])
fi
AC_MSG_CHECKING([for SCM_SMOB_DATA])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_SMOB_DATA(x);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_SMOB_DATA, 1, [define if we have SCM_SMOB_DATA])
fi
AC_MSG_CHECKING([for SCM_NEWSMOB])
AC_TRY_LINK([#include <guile/gh.h>
#include <libguile/smob.h>],
                    [int x; SCM_NEWSMOB(x,1,0);],
                    ok=yes, ok=no)
AC_MSG_RESULT($ok)
if test $ok = yes; then
	AC_DEFINE(HAVE_SCM_NEWSMOB, 1, [define if we have SCM_NEWSMOB])
fi

# Check how to activate Guile readline support:

ACTIVATE_READLINE=""
if test "x$cross_compiling" = xno; then
AC_MSG_CHECKING(how to activate readline in Guile)
ractivate="(use-modules (readline-activator)) (activate-readline)"
if guile -c "$ractivate" > /dev/null 2>&1; then
	AC_MSG_RESULT(readline-activator)
	ACTIVATE_READLINE="$ractivate"
else
	ractivate="(use-modules (ice-9 readline)) (activate-readline)"
	if guile -c "$ractivate" >& /dev/null; then
	        AC_MSG_RESULT(ice-9 readline)
	        ACTIVATE_READLINE="$ractivate"
	else
		AC_MSG_RESULT(cannot)
		ACTIVATE_READLINE=""
	fi
fi
fi
AC_SUBST(ACTIVATE_READLINE)

fi # if with_libctl

##############################################################################
# Check for libctl library and files

if test "x$with_libctl" != xno; then

AC_MSG_CHECKING([for libctl dir])
if test x != x"$LIBCTL_DIR" -a ! -r "$LIBCTL_DIR/base/ctl.scm"; then
	LIBCTL_DIR=""
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/local/share/libctl/base/ctl.scm; then
	LIBCTL_DIR="/usr/local/share/libctl"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl3/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl3"
fi
if test x = x"$LIBCTL_DIR" -a -r /usr/share/libctl/base/ctl.scm; then
        LIBCTL_DIR="/usr/share/libctl"
fi
if test x = x"$LIBCTL_DIR"; then
	AC_MSG_ERROR([could not find libctl files; use --with-libctl=<dir>])
fi
AC_MSG_RESULT($LIBCTL_DIR)
AC_SUBST(LIBCTL_DIR)

AC_CHECK_PROGS(GEN_CTL_IO, gen-ctl-io, echo)
if test x"$GEN_CTL_IO" = xecho; then
	AC_MSG_ERROR([could not find gen-ctl-io program; check your PATH])
fi

# check for -lctl:
AC_CHECK_LIB(ctl, ctl_get_vector3, [], [AC_MSG_ERROR([Couldn't find the required libctl library.])])

# check for ctl.h (Fedora puts it in /usr/include/ctl, grrr)
CTL_H_CPPFLAG=""
AC_CHECK_HEADER(ctl.h, [CTL_H_CPPFLAG=" "], [
  save_CPPFLAGS=$CPPFLAGS
  AC_MSG_CHECKING([for alternate ctl.h location])
  for ctl_h_path in /usr/include/ctl /usr/local/include/ctl; do
     CPPFLAGS="$save_CPPFLAGS -I$ctl_h_path"
     AC_COMPILE_IFELSE(AC_LANG_PROGRAM([[#include <ctl.h>]],[]), [CTL_H_CPPFLAG=-I$ctl_h_path; break])
  done
  AC_MSG_RESULT([$CTL_H_CPPFLAG])
  CPPFLAGS=$save_CPPFLAGS])
if test "x$CTL_H_CPPFLAG" = "x"; then
  AC_MSG_ERROR([Couldn't find the ctl.h header file for libctl.])
fi
AC_SUBST(CTL_H_CPPFLAG)

# Check libctl version >= LIBCTL_MAJOR.LIBCTL_MINOR.LIBCTL_BUGFIX
LIBCTL_MAJOR=3; LIBCTL_MINOR=2; LIBCTL_BUGFIX=0
AC_MSG_CHECKING([whether libctl version is at least ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX}])
AC_EGREP_CPP(yes, [[
#include <ctl.h>
#if LIBCTL_MAJOR_VERSION > $LIBCTL_MAJOR || (LIBCTL_MAJOR_VERSION == $LIBCTL_MAJOR && (LIBCTL_MINOR_VERSION > $LIBCTL_MINOR || (LIBCTL_MINOR_VERSION == $LIBCTL_MINOR && LIBCTL_BUGFIX_VERSION >= $LIBCTL_BUGFIX)))
  yes
#endif
]], [AC_MSG_RESULT(ok)], [AC_MSG_ERROR([libctl version ${LIBCTL_MAJOR}.${LIBCTL_MINOR}.${LIBCTL_BUGFIX} or later is required])])

# Check for libctl_quiet feature (in libctl 3.1.1 or later):
AC_MSG_CHECKING([for libctl_quiet feature])
if grep libctl_quiet "$LIBCTL_DIR/base/main.c" 2>&1 > /dev/null; then
   AC_MSG_RESULT(yes)
   AC_DEFINE([HAVE_LIBCTL_QUIET], [1], [If we have the libctl_quiet variable])
else
   AC_MSG_RESULT(no)
fi

# On IRIX, basename/dirname functions in libctl/main.c require -lgen
AC_CHECK_LIB(gen, basename)

fi # if with_libctl

##############################################################################

LIBCTL_LIBS=$LIBS
AC_SUBST(LIBCTL_LIBS)

##############################################################################

# checks for python
##############################################################################
AC_ARG_WITH(python,
	[AC_HELP_STRING([--without-python], [compile without Python interface])],
        with_python=$withval,with_python=yes)

if test "x$with_python" = xno; then
  have_python=no
else
  dnl if test "x$have_libctlgeom" = xno; then
  dnl   AC_MSG_ERROR([libctlgeom was not found, is required for Python interface])
  dnl fi

  if test "$enable_shared" = no; then
     AC_MSG_WARN([Python interface requires --enable-shared; disabling])
     have_python=no
  else

    AM_PATH_PYTHON([],[have_python=yes],[have_python=no])
    if test $have_python = yes; then

      AC_MSG_CHECKING([for Python include directory])
      pinc=`echo "import distutils.sysconfig; print (distutils.sysconfig.get_python_inc())" | $PYTHON - 2>/dev/null`
      AC_MSG_RESULT([${pinc:-unknown}])
      test -n "$pinc" && PYTHON_INCLUDES="-I$pinc"
      save_CPPFLAGS=$CPPFLAGS
      CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
      AC_CHECK_HEADER([Python.h], [], [AC_MSG_WARN([disabling Python wrappers])
                                       have_python=no])

      if test $have_python = yes; then
        AC_MSG_CHECKING([for Numpy include directory])
        pinc=`echo "import numpy; print (numpy.get_include())" | $PYTHON - 2>/dev/null`
        AC_MSG_RESULT([${pinc:-unknown}])
        test -n "$pinc" && PYTHON_INCLUDES="$PYTHON_INCLUDES -I$pinc"
        CPPFLAGS="$save_CPPFLAGS $PYTHON_INCLUDES"
        AC_CHECK_HEADER([numpy/arrayobject.h],[],[
          AC_MSG_WARN([disabling Python wrappers])
          have_python=no],[#include <Python.h>])

        AC_MSG_CHECKING([for coverage module])
        $PYTHON -c 'import coverage' 2>/dev/null
        if test $? = 0; then
          AC_MSG_RESULT([yes])
          have_coverage=yes
        else
          AC_MSG_RESULT([no])
          have_coverage=no
        fi
      fi

      CPPFLAGS=$save_CPPFLAGS
###	jm test
	  dnl PY_PREFIX=`$PYTHON -c 'import sys ; print(sys.prefix)'`
	  dnl PYTHON_LIBS="-lpython$PYTHON_VERSION"
	  dnl PYTHON_CFLAGS="-I$PY_PREFIX/include/python$PYTHON_VERSION"
	dnl pythondir=`echo "$pythondir" | sed "s,$PYTHON_PREFIX,$myown_pfx,"` 
### 
    fi # have_python

  fi # enabled_shared
fi # with_python

AC_SUBST(PYTHON_INCLUDES)
AM_CONDITIONAL(WITH_PYTHON, test x"$have_python" = "xyes")
AM_CONDITIONAL(WITH_COVERAGE, test x"$have_coverage" = "xyes")

### jm test
dnl AC_SUBST([PYTHON_LIBS])
dnl AC_SUBST([PYTHON_CFLAGS])
###
##############################################################################
# Checks for SWIG
AX_PKG_SWIG(1.3.21, [], [AC_MSG_ERROR([SWIG is required to build..]) ])
dnl AX_SWIG_ENABLE_CXX
dnl AX_SWIG_PYTHON
dnl AC_CHECK_PROG(SWIG, swig, swig, echo)
dnl if test "$USE_MAINTAINER_MODE" = yes && (test "x$with_python" = "xyes" || test "x$with
dnl _scheme" = "xyes") && test "x$SWIG" = xecho; then
dnl   AC_MSG_ERROR([SWIG not found; configure --without-python --without-scheme or use a r
dnl elease tarball])
dnl fi

##############################################################################
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

if test "x$with_libctl" = xyes; then
   AC_CONFIG_FILES([libctl/Makefile libctl/fetmodel.scm])
fi

if test "x$with_python" = xyes; then
   AC_CONFIG_FILES([python/Makefile])
fi

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
