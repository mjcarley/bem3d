AC_INIT([bem3d], 0.2)
AC_CONFIG_SRCDIR([src/bem3d.h])
AC_CONFIG_MACRO_DIR([m4])

# Making releases:
#   BEM3D_MICRO_VERSION += 1;
#   BEM3D_INTERFACE_AGE += 1;
#   BEM3D_BINARY_AGE += 1;
# if any functions have been added, set BEM3D_INTERFACE_AGE to 0.
# if backwards compatibility has been broken,
# set BEM3D_BINARY_AGE and BEM3D_INTERFACE_AGE to 0.
#
BEM3D_MAJOR_VERSION=0
BEM3D_MINOR_VERSION=2
BEM3D_MICRO_VERSION=1
BEM3D_INTERFACE_AGE=1
BEM3D_BINARY_AGE=1
BEM3D_VERSION=$BEM3D_MAJOR_VERSION.$BEM3D_MINOR_VERSION.$BEM3D_MICRO_VERSION

#AC_SUBST(BEM3D_MAJOR_VERSION)
#AC_SUBST(BEM3D_MINOR_VERSION)
#AC_SUBST(BEM3D_MICRO_VERSION)
AC_SUBST(BEM3D_VERSION)

# libtool versioning
LT_RELEASE=$BEM3D_MAJOR_VERSION.$BEM3D_MINOR_VERSION
LT_CURRENT=`expr $BEM3D_MICRO_VERSION - $BEM3D_INTERFACE_AGE`
LT_REVISION=$BEM3D_INTERFACE_AGE
LT_AGE=`expr $BEM3D_BINARY_AGE - $BEM3D_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$BEM3D_VERSION
PACKAGE=bem3d

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(BEM3D_INTERFACE_AGE, $BEM3D_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(BEM3D_BINARY_AGE, $BEM3D_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_F77

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

PKG_PROG_PKG_CONFIG()

PKG_CHECK_MODULES(gts, gts, [gts_cflags=`pkg-config --cflags-only-I gts`
		        gts_libs=`pkg-config --libs gts`],			
			[AC_CHECK_PROG(GTS, gts-config, yes, no)
			if test x$GTS = xyes ; then
   gts_cflags=`gts-config --cflags`
   gts_libs=`gts-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# GNU Triangulated Surface Library not available or not
  echo \# properly installed. You can get the most recent version from 
  echo \#      http://gts.sourceforge.net/
  echo \#
  echo \#################################################################
  exit
fi
])

AC_CHECK_PROG(GSL, gsl-config, yes, no)
if test x$GSL = xyes ; then
   gsl_cflags=`gsl-config --cflags`
   gsl_libs=`gsl-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# GNU Scientific Library not available or not properly installed.
  echo \# You can get the most recent version from 
  echo \#      http://www.gnu.org/software/gsl/
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(GQR, gqr-config, yes, no)
if test x$GQR = xyes ; then
   AC_DEFINE(HAVE_GQR, 1, [Gauss quadrature library available])
   gqr_cflags=`gqr-config --cflags`
   gqr_libs=`gqr-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# Gauss Quadrature Rule library not installed or not properly
  echo \# installed. You can get it from
  echo \# http://www.paraffinalia.co.uk/Software/downloads.shtml
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_PROG(SISL, sisl-config, yes, no)
if test x$SISL = xyes ; then
   AC_DEFINE(HAVE_SISL, 1, [Simple iterative solver library available])
   sisl_cflags=`sisl-config --cflags`
   sisl_libs=`sisl-config --libs`
else
  AC_DEFINE(HAVE_SISL, 0, [Simple iterative solver library available])
  echo \#################################################################
  echo \#
  echo \# SISL not found: some test programs will not be compiled.
  echo \#
  echo \#################################################################
fi

AC_CHECK_PROG(WMPI, wmpi-config, yes, no)
if test x$WMPI = xyes ; then
   AC_DEFINE(HAVE_WMPI, 1, [MPI wrapper library available])
   wmpi_cflags=`wmpi-config --cflags`
   wmpi_libs=`wmpi-config --libs`
   wmpi-config --cc --mpi >& /dev/null
   status=$?
   if test $status -eq 0  ; then
  echo \#################################################################
  echo \#
  echo \# WMPI parallel support: parallel programs will be compiled
  echo \#
  echo \#################################################################
	mpi_tools=mpitools
   else
  echo \#################################################################
  echo \#
  echo \# WMPI has no parallel support: only serial programs will be compiled
  echo \#
  echo \#################################################################
   fi
else
  AC_DEFINE(HAVE_WMPI, 0, [MPI wrapper library available])
  echo \#################################################################
  echo \#
  echo \# WMPI not found: some test programs will not be compiled.
  echo \#
  echo \#################################################################
fi

AC_CHECK_FILES(fmmlib3d-1.2)

AM_CONDITIONAL([USE_FMMLIB3D_1_2], [test x"$ac_cv_file_fmmlib3d_1_2" = xyes])
##AM_CONDITIONAL([USE_FMMLIB3D_1_2], [test -d fmmlib3d-1.2])

PKG_CHECK_MODULES(libmatheval,libmatheval,
	[AC_DEFINE(HAVE_LIBMATHEVAL, 1, [libmatheval available])
#	matheval_defines="#define BEM3D_HAVE_MATHEVAL"],
	echo libmatheval not available
	)

AC_SUBST(gsl_cflags)
AC_SUBST(gsl_libs)
AC_SUBST(gqr_cflags)
AC_SUBST(gqr_libs)
AC_SUBST(gts_cflags)
AC_SUBST(gts_libs)
AC_SUBST(libmatheval_CFLAGS)
AC_SUBST(libmatheval_LIBS)
AC_SUBST(sisl_cflags)
AC_SUBST(sisl_libs)
AC_SUBST(wmpi_cflags)
AC_SUBST(wmpi_libs)
AC_SUBST(mpi_tools)
AC_SUBST(amm_tools)

CFLAGS="$CFLAGS $libmatheval_CFLAGS $gts_cflags $gqr_cflags $gsl_cflags"
LDFLAGS="$LDFLAGS $gts_libs $gqr_libs $gsl_libs"

AC_SUBST(CFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(LDFLAGS)

AC_CHECK_LIB(m, sqrt)
#AC_CHECK_LIB(matheval, evaluator_evaluate)

LIBS="$gqr_libs $gsl_libs $libmatheval_LIBS $gts_libs $LIBS"

dnl generate `src/bem3dconfig.h'
AC_CONFIG_COMMANDS([default-1],[[
	echo creating src/bem3dconfig.h
	rm -f src/bem3dconfig.h
	cat >src/bem3dconfig.h <<\__EOF
/* bem3dconfig.h
 * 
 * This is a generated file.  Modify `configure.ac'
 */

#ifndef BEM3DCONFIG_H
#define BEM3DCONFIG_H

__EOF
	case "$build" in
	*-cray-unicos*)
	  echo "/* CRAY inlining directives */" >> src/bem3dconfig.h
	  awk '
BEGIN {
  nc = 10;
  nf = 0;
}
{
  if ($1 == "G_INLINE_FUNC") {
    nc = 0;
    start = 2;
  }
  else
    start = 1;
  for (i = start; i <= NF; i++) {
    if (nc == 1) {
      if (substr ($i, 1, 1) != "*") {
	in_list = 0;
	for (j = 0; j < nf && !inlist; j++)
	  if ($i == f[j])
	    inlist = 1;
	if (!inlist)
          f[nf++] = $i;
	nc = 10;
      }
    }
    else
      nc++;
  }
}
END {
  if (nf > 0) {
    printf ("#pragma _CRI inline %s", f[0]);
    for (i = 1; i < nf; i++)
      printf (", %s", f[i]);
    printf ("\n");
  }
}
' < src/bem3d.h >> src/bem3dconfig.h
	;;
	esac
	cat >>src/bem3dconfig.h <<__EOF
$bem3d_defines 
$matheval_defines

#endif /* BEM3DCONFIG_H */
__EOF
]],[[
	bem3d_defines="
#define BEM3D_MAJOR_VERSION $BEM3D_MAJOR_VERSION
#define BEM3D_MINOR_VERSION $BEM3D_MINOR_VERSION
#define BEM3D_MICRO_VERSION $BEM3D_MICRO_VERSION"
	matheval_defines="$matheval_defines"
	build="$build"
]])

AC_CONFIG_FILES([
Makefile
src/Makefile
src/bem3d-config
tools/Makefile
doc/Makefile
doc/html/Makefile
mpitools/Makefile
mpitools/bem3d-pbs-run])

AC_OUTPUT
