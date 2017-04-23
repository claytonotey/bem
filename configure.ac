AC_INIT([bem],[0.0.0])

dnl Use the m4/ directory to contain libtool macros that will be needed to 
dnl go with the ltmain.sh script (both will be created/updated via the command
dnl libtoolize --copy --force
AC_CONFIG_MACRO_DIR([m4])

#------------------------------------------------------------------------------------
# Library's "official" release version:

LIBBEM_RELEASE_MAJOR=0
LIBBEM_RELEASE_MINOR=0
LIBBEM_RELEASE_BUILD=0
AC_CONFIG_SRCDIR([src/transfer.cpp])
AC_CANONICAL_TARGET([])
AM_CONFIG_HEADER(src/config.h)

AC_DISABLE_SHARED dnl allows for optimizations

#------------------------------------------------------------------------------------
# The following is the libtool / shared library version. This doesn't have to
# do anything with the release version. It MUST conform to the following rules:
#
#  1. Start with version information of `0:0:0' for each libtool library.
#  2. Update the version information only immediately before a public release of
#     your software. More frequent updates are unnecessary, and only guarantee
#     that the current interface number gets larger faster.
#  3. If the library source code has changed at all since the last update, then
#     increment revision (`c:r:a' becomes `c:r+1:a').
#  4. If any interfaces have been added, removed, or changed since the last update,
#     increment current, and set revision to 0.
#  5. If any interfaces have been added since the last public release, then increment
#     age.
#  6. If any interfaces have been removed since the last public release, then set age
#     to 0.

SHARED_VERSION_INFO="0:0:0"

AC_PROG_LIBTOOL
AM_PROG_LIBTOOL
AC_LANG([C++])

# save $CXXFLAGS etc. since AC_PROG_CXX likes to insert "-g -O2"
# if $CXXFLAGS is blank and it finds GCC
cflags_save="$CFLAGS"
cppflags_save="$CPPFLAGS"
cxxflags_save="$CXXFLAGS"
AC_LANG([C++])
AC_PROG_CXX

# restore those variables back
CFLAGS="$cflags_save"
CPPFLAGS="$cppflags_save"
CXXFLAGS="$cxxflags_save"

AC_SUBST(SHLIB_VERSION_ARG)
AC_SUBST(SHARED_VERSION_INFO)

AC_ARG_ENABLE(static,[  --enable-static enable static build: (default: disable)],[enable_static=$enableval],[enable_static=no])

AC_ARG_ENABLE(debug,[  --enable-debug enable debug build: (default: disable)],[enable_debug=$enableval],[enable_debug=no])

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
         [compile with MPI (parallelization) support. If none is found,
             MPI is not used. Default: auto])
        ],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[have_mpi=yes],[
        have_mpi=no
        if test x"$with_mpi" = xyes; then
            AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
        else
          AC_MSG_WARN([No MPI compiler found, won't use MPI.])
        fi
    ])

AC_ARG_WITH(omp, [AS_HELP_STRING([--with-omp],
            [compile with OpenMP support. If none is found,
                OpenMP is not used. Default: auto])
        ],,[with_omp=auto])

AX_OPENMP([have_omp=yes],[
        have_omp=no
        if test x"$with_omp" = xyes; then
          AC_MSG_FAILURE([OMP requested, but couldn't use OMP.])
        else
          AC_MSG_WARN([No OMP found, won't use OMP.])
        fi
    ])

PKG_CHECK_MODULES(XERCES, xerces-c >= 3.1.0, ac_cv_xerces=1, ac_cv_xerces=0)


if test x$enable_sndfile = xyes; then
   # check for presence of xerces-c
   if test "$ac_cv_xerces" = "0"; then
      echo "Sorry, you need libxerces-c";
      exit -1;
   else
      AC_SUBST(XERCES_CFLAGS)
      AC_SUBST(XERCES_LIBS)
      AC_DEFINE(HAVE_XERCES,1,[Define to use xerces-c])
   fi
fi

if test x$enable_static = xno; then
  case "$target_os" in
     darwin*)
	BEM_CFLAGS="$BEM_CFLAGS -mdynamic-no-pic"
	;;
	*)
	;;
   esac
fi

if test x$enable_debug = xyes; then
   BEM_CFLAGS="$BEM_CFLAGS -g"
else
   BEM_CFLAGS="$BEM_CFLAGS -O3 -fomit-frame-pointer"
fi

BEM_CFLAGS="$BEM_CFLAGS -I./util"

if test x$have_omp = xyes; then
   BEM_CFLAGS="$BEM_CFLAGS -fopenmp"
   AC_DEFINE(HAVE_OMP,1,[Define to use OpenMP])
fi

if test x$have_mpi = xyes; then
   AC_DEFINE(HAVE_MPI,1,[Define to use MPI])
fi

AM_CONDITIONAL(HAVE_MPI, test x$have_mpi = xyes)

AC_SUBST(BEM_CFLAGS)

BEM_PROGS=transfer
if test x$use_mpi = xyes; then
   BEM_PROGS="$BEM_PROGS transfermpi"
fi

AC_SUBST(BEM_PROGS)

AM_INIT_AUTOMAKE([])

AM_MAINTAINER_MODE()

AC_CONFIG_FILES([Makefile src/Makefile src/util/Makefile libbem.spec bem.pc])

AC_OUTPUT
