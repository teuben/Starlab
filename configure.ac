# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
m4_include([version.m4])
AC_INIT(Starlab, VERSION_NUMBER, [steve@physics.drexel.edu])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_LIBOBJ_DIR(lib)
AC_CONFIG_MACRO_DIR(local)
AC_CONFIG_HEADERS(include/config.h)
AC_CONFIG_SRCDIR(src/std/util_math.C)
AC_LANG([C++])
AC_PREFIX_DEFAULT([$PWD/usr])
AM_INIT_AUTOMAKE([1.8 gnu no-exeext nostdinc std-options])

# HUGE FIXME - $PWD will fail if building outside the source directory
# but I don't quite understand how ac_top_srcdir/top_srcdir work
CPPFLAGS="$CPPFLAGS -I$PWD/include -I$PWD/include/star"

AC_ARG_ENABLE(xreal,
    AS_HELP_STRING([--enable-xreal], use the xreal type [[[default=no]]]),
    [AC_DEFINE(USE_XREAL, 1,
        [Define to 1 if you want to enable the xreal type. Although this significantly increases precision, it also incurs a substantial performance penalty.])])

AC_ARG_ENABLE(sdyn3,
    AS_HELP_STRING([--enable-sdyn3], compile the sdyn3 tree [[[default=no]]]),
    sdyn3=yes, sdyn3=no)
AM_CONDITIONAL(USE_SDYN3, [test "x$sdyn3" = "xyes"])

AC_ARG_ENABLE(sdyn,
    AS_HELP_STRING([--enable-sdyn], compile the sdyn tree [[[default=no]]]),
    sdyn=yes, sdyn=no)
AM_CONDITIONAL(USE_SDYN, [test "x$sdyn" = "xyes"])

# Another HUGE FIXME - pgetopt() must be replaced by getopt_long()
CPPFLAGS="$CPPFLAGS -D_SRC_='\"'\`pwd\`/'\"'__FILE__"

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_F77_WRAPPERS
AC_PROG_CXX
# disable shared libs for now - I'll work this out later
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
LIBTOOL="$LIBTOOL --preserve-dup-deps"
AC_SUBST(LIBTOOL)
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_INSTALL

# Checks for libraries.
m4_include([local/grape.m4])
AC_ARG_WITH(grape,
    AS_HELP_STRING([--with-grape], use the GRAPE [[[default=yes]]]),
    [grape="$withval"])
if test -z "$grape" -o "x$grape" = "xyes"; then
    AC_SUBST(GRAPE_LDFLAGS, "GRAPE_LDFLAGS_")
    AC_SUBST(GRAPE_LIBS, "GRAPE_LIBS_")
    OLDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $GRAPE_LDFLAGS"
    AC_CHECK_LIB(g6lx, g6_open_, grape6=yes, grape6=no)
    if test "x$grape6" = "xno"; then
        AC_CHECK_LIB(g4x, h3open_, grape4=yes, grape4=no)
	if test "x$grape4" = "xyes"; then
	    AC_DEFINE(USE_HALF_CHIP, 1, [Define to 1 to use half chip and preserve neighbor lists.])
	fi
    fi
    if test "x$grape6" = "xno" -a "x$grape4" = "xno"; then
        if test "x$grape" = "xyes"; then
            AC_MSG_ERROR([GRAPE was selected but a test program failed to compile.
Please make sure correct paths are specified in `local/grape.m4'.
See `config.log' for more details.])
        fi
    else
        AC_DEFINE(SHARE_GRAPE, 1, [Define to 1 if you want to release the GRAPE hardware periodically.])
	grape=yes
    fi
    LDFLAGS="$OLDFLAGS"
fi
AM_CONDITIONAL(HAVE_GRAPE, [test "x$grape" = "xyes"])
if test x"OLD_READ_NEIGHBOUR_LIST" = "xyes"; then
    AC_DEFINE(G6_OLD_READ_NEIGHBOUR_LIST, 1, [Define to 1 to use the old (PIO) neighbor list interface.])
fi

AC_ARG_WITH(png, AS_HELP_STRING([--with-png], use PNG format [[[default=yes]]]),
    [png="$withval"])
if test -z "$png" -o "x$png" = "xyes"; then
    AC_CHECK_LIB(png, main, [AC_DEFINE(HAVE_LIBPNG, 1,
        [Define to 1 if you have the `png' library (-lpng).])])
fi
AM_CONDITIONAL(HAVE_LIBPNG, [test "x$ac_cv_lib_png_main" = "xyes"])

AC_ARG_WITH(qt, AS_HELP_STRING([--with-qt], use Qt [[[default=yes]]]),
    [qt="$withval"])
AC_ARG_WITH(qt-dir,
    AS_HELP_STRING([--with-qt-dir=QTDIR], directory where Qt is installed),
    [qt_dir="$withval"])
if test -z "$qt" -o "x$qt" = "xyes"; then
    if test -z "$qt_dir"; then
        AC_MSG_CHECKING(whether QTDIR environment variable is set)
	if test "x$QTDIR" = "x"; then
	    AC_MSG_RESULT(no)
	    if test "x$qt" = "xyes"; then
	       AC_MSG_ERROR([Qt was selected with --with-qt but QTDIR isn't set.
Please set QTDIR or use --with-qt-dir.])
	    fi
	else
	    AC_MSG_RESULT([$QTDIR])
	    qt_dir="$QTDIR"
        fi
    fi
    CPPFLAGS="$CPPFLAGS -I$qt_dir/include"
    LDFLAGS="$LDFLAGS -L$qt_dir/lib"
    AC_SUBST(MOC, [$qt_dir/bin/moc])
    AC_CHECK_LIB(qt, main, [AC_DEFINE(HAVE_LIBQT, 1,
   	[Define to 1 if you have the `qt' library (-lqt).])])
    if test "x$ac_cv_lib_qt_main" = "xyes"; then
	AC_MSG_CHECKING(whether Qt is version 3 or newer)
        AC_PREPROC_IFELSE(
	    [AC_LANG_PROGRAM(
[#include <qglobal.h>
#if QT_VERSION < 0x030000
#  error Qt is too old
#endif])],
	    [AC_MSG_RESULT(yes); good_qt=yes],
	    [AC_MSG_RESULT(no);
	        if test "x$qt" = "xyes"; then
	            AC_MSG_ERROR([Qt was selected but it is older than version 3
Please update Qt, or do not enable it with --with-qt.])
	        else
		    good_qt=no
		fi])
    elif test "x$qt" = "xyes"; then
    	AC_MSG_ERROR([Qt was selected but a test program failed to compile
See `config.log' for more details.])
    fi
fi
AM_CONDITIONAL(HAVE_LIBQT, [test "x$good_qt" = "xyes"])

# Checks for header files.
AC_LANG_PUSH(C)
AC_PATH_XTRA
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_LANG_POP

# Checks for library functions.
AC_REPLACE_FUNCS(stpcpy)

AC_CONFIG_FILES([Makefile src/Makefile src/std/Makefile src/gfx/Makefile \
    src/gfx/lux/Makefile src/gfx/lux++/Makefile src/gfx/util/Makefile \
    src/node/Makefile src/node/util/Makefile src/node/dyn/Makefile \
    src/node/dyn/util/Makefile src/node/dyn/kepler/Makefile \
    src/node/dyn/init/Makefile src/node/dyn/evolve/Makefile \
    src/node/dyn/_dyn_/Makefile src/node/dyn/xutil/Makefile \
    src/node/dyn/hdyn/Makefile src/node/dyn/hdyn/util/Makefile \
    src/node/dyn/hdyn/evolve/Makefile src/node/dyn/hdyn/stubs/Makefile \
    src/node/dyn/hdyn/xutil/Makefile src/node/dyn/hdyn/sdyn/Makefile \
    src/node/dyn/hdyn/sdyn/util/Makefile src/node/dyn/hdyn/sdyn/init/Makefile \
    src/node/dyn/hdyn/sdyn/evolve/Makefile src/node/dyn/hdyn/sdyn3/Makefile \
    src/node/dyn/hdyn/sdyn3/util/Makefile \
    src/node/dyn/hdyn/sdyn3/evolve/Makefile src/node/dyn/tdyn/Makefile \
    src/star/Makefile src/star/sstar/Makefile \
    src/star/sstar/starclass/Makefile src/star/sstar/io/Makefile \
    src/star/sstar/init/Makefile src/star/sstar/stardyn/Makefile \
    src/star/sstar/evolve/Makefile src/star/dstar/Makefile \
    src/star/dstar/starclass/Makefile src/star/dstar/io/Makefile \
    src/star/dstar/init/Makefile src/star/dstar/stardyn/Makefile \
    src/star/dstar/util/Makefile src/star/dstar/evolve/Makefile \
    src/star/util/Makefile src/star/xutil/Makefile src/hydro/Makefile \
    src/hydro/util/Makefile src/hydro/init/Makefile src/hydro/evolve/Makefile \
    src/harpsims/Makefile])
AC_OUTPUT
