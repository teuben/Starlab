# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Starlab, 4.3, [steve@physics.drexel.edu])
AC_CONFIG_AUX_DIR(config)
#AC_CONFIG_MACRO_DIR(local)    # not implemented yet as of 2.59
AC_CONFIG_HEADERS(include/config.h)
AC_CONFIG_SRCDIR(src/std/util_math.C)
AC_LANG([C++])
AC_PREFIX_DEFAULT([$PWD/usr])
AM_INIT_AUTOMAKE([1.8 gnu no-exeext nostdinc std-options])

# HUGE FIXME - $PWD will fail if building outside the source directory
# but I don't quite understand how ac_top_srcdir/top_srcdir work
CPPFLAGS="$CPPFLAGS -I$PWD/include -I$PWD/include/star"
AC_ARG_ENABLE(xreal, AS_HELP_STRING([--enable-xreal], use the xreal type),
	[AC_DEFINE(USE_XREAL, 1,
		[Define to 1 if you want to enable the xreal type. Although this significantly increases precision, it also incurs a substantial performance penalty.])])
# Another HUGE FIXME - pgetopt() must be replaced by getopt_long()
CPPFLAGS="$CPPFLAGS -D_SRC_='\"'\`pwd\`/'\"'__FILE__ -D_COMPILE_DATE_=__DATE__"

# Checks for programs.
AC_PROG_CC
AC_PROG_F77
AC_F77_WRAPPERS
AC_PROG_CXX
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
LIBTOOL="$LIBTOOL --preserve-dup-deps"
AC_SUBST(LIBTOOL)
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_INSTALL

# Checks for libraries.
AC_ARG_WITH(qt-dir,
        AS_HELP_STRING([--with-qt-dir], directory where Qt is installed),
        [qt_dir="$withval"])
if test -z "$qt_dir"; then
   AC_MSG_CHECKING(whether QTDIR environment variable is set)
   if test "x$QTDIR" = "x"; then
      AC_MSG_RESULT(no)
      AC_MSG_ERROR([QTDIR or --with-qt-dir must be properly set])
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
AC_CHECK_LIB(png, main, [AC_DEFINE(HAVE_LIBPNG, 1,
	[Define to 1 if you have the `png' library (-lpng).])])

# Checks for header files.
AC_LANG_PUSH(C)
AC_PATH_XTRA
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_CONST
AC_C_INLINE
AC_LANG_POP

# Checks for library functions.
AC_CHECK_FUNCS(stpcpy)

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
