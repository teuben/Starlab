# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
m4_include([version.m4])
AC_INIT(Starlab, VERSION_NUMBER, [steve@physics.drexel.edu])

all_args="$@"		# seems to pick up the proper args
			# too soon, ignired, too late, args changed

AC_CONFIG_AUX_DIR(config)
# doesn't really matter because there is a bug in Autoconf that disregards this
#AC_CONFIG_LIBOBJ_DIR(lib)
AC_CONFIG_MACRO_DIR(local)
AC_CONFIG_HEADERS(include/config.h)
AC_CONFIG_SRCDIR(src/std/util_math.C)
AC_LANG([C++])
AC_PREFIX_DEFAULT([$PWD/usr])
AM_INIT_AUTOMAKE([1.8 gnu no-exeext nostdinc std-options])

#--------------------------------------------------------------------------
# HUGE FIXME - $PWD will fail if building outside the source directory
# but I don't quite understand how ac_top_srcdir/top_srcdir work
CPPFLAGS="$CPPFLAGS -I$PWD/include -I$PWD/include/star"

# Another HUGE FIXME - pgetopt() must be replaced by getopt_long()
CPPFLAGS="$CPPFLAGS -D_SRC_='\"'\`pwd\`/'\"'__FILE__"

#--------------------------------------------------------------------------
# Checks for programs.
AC_PROG_CC
AC_PROG_F77
if test "x$ac_cv_prog_f77" = "xyes"; then
    AC_F77_WRAPPERS
fi
AM_CONDITIONAL(HAVE_F77, [test "x$ac_cv_prog_f77" = "xyes"])
AC_PROG_CXX
# disable shared libs for now -- work this out later
AC_DISABLE_SHARED
AC_PROG_LIBTOOL
LIBTOOL="$LIBTOOL --preserve-dup-deps"
AC_SUBST(LIBTOOL)
AC_PROG_LEX
AC_PROG_YACC
AC_PROG_INSTALL

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

#==========================================================================
#==========================================================================
# Add our own conditional stuff **after** the standard checks.  Some of
# our conditionals may force other checks to be performed, but the parser
# isn't smart enough to realize that these checks are only conditionally
# carried out, and so the script may fail to perform them.  (Steve, 9/04)

#--------------------------------------------------------------------------
# Source and install paths will be used in starlab_setup.* scripts.

AC_SUBST(STARLAB_SOURCE_PATH, $PWD)
AC_SUBST(STARLAB_VERSION_NUMBER, VERSION_NUMBER)

# From Steve, 09/03/04:  Configure will have set prefix=NONE by the time
# we get here if we use the default.
if test "x$prefix" = "xNONE" ; then
    AC_SUBST(STARLAB_INSTALL_PATH, $ac_default_prefix)
else
    AC_SUBST(STARLAB_INSTALL_PATH, $prefix)
fi

#--------------------------------------------------------------------------
# xreal setup.

AC_ARG_ENABLE(xreal,
    AS_HELP_STRING([--enable-xreal], use the xreal type [[[default=no]]]),
    xreal=yes, xreal=no)

if test "x$xreal" = "xyes" ; then

#   xreal is desired, but we can't do it if long long isn't defined or
#   has length other than 8 bytes.  Check these here and modify xreal
#   if necessary.

#   First check that long long exists.

    AC_CHECK_TYPE(long long, llong=yes, llong=no)

    if test "x$llong" != "xyes" ; then
	AC_MSG_WARN([ ])
	AC_MSG_WARN([* unsetting xreal because long long is not defined *])
	AC_MSG_WARN([ ])
 	xreal=no
    else
	AC_DEFINE(HAVE_LONG_LONG, 1,
		  [Define iff the long long type exists.])
    fi
    if test "x$xreal" = "xyes" ; then
	AC_CHECK_SIZEOF(long long)
	if test "x$ac_cv_sizeof_long_long" != "x8" ; then
	    AC_MSG_WARN([ ])
	    AC_MSG_WARN([* unsetting xreal because long long is not 8 bytes *])
	    AC_MSG_WARN([ ])
 	    xreal=no
	fi
    fi
    if test "x$xreal" = "xyes" ; then
        AC_DEFINE(USE_XREAL, 1,
        [Define to 1 if you want to enable the xreal type. Although this significantly increases precision, it also incurs a substantial performance penalty.])
    fi
fi

#--------------------------------------------------------------------------
# sdyn* setup.

AC_ARG_ENABLE(sdyn3,
    AS_HELP_STRING([--enable-sdyn3], compile the sdyn3 tree [[[default=no]]]),
    sdyn3=yes, sdyn3=no)
AM_CONDITIONAL(USE_SDYN3, [test "x$sdyn3" = "xyes"])

AC_ARG_ENABLE(sdyn,
    AS_HELP_STRING([--enable-sdyn], compile the sdyn tree [[[default=no]]]),
    sdyn=yes, sdyn=no)
AM_CONDITIONAL(USE_SDYN, [test "x$sdyn" = "xyes"])

#--------------------------------------------------------------------------
# Check for GRAPE libraries.

AC_ARG_WITH(grape,
    AS_HELP_STRING([--with-grape], use the GRAPE [[[default=yes]]]),
    [grape="$withval"])
if test -z "$grape" -o "x$grape" = "xyes"; then
    if test -f local/grape.sh; then
        . local/grape.sh
    elif test "x$grape" = "xyes"; then
        AC_MSG_ERROR([GRAPE was selected but configuration is missing.
Please put relevant information in `local/grape.sh'
as shown in `local/grape.sh.examples'.])
    else
        GRAPE_LDFLAGS_=
	GRAPE_LIBS_=
	OLD_READ_NEIGHBOUR_LIST=
    fi
    AC_SUBST(GRAPE_LDFLAGS, "$GRAPE_LDFLAGS_")
    AC_SUBST(GRAPE_LIBS, "$GRAPE_LIBS_")
    OLDFLAGS="$LDFLAGS"
    LDFLAGS="$LDFLAGS $GRAPE_LDFLAGS_"
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
Please make sure correct paths are specified in `local/grape.sh'.
See `config.log' for more details.])
        fi
    else
        AC_DEFINE(SHARE_GRAPE, 1, [Define to 1 if you want to release the GRAPE hardware periodically.])
	grape=yes
    fi
    LDFLAGS="$OLDFLAGS"
fi
AM_CONDITIONAL(HAVE_GRAPE, [test "x$grape" = "xyes"])
if test x"$OLD_READ_NEIGHBOUR_LIST" = "xyes"; then
    AC_DEFINE(G6_OLD_READ_NEIGHBOUR_LIST, 1, [Define to 1 to use the old (PIO) neighbor list interface.])
fi

#--------------------------------------------------------------------------
# Check for png.

AC_ARG_WITH(png, AS_HELP_STRING([--with-png], use PNG format [[[default=yes]]]),
    [png="$withval"])
if test -z "$png" -o "x$png" = "xyes"; then
    AC_CHECK_LIB(png, main, [AC_DEFINE(HAVE_LIBPNG, 1,
        [Define to 1 if you have the `png' library (-lpng).])])
fi
AM_CONDITIONAL(HAVE_LIBPNG, [test "x$ac_cv_lib_png_main" = "xyes"])

#--------------------------------------------------------------------------
# Check for Qt.

AC_ARG_WITH(qt, AS_HELP_STRING([--with-qt], use Qt [[[default=yes]]]),
    [qt="$withval"])
AC_ARG_WITH(qt-dir,
    AS_HELP_STRING([--with-qt-dir=QTDIR], directory where Qt is installed),
    [qt_dir="$withval"])

if test -z "$qt" -o "x$qt" = "xyes"; then

    good_qt=yes

    if test -z "$qt_dir"; then
        AC_MSG_CHECKING(whether QTDIR environment variable is set)
	if test "x$QTDIR" = "x"; then
	    AC_MSG_RESULT(no)
	    if test "x$qt" = "xyes"; then
	       AC_MSG_ERROR([Qt selected with --with-qt but QTDIR not set.
Please set environment QTDIR or use --with-qt-dir=QTDIR.])
	    fi
	    AC_MSG_NOTICE(*** disabling Qt ***)
	    good_qt=no
	else
	    AC_MSG_RESULT([$QTDIR])
	    qt_dir="$QTDIR"
        fi
    fi

    if test "x$good_qt" = "xyes"; then

#	We want Qt and have a Qt directory name.  Check that the
#	library is usable.  Look in both lib and lib64 if necessary.

	AC_SUBST(MOC, [$qt_dir/bin/moc])
	CPPFLAGS="$CPPFLAGS -I$qt_dir/include"

	found_qt=no
	save_LDFLAGS="$LDFLAGS"
	if test -d "$qt_dir/lib"; then
	    LDFLAGS="$LDFLAGS -L$qt_dir/lib"
	    AC_CHECK_LIB(qt, main, [found_qt="yes"], [found_qt="no"])
	fi

	if test "x$found_qt" = "xno"; then
	    if test -d "$qt_dir/lib64"; then
		LDFLAGS="$save_LDFLAGS -L$qt_dir/lib64"
		AC_MSG_NOTICE(trying lib64...)
		unset ac_cv_lib_qt_main		# (avoid possible caching)
		AC_CHECK_LIB(qt, main, [found_qt="yes"], [found_qt="no"])
	    fi
	fi

	good_qt=$found_qt

	if test "x$good_qt" = "xyes"; then

#	    Library is valid, but we need version 3 of libqt.

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
			AC_MSG_ERROR([Qt selected, but older than version 3
Please update Qt, or do not enable it with --with-qt.])
		    else
			AC_MSG_NOTICE(*** disabling Qt ***)
			good_qt=no
		    fi])

	    if test "x$good_qt" = "xyes"; then

#		Everything seems OK, finally add the Qt configuration.

		AC_DEFINE(HAVE_LIBQT, 1,
		      [Define to 1 if you have the `qt' library (-lqt).])
	    fi

	elif test "x$qt" = "xyes"; then

    	    AC_MSG_ERROR([Qt was selected but a test program failed to compile.
See `config.log' for more details.])

	else

	    AC_MSG_NOTICE(*** unable to compile test program: disabling Qt ***)
	    good_qt=no

	fi
    fi
fi

# For conditional builds in various Makefiles.
AM_CONDITIONAL(HAVE_LIBQT, [test "x$good_qt" = "xyes"])

#--------------------------------------------------------------------------
# Checks for library functions.
AC_REPLACE_FUNCS(stpcpy)

#--------------------------------------------------------------------------
AC_CONFIG_FILES([Makefile starlab_start.csh starlab_start.sh \
    starlab.pl VERSION include/Makefile \
    src/Makefile src/std/Makefile src/gfx/Makefile \
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
    src/star/dstar/util/Makefile src/star/dstar/sdyn3/Makefile \
    src/star/dstar/evolve/Makefile \
    src/star/util/Makefile src/star/xutil/Makefile src/hydro/Makefile \
    src/hydro/util/Makefile src/hydro/init/Makefile src/hydro/evolve/Makefile \
    src/harpsims/Makefile])

#--------------------------------------------------------------------------
# From Steve, 09/03/04:  We would like to make a handy link to the appropriate
# starlab_start file, but it is unclear how to distinguish between various
# shells.  The $SHELL variable will apparently always be /bin/sh here.  We need
# something that identifies the shell used by the user running configure.  If
# such a variable could be identified, then the following code would work...

# Assume that if $SHELL is defined and contains csh we are a C-shell or
# something similar.  Anything else is treated as a Bourne shell.
#
#shell=sh
#if test "x$SHELL" != "x" ; then
#    csh=`echo $SHELL | fgrep csh`
#    if test "x$csh" != "x" ; then
#	shell=csh
#    fi
#fi
#
# ...but it doesn't, so just link to the csh version (which is consistent
# with previous usage in any case.)

shell=csh

AC_MSG_NOTICE([linking starlab_start to starlab_start.$shell])

cd $STARLAB_SOURCE_PATH
ln -f -s starlab_start.$shell starlab_start

# Convenient to save a record of where configure was run.

AC_MSG_NOTICE([saving hostname and configure arguments in CONFIGURE_HOST])
echo Starlab configured on `hostname` on `date`		> CONFIGURE_HOST
echo "configure arguments: $all_args"		       >> CONFIGURE_HOST

#==========================================================================

AC_OUTPUT
