# 
#        #=======================================================#     _\|/_
#       #  __  _____           ___                    ___       #       /|\
#      #  /      |      ^     |   \  |         ^     |   \     #           _\|/_
#     #   \__    |     / \    |___/  |        / \    |___/    #             /|\
#    #       \   |    /___\   |  \   |       /___\   |   \   #  _\|/_
#   #     ___/   |   /     \  |   \  |____  /     \  |___/  #    /|\
#  #                                                       #             _\|/_
# #=======================================================#               /|\
#
##
##	src/gfx/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util
if !NO_X
    SUBDIRS += lux lux++ #libs draw
endif

pkglib_LTLIBRARIES = libgfx.la

if !NO_X
    libgfx_la_LDFLAGS = $(X_PRE_LIBS) -lX11 $(X_EXTRA_LIBS) -lm
else
    libgfx_la_LDFLAGS = -lm
endif

if HAVE_LIBPNG
    libgfx_la_LDFLAGS += -lpng
endif

libgfx_la_LIBADD = util/libgfx-1.la
if !NO_X
    libgfx_la_LIBADD += lux/libgfx-2.la lux++/libgfx-3.la
endif
libgfx_la_SOURCES =

AM_LDFLAGS = $(X_LIBS)

if !NO_X
    CLIBSX = clibs23
endif

.PHONY: libs clibs1 $(CLIBSX)
libs: clibs1 $(CLIBSX) libgfx.la

clibs1:
	$(MAKE) $(AM_MAKEFLAGS) -C util libgfx-1.la

clibs23:
	$(MAKE) $(AM_MAKEFLAGS) -C lux libgfx-2.la
	$(MAKE) $(AM_MAKEFLAGS) -C lux++ libgfx-3.la
