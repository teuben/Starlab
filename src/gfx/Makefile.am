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

SUBDIRS = lux lux++ util #libs draw

pkglib_LTLIBRARIES = liblux.la

liblux_la_LDFLAGS = $(X_PRE_LIBS) -lX11 $(X_EXTRA_LIBS) -lm
if HAVE_LIBPNG
  liblux_la_LDFLAGS += -lpng
endif
liblux_la_LIBADD = lux/liblux-1.la lux++/liblux-2.la util/liblux-3.la
liblux_la_SOURCES =

AM_LDFLAGS = $(X_LIBS)

.PHONY: libs clibs
libs: clibs liblux.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C lux liblux-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C lux++ liblux-2.la
	$(MAKE) $(AM_MAKEFLAGS) -C util liblux-3.la
