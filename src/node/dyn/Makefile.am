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
##	src/node/dyn/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util kepler init evolve _dyn_ hdyn tdyn

if !NO_X
   SUBDIRS += xutil
endif

pkglib_LTLIBRARIES = libdyn.la

libdyn_la_LDFLAGS = ../libnode.la
libdyn_la_LIBADD = util/libdyn-1.la kepler/libdyn-2.la init/libdyn-3.la \
    evolve/libdyn-4.la $(LTLIBOBJS)
libdyn_la_SOURCES =

.PHONY: libs clibs
libs: clibs libdyn.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libdyn-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C kepler libdyn-2.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libdyn-3.la
	$(MAKE) $(AM_MAKEFLAGS) -C evolve libdyn-4.la
