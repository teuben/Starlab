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
##	src/star/dstar/Makefile.am
##..............................................................................
##	Mar 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##	Sep 2004		Annotated by Steve McMillan
##..............................................................................
##

SUBDIRS = starclass init io stardyn util sdyn3 evolve

# Library to build and install:

pkglib_LTLIBRARIES = libdstar.la

# Local sub-libraries used to build library:

libdstar_la_LIBADD =	starclass/libdstar-starclass.la \
			init/libdstar-init.la \
			io/libdstar-io.la \
			stardyn/libdstar-stardyn.la \
			util/libdstar-util.la \
			sdyn3/libdstar-sdyn3.la

libdstar_la_LDFLAGS = ../sstar/libsstar.la
libdstar_la_SOURCES =

# Standard targets:

.PHONY: libs clibs
libs: clibs libdstar.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C starclass libdstar-starclass.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libdstar-init.la
	$(MAKE) $(AM_MAKEFLAGS) -C io libdstar-io.la
	$(MAKE) $(AM_MAKEFLAGS) -C stardyn libdstar-stardyn.la
	$(MAKE) $(AM_MAKEFLAGS) -C util libdstar-util.la
	$(MAKE) $(AM_MAKEFLAGS) -C sdyn3 libdstar-sdyn3.la
