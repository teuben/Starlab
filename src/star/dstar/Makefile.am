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

SUBDIRS = starclass init io stardyn util evolve

if USE_SDYN3
    SUBDIRS += sdyn3
endif

# Library to build and install:

pkglib_LTLIBRARIES = libdstar.la

# Local sub-libraries used to build library:

libdstar_la_LIBADD =	starclass/libdstar-starclass.la \
			init/libdstar-init.la \
			io/libdstar-io.la \
			stardyn/libdstar-stardyn.la \
			util/libdstar-util.la
CLIBS = clibs1

if USE_SDYN3
    libdstar_la_LIBADD += sdyn3/libdstar-sdyn3.la
    CLIBS += clibs2
endif

libdstar_la_LDFLAGS = ../sstar/libsstar.la
libdstar_la_SOURCES =

# Standard targets:

.PHONY: libs clibs1 clibs2
libs: $(CLIBS) libdstar.la

clibs1:
	$(MAKE) $(AM_MAKEFLAGS) -C starclass libdstar-starclass.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libdstar-init.la
	$(MAKE) $(AM_MAKEFLAGS) -C io libdstar-io.la
	$(MAKE) $(AM_MAKEFLAGS) -C stardyn libdstar-stardyn.la
	$(MAKE) $(AM_MAKEFLAGS) -C util libdstar-util.la

# The following needs sdyn3, hence the separation into a separate rule.

clibs2:
	$(MAKE) $(AM_MAKEFLAGS) -C sdyn3 libdstar-sdyn3.la
