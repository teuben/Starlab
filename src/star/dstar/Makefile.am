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
##..............................................................................
##

SUBDIRS = starclass init io stardyn util evolve

pkglib_LTLIBRARIES = libdstar.la

libdstar_la_LDFLAGS = ../sstar/libsstar.la
libdstar_la_LIBADD = starclass/libdstar-1.la init/libdstar-2.la \
    io/libdstar-3.la stardyn/libdstar-4.la util/libdstar-5.la
libdstar_la_SOURCES =

.PHONY: libs clibs
libs: clibs libdstar.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C starclass libdstar-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libdstar-2.la
	$(MAKE) $(AM_MAKEFLAGS) -C io libdstar-3.la
	$(MAKE) $(AM_MAKEFLAGS) -C stardyn libdstar-4.la
	$(MAKE) $(AM_MAKEFLAGS) -C util libdstar-5.la
