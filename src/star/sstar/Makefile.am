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
##	src/star/sstar/Makefile.am
##..............................................................................
##	Mar 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = starclass init io stardyn evolve

pkglib_LTLIBRARIES = libsstar.la

libsstar_la_LIBADD = starclass/libsstar-1.la init/libsstar-2.la \
    io/libsstar-3.la stardyn/libsstar-4.la ../util/libsstar-5.la
libsstar_la_SOURCES =

.PHONY: libs clibs
libs: clibs libsstar.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C starclass libsstar-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libsstar-2.la
	$(MAKE) $(AM_MAKEFLAGS) -C io libsstar-3.la
	$(MAKE) $(AM_MAKEFLAGS) -C stardyn libsstar-4.la
	$(MAKE) $(AM_MAKEFLAGS) -C ../util libsstar-5.la
