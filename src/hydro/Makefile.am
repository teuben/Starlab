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
##	src/hydro/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util init evolve

pkglib_LTLIBRARIES = libhydro.la

libhydro_la_LIBADD = util/libhydro-1.la init/libhydro-2.la
libhydro_la_SOURCES =

.PHONY: libs clibs
libs: clibs libhydro.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libhydro-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libhydro-2.la
