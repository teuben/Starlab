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
##	src/node/dyn/hdyn/sdyn/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util init evolve

pkglib_LTLIBRARIES = libsdyn.la

libsdyn_la_LDFLAGS = ../libhdyn.la
libsdyn_la_LIBADD = util/libsdyn-1.la init/libsdyn-2.la
libsdyn_la_SOURCES =

.PHONY: libs clibs
libs: clibs libsdyn.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libsdyn-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C init libsdyn-2.la
