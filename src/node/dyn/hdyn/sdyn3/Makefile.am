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
##	src/node/dyn/hdyn/sdyn3/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util evolve

pkglib_LTLIBRARIES = libsdyn3.la

libsdyn3_la_LDFLAGS = ../../_dyn_/lib_dyn_.la
libsdyn3_la_LIBADD = util/libsdyn3-1.la evolve/libsdyn3-2.la
libsdyn3_la_SOURCES =

.PHONY: libs clibs
libs: clibs libsdyn3.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libsdyn3-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C evolve libsdyn3-2.la
