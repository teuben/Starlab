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
##	src/node/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

SUBDIRS = util dyn

pkglib_LTLIBRARIES = libnode.la

libnode_la_LDFLAGS = ../std/libstd.la
libnode_la_LIBADD = util/libnode-1.la
libnode_la_SOURCES =

.PHONY: libs clibs
libs: clibs libnode.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libnode-1.la
