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
##	src/node/dyn/hdyn/sdyn/util/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libsdyn-1.la
bin_PROGRAMS = make_tree sdyn_io

libsdyn_1_la_SOURCES = make_tree.C sdyn_io.C

make_tree_CPPFLAGS = -DTOOLBOX
sdyn_io_CPPFLAGS = -DTOOLBOX

make_tree_LDADD = ../libsdyn.la
sdyn_io_LDADD = ../libsdyn.la

make_tree_SOURCES = make_tree.C
sdyn_io_SOURCES = sdyn_io.C
