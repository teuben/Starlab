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
##	src/hydro/init/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libhydro-2.la
bin_PROGRAMS = addchydro addhydro

libhydro_2_la_SOURCES = addchydro.C addhydro.C

addchydro_CPPFLAGS = -DTOOLBOX
addhydro_CPPFLAGS = -DTOOLBOX

addchydro_LDADD = ../libhydro.la ../../node/libnode.la
addhydro_LDADD = ../libhydro.la ../../node/libnode.la

addchydro_SOURCES = addchydro.C
addhydro_SOURCES = addhydro.C
