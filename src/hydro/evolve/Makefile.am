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
##	src/hydro/evolve/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

bin_PROGRAMS = chydro_leapfrog hydro_leapfrog

chydro_leapfrog_CPPFLAGS = -DTOOLBOX
hydro_leapfrog_CPPFLAGS = -DTOOLBOX

chydro_leapfrog_LDADD = ../libhydro.la ../../node/dyn/libdyn.la
hydro_leapfrog_LDADD = ../libhydro.la ../../node/dyn/libdyn.la

chydro_leapfrog_SOURCES = chydro_leapfrog.C
hydro_leapfrog_SOURCES = hydro_leapfrog.C
