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
##	src/star/sstar/evolve/Makefile.am
##..............................................................................
##	Mar 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

bin_PROGRAMS = starev nstarev

starev_CPPFLAGS = -DTOOLBOX
nstarev_CPPFLAGS = -DTOOLBOX

starev_LDADD = ../libsstar.la ../../../node/dyn/libdyn.la
nstarev_LDADD = ../libsstar.la ../../../node/dyn/libdyn.la

starev_SOURCES = starev.C
nstarev_SOURCES = nstarev.C
