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
##	src/star/sstar/init/Makefile.am
##..............................................................................
##	Dec 1993		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libsstar-2.la
bin_PROGRAMS = add_star

libsstar_2_la_SOURCES = add_star.C

add_star_CPPFLAGS = -DTOOLBOX

add_star_LDADD = ../libsstar.la ../../../node/dyn/hdyn/libhdyn.la

add_star_SOURCES = add_star.C
