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
##	src/star/sstar/stardyn/Makefile.am
##..............................................................................
##	Dec 1993		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libsstar-4.la
bin_PROGRAMS = sstar_stats

libsstar_4_la_SOURCES = sstar_to_dyn.C sstar_state.C stardyn_story.C \
    addstar_dyn.C sstar_stats.C proto_star_dyn.C

sstar_stats_CPPFLAGS = -DTOOLBOX

sstar_stats_LDADD = ../libsstar.la ../../../node/dyn/hdyn/libhdyn.la

sstar_stats_SOURCES = sstar_stats.C
