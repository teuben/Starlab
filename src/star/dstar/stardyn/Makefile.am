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
##	src/star/dstar/stardyn/Makefile.am
##..............................................................................
##	Oct 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libdstar-4.la
bin_PROGRAMS = adddouble dstar_state dstar_stats

libdstar_4_la_SOURCES = adddouble.C dstar_state.C dstar_stats.C dstar_to_kira.C

adddouble_CPPFLAGS = -DTOOLBOX
dstar_state_CPPFLAGS = -DTOOLBOX
dstar_stats_CPPFLAGS = -DTOOLBOX

adddouble_LDADD = ../libdstar.la ../../../node/dyn/hdyn/libhdyn.la
dstar_state_LDADD = ../libdstar.la ../../../node/dyn/hdyn/libhdyn.la
dstar_stats_LDADD = ../libdstar.la ../../../node/dyn/hdyn/libhdyn.la

adddouble_SOURCES = adddouble.C
dstar_state_SOURCES = dstar_state.C
dstar_stats_SOURCES = dstar_stats.C
