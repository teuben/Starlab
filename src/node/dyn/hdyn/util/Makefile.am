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
##	src/node/dyn/hdyn/util/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libhdyn-1.la
bin_PROGRAMS = hdyn_pp3 reset_counters sys_stats scale get_densities

libhdyn_1_la_SOURCES = hdyn_io.C hdyn_init.C hdyn_kepler.C hdyn_pp3.C \
    hdyn_tt.C kira_diag.C kira_options.C reset_counters.C

hdyn_pp3_CPPFLAGS = -DTOOLBOX
reset_counters_CPPFLAGS = -DTOOLBOX
sys_stats_CPPFLAGS = -DTOOLBOX
scale_CPPFLAGS = -DTOOLBOX
get_densities_CPPFLAGS = -DTOOLBOX

hdyn_pp3_LDADD = ../libhdyn.la
reset_counters_LDADD = ../libhdyn.la
sys_stats_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la \
    ../stubs/libstubs.la
scale_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la \
    ../stubs/libstubs.la
get_densities_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la \
    ../stubs/libstubs.la

hdyn_pp3_SOURCES = hdyn_pp3.C
reset_counters_SOURCES = reset_counters.C
sys_stats_SOURCES = sys_stats.C
scale_SOURCES = scale.C
get_densities_SOURCES = get_densities.C
