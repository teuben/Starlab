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
##	src/node/dyn/hdyn/evolve/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libhdyn-2.la
bin_PROGRAMS = kira

libhdyn_2_la_SOURCES = hdyn_unpert.C hdyn_ev.C hdyn_tree.C hdyn_merge.C \
    hdyn_schedule.C hdyn_schedule2.C hdyn_slow.C kira_init.C kira_config.C \
    kira_runtime.C perturbed_list.C correct_perturbers.C kira_check.C \
    kira_counters.C kira_energy.C kira_escape.C kira_ev.C kira_log.C \
    kira_stellar.C kira_approx.C kira_encounter.C kira_top_level.C kira_int.C \
    kira_density.C hdyn_grape4.C hdyn_grape6.C kira.C

kira_CPPFLAGS = -DTOOLBOX

kira_LDADD = ../libhdyn.la ../xutil/libxhdyn.la \
	     ../../../../star/dstar/libdstar.la

kira_SOURCES = kira.C

if HAVE_GRAPE
  kira_LDADD += $(GRAPE_LIBS)
  AM_LDFLAGS = $(GRAPE_LDFLAGS)
endif

kira_LDADD += ../stubs/libstubs.la
