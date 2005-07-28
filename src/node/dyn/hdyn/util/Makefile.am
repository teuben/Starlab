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
bin_PROGRAMS = hdyn_pp3 reset_counters hsys_stats hscale get_densities \
		hdyn_io

libhdyn_1_la_SOURCES = hdyn_io.C hdyn_init.C hdyn_kepler.C hdyn_pp3.C \
    hdyn_tt.C kira_diag.C kira_options.C reset_counters.C

hdyn_pp3_CPPFLAGS = -DTOOLBOX
reset_counters_CPPFLAGS = -DTOOLBOX
hsys_stats_CPPFLAGS = -DTOOLBOX
hscale_CPPFLAGS = -DTOOLBOX
get_densities_CPPFLAGS = -DTOOLBOX
hdyn_io_CPPFLAGS = -DTOOLBOX

hdyn_pp3_LDADD = ../libhdyn.la
reset_counters_LDADD = ../libhdyn.la
hsys_stats_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la
hscale_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la
get_densities_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la
hdyn_io_LDADD = ../libhdyn.la ../../../../star/dstar/libdstar.la

hdyn_pp3_SOURCES = hdyn_pp3.C
reset_counters_SOURCES = reset_counters.C
hsys_stats_SOURCES = hsys_stats.C
hscale_SOURCES = hscale.C
get_densities_SOURCES = get_densities.C
hdyn_io_SOURCES = hdyn_io.C

if HAVE_GRAPE
  hsys_stats_LDADD += $(GRAPE_LIBS)
  hscale_LDADD += $(GRAPE_LIBS)
  get_densities_LDADD += $(GRAPE_LIBS)
  AM_LDFLAGS = $(GRAPE_LDFLAGS)
endif
hsys_stats_LDADD += ../stubs/libstubs.la
hscale_LDADD += ../stubs/libstubs.la
get_densities_LDADD += ../stubs/libstubs.la
