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
##	src/node/dyn/hdyn/sdyn/evolve/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

pkglib_LTLIBRARIES = libscat.la
bin_PROGRAMS = low_n_evolve scatter scatter_exp sigma tree_evolve

libscat_la_LIBADD = ../libsdyn.la
libscat_la_SOURCES = low_n_evolve.C sdyn_ev.C tree_evolve.C sdyn_ev.C \
    scatter_exp.C scatter_hist.C tree_evolve.C sig_helper.C scatter_merge.C \
    extend_or_end.C specific_function.C scatter.C multiscatter_MPI.C MPI_sigma.C

low_n_evolve_CPPFLAGS = -DTOOLBOX
scatter_CPPFLAGS = -DTOOLBOX
scatter_exp_CPPFLAGS = -DTOOLBOX
sigma_CPPFLAGS = -DTOOLBOX
tree_evolve_CPPFLAGS = -DTOOLBOX

low_n_evolve_LDADD = libscat.la
scatter_LDADD = libscat.la
scatter_exp_LDADD = libscat.la
sigma_LDADD = libscat.la
tree_evolve_LDADD = libscat.la

low_n_evolve_SOURCES = low_n_evolve.C
scatter_SOURCES = scatter.C
scatter_exp_SOURCES = scatter_exp.C
sigma_SOURCES = sigma.C
tree_evolve_SOURCES = tree_evolve.C
