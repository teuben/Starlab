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
##	src/star/util/Makefile.am
##..............................................................................
##	Dec 1993		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libsstar-5.la
bin_PROGRAMS =  kingfit red_stellar_system rcore2 coagulation rearrange_mf \
    red_star_cluster encountervelocity imf_evolve proj_lagr lmt_to_ubvri \
    makedyn_from_input nn_distribution

libsstar_5_la_SOURCES = kingfit.C lmt_to_ubvri.C proj_lagr.C

kingfit_CPPFLAGS = -DTOOLBOX
red_stellar_system_CPPFLAGS = -DTOOLBOX
rcore2_CPPFLAGS = -DTOOLBOX
coagulation_CPPFLAGS = -DTOOLBOX
rearrange_mf_CPPFLAGS = -DTOOLBOX
red_star_cluster_CPPFLAGS = -DTOOLBOX
encountervelocity_CPPFLAGS = -DTOOLBOX
imf_evolve_CPPFLAGS = -DTOOLBOX
proj_lagr_CPPFLAGS = -DTOOLBOX
lmt_to_ubvri_CPPFLAGS = -DTOOLBOX
makedyn_from_input_CPPFLAGS = -DTOOLBOX
nn_distribution_CPPFLAGS = -DTOOLBOX

kingfit_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
red_stellar_system_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
rcore2_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
coagulation_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
rearrange_mf_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
red_star_cluster_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
encountervelocity_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
imf_evolve_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
proj_lagr_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
lmt_to_ubvri_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
makedyn_from_input_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la
nn_distribution_LDADD = ../dstar/libdstar.la ../../node/dyn/libdyn.la

kingfit_SOURCES = kingfit.C
red_stellar_system_SOURCES = red_stellar_system.C
rcore2_SOURCES = rcore2.C
coagulation_SOURCES = coagulation.C
rearrange_mf_SOURCES = rearrange_mf.C
red_star_cluster_SOURCES = red_star_cluster.C
encountervelocity_SOURCES = encountervelocity.C
imf_evolve_SOURCES = imf_evolve.C
proj_lagr_SOURCES = proj_lagr.C
lmt_to_ubvri_SOURCES = lmt_to_ubvri.C
makedyn_from_input_SOURCES = makedyn_from_input.C
nn_distribution_SOURCES = nn_distribution.C
