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
##	src/star/rdc/Makefile.am
##..............................................................................
##	Mar 2008	Ernest Mamikonyan
##..............................................................................
##

pkglib_LTLIBRARIES = librdc.la
bin_PROGRAMS = list_sdumb makepovfile rdc_SeBa rs_birthrate rs_countsn \
    rs_findtype rs_population rs_reorder rs_scenarios rs_snapshot
#make_ccd

librdc_la_SOURCES = SeBa_hist.C makepovfile.C rdc_SeBa.C

list_sdumb_CPPFLAGS = -DTOOLBOX
#make_ccd_CPPFLAGS = -DTOOLBOX
makepovfile_CPPFLAGS = -DTOOLBOX
rdc_SeBa_CPPFLAGS = -DTOOLBOX
rs_birthrate_CPPFLAGS = -DTOOLBOX
rs_countsn_CPPFLAGS = -DTOOLBOX
rs_findtype_CPPFLAGS = -DTOOLBOX
rs_population_CPPFLAGS = -DTOOLBOX
rs_reorder_CPPFLAGS = -DTOOLBOX
rs_scenarios_CPPFLAGS = -DTOOLBOX
rs_snapshot_CPPFLAGS = -DTOOLBOX


list_sdumb_LDADD = ../sstar/libsstar.la ../../node/dyn/hdyn/libhdyn.la
#make_ccd_LDADD = ../sstar/libsstar.la ../../node/dyn/libdyn.la
makepovfile_LDADD = librdc.la ../sstar/libsstar.la ../../node/dyn/libdyn.la
rdc_SeBa_LDADD = librdc.la ../dstar/libdstar.la ../../node/dyn/libdyn.la
rs_birthrate_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_countsn_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_findtype_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_population_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_reorder_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_scenarios_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la
rs_snapshot_LDADD = librdc.la ../dstar/libdstar.la ../../node/libnode.la

list_sdumb_SOURCES = list_sdumb.C
#make_ccd_SOURCES = make_ccd.C
makepovfile_SOURCES = makepovfile.C
rdc_SeBa_SOURCES = makepovfile.C
rs_birthrate_SOURCES = rs_birthrate.C
rs_countsn_SOURCES = rs_countsn.C
rs_findtype_SOURCES = rs_findtype.C
rs_population_SOURCES = rs_population.C
rs_reorder_SOURCES = rs_reorder.C
rs_scenarios_SOURCES = rs_scenarios.C
rs_snapshot_SOURCES = rs_snapshot.C
