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
##	Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2001		Steve McMillan & Peter Teuben
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

pkginclude_HEADERS = include/bound3.h include/chydro.h include/config.h \
    include/c_stdinc.h include/diag.h include/_dyn_.h include/dyn.h \
    include/dyn_kepler.h include/grape4.h include/grape6.h include/hdyn.h \
    include/hdyn_kepler.h include/hydrobase.h include/hydro.h include/kepler.h \
    include/kira.h include/node.h include/pdyn.h include/pp3_widget.h \
    include/scatter3.h include/scatter.h include/sdyn3.h include/sdyn.h \
    include/sigma3.h include/slow_binary.h include/stamp-h.in \
    include/starbase.h include/star/black_hole.h include/star/brown_dwarf.h \
    include/star/cluster_support.h include/star/constants.h \
    include/star/constants_new.h include/star/disintegrated.h \
    include/star/double_star.h include/star/double_support.h \
    include/star/dstar_to_dyn.h include/star/dstar_to_kira.h \
    include/star/general_support.h include/star/giant.h \
    include/star/helium_giant.h include/star/helium_star.h \
    include/star/hertzsprung_gap.h include/star/horizontal_branch.h \
    include/star/hyper_giant.h include/starlab_vector.h \
    include/star/main_sequence.h include/star/neutron_star.h \
    include/star/planet.h include/star/proto_star.h include/star/seba.h \
    include/star/single_star.h include/star/SPZDCH_star.h \
    include/star/sstar_to_dyn.h include/star/star_cluster.h \
    include/star/stardyn_util.h include/star/star.h include/star/star_state.h \
    include/star/star_support.h include/star/stdfunc.h \
    include/star/stellar_type.h include/star/sub_giant.h \
    include/star/super_giant.h include/star/thorne_zytkow.h \
    include/star/white_dwarf.h include/static_star.h include/stdinc.h \
    include/story.h include/tdyn.h include/util_io.h include/util_math.h \
    include/worldline.h include/xreal.h

EXTRA_DIST = doc/Primer.tex doc/README doc/Testing.txt local sbin/preconf \
    scripts/kira_test scripts/kira_xdemo scripts/make_binary_demo \
    scripts/make_stellar_demo

SUBDIRS = src

all-recursive: libs

.PHONY: libs
libs:
	for dir in src/std src/gfx src/node src/node/dyn \
	    src/node/dyn/_dyn_ src/node/dyn/hdyn src/node/dyn/hdyn/stubs \
	    src/node/dyn/hdyn/sdyn src/node/dyn/hdyn/sdyn3 \
	    src/node/dyn/hdyn/xutil src/node/dyn/tdyn \
	    src/star/sstar src/star/dstar src/hydro; do \
	    $(MAKE) $(AM_MAKEFLAGS) -C $$dir libs; \
	done
