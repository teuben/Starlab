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

# Extra material for inclusion in the tarball, not explicitly mentioned
# in the various Makefiles.

EXTRA_DIST =	doc/Primer.tex doc/README doc/Testing.txt sbin/preconf \
		scripts/kira_test scripts/kira_xdemo scripts/make_binary_demo \
		scripts/make_stellar_demo \
		local/grape.sh.examples local/grape.sh.template \
		install-starlab \
		etc \
		EXAMPLES \
		sbin/add_options \
		sbin/beautify-install \
		sbin/build_report \
		sbin/bundle \
		sbin/cath \
		sbin/cath1 \
		sbin/change_option.awk \
		sbin/change_shell \
		sbin/change_shells_list \
		sbin/check_autosave0 \
		sbin/check_dep.c \
		sbin/check_display \
		sbin/check_make0 \
		sbin/create_system \
		sbin/diff_loop.sh \
		sbin/distribute_tar \
		sbin/do_big_tar \
		sbin/earliest_time.c \
		sbin/enforce_dep \
		sbin/enforce_dep.c \
		sbin/enforce_dep.csh \
		sbin/extrange \
		sbin/file_dir \
		sbin/get_tar \
		sbin/gpl \
		sbin/hide \
		sbin/index \
		sbin/kiraformat.c \
		sbin/kiraindex.c \
		sbin/kn \
		sbin/komaba_get_tar \
		sbin/largest_int.c \
		sbin/latest_time.c \
		sbin/latest_time.save \
		sbin/length \
		sbin/lengthof \
		sbin/loop \
		sbin/make_ar_cleanup \
		sbin/make_cond_build \
		sbin/make_cond_compile \
		sbin/make_db \
		sbin/make_dyn_select \
		sbin/make_export \
		sbin/Makefile \
		sbin/make_log \
		sbin/make_SeBa \
		sbin/make_touch_src \
		sbin/makewalk \
		sbin/makewalk.quiet \
		sbin/mkshadow \
		sbin/mkstarlab \
		sbin/mkwsl \
		sbin/mon2dmem \
		sbin/mon2mem \
		sbin/mtail \
		sbin/mytime.c \
		sbin/png2mng.pl \
		sbin/procmon \
		sbin/qt_major_version \
		sbin/README \
		sbin/readlink.c \
		sbin/remove_options \
		sbin/remove_starlab_bin \
		sbin/reset_list \
		sbin/rget \
		sbin/routines.c \
		sbin/send_tar \
		sbin/shadowdir \
		sbin/sinstall \
		sbin/show \
		sbin/skn \
		sbin/slow_tail \
		sbin/smake \
		sbin/sqrt.c \
		sbin/sremake \
		sbin/starlab_compare \
		sbin/starlab-do2.cfg \
		sbin/strindex \
		sbin/strings.c \
		sbin/strip_starlab_path.c \
		sbin/substr \
		sbin/substring \
		sbin/test.c \
		sbin/testsuite \
		sbin/trn.c \
		sbin/untar \
		sbin/update \
		sbin/vcp \
		sbin/view_core \
		sbin/view_cpu \
		sbin/view_dcpu \
		sbin/view_dde \
		sbin/view_de \
		sbin/view_lagr \
		sbin/view_n \
		sbin/wsl \
		sbin/wsll \
		sbin/wsll_node \
		sbin/wsll_star \
		sbin/wsl_node \
		sbin/wsl_star \
		src/gfx/draw \
		src/gfx/libs \
		src/gfx/Local_setup.inc \
		src/gfx/lux/data.c \
		src/gfx/lux/demo.c \
		src/gfx/lux/hydro2.c \
		src/gfx/lux/make_data.c \
		src/gfx/lux/simple_plot.h \
		src/gfx/lux/win.h \
		src/gfx/lux/ximage_compiled.c \
		src/gfx/lux/Xnames.h \
		src/gfx/lux++/c_lux.h \
		src/gfx/lux++/c_stdinc.h \
		src/gfx/lux++/local.h \
		src/gfx/Makefile.standalone \
		src/gfx/Make_tar \
		src/gfx/palettes \
		src/gfx/STANDALONE_README \
		src/gfx/util/image_fmt.h \
		src/gfx/util/make_data.C \
		src/gfx/util/make_header.C \
		src/gfx/util/old_make_image.C \
		src/gfx/util/test_byte.C \
		src/gfx/util/whirlgif.h \
		src/gfx/util/write_image.C \
		src/harpsims/harp3fastsim.cmp.c \
		src/harpsims/harp3fastsim-ok.c \
		src/harpsims/harp3futil.c \
		src/hydro/util/dummy.C \
		src/node/dyn/hdyn/evolve/Binaries.tex \
		src/node/dyn/hdyn/evolve/flat_hermite.C \
		src/node/dyn/hdyn/evolve/hdyn_flat.C \
		src/node/dyn/hdyn/evolve/hdyn_inline.C \
		src/node/dyn/hdyn/evolve/hpsort.c \
		src/node/dyn/hdyn/evolve/kira_debug.h \
		src/node/dyn/hdyn/evolve/kira_defaults.h \
		src/node/dyn/hdyn/evolve/kira_grape_README.txt \
		src/node/dyn/hdyn/evolve/kira_post_step_debug.C \
		src/node/dyn/hdyn/evolve/kira_pre_step_debug.C \
		src/node/dyn/hdyn/evolve/kira_sample_debug.C \
		src/node/dyn/hdyn/evolve/kira_timing.h \
		src/node/dyn/hdyn/evolve/piksrt.c \
		src/node/dyn/hdyn/evolve/README.steve \
		src/node/dyn/hdyn/evolve/SNAP_PARAMS.txt \
		src/node/dyn/hdyn/evolve/sort.c \
		src/node/dyn/hdyn/evolve/sqrt_tst.C \
		src/node/dyn/hdyn/evolve/stest.C \
		src/node/dyn/hdyn/sdyn/evolve/localmpi++.h \
		src/node/dyn/hdyn/sdyn/evolve/scatter_exp.h \
		src/node/dyn/hdyn/sdyn/evolve/scatter_hist.h \
		src/node/dyn/hdyn/sdyn/evolve/scatter_MPI.C \
		src/node/dyn/hdyn/sdyn/evolve/scatter_MPI_old.C \
		src/node/dyn/hdyn/sdyn/evolve/sigma.h \
		src/node/dyn/hdyn/sdyn/evolve/sigma_MPI.h \
		src/node/dyn/hdyn/sdyn/evolve/test_low_n_evolve.C \
		src/node/dyn/hdyn/sdyn3/evolve/all_sig.C \
		src/node/dyn/hdyn/sdyn3/evolve/bound3.C \
		src/node/dyn/hdyn/sdyn3/evolve/debug3.C \
		src/node/dyn/hdyn/sdyn3/evolve/d_ecc.C \
		src/node/dyn/hdyn/sdyn3/evolve/d_ecc_test.C \
		src/node/dyn/hdyn/sdyn3/evolve/delta_a.C \
		src/node/dyn/hdyn/sdyn3/evolve/delta_e.C \
		src/node/dyn/hdyn/sdyn3/evolve/flyby3.C \
		src/node/dyn/hdyn/sdyn3/evolve/hier3.C \
		src/node/dyn/hdyn/sdyn3/evolve/interfaces \
		src/node/dyn/hdyn/sdyn3/evolve/multiscatter3_pvm.C \
		src/node/dyn/hdyn/sdyn3/evolve/new_stats.C \
		src/node/dyn/hdyn/sdyn3/evolve/out3.C \
		src/node/dyn/hdyn/sdyn3/evolve/perturb.C \
		src/node/dyn/hdyn/sdyn3/evolve/pscat.C \
		src/node/dyn/hdyn/sdyn3/evolve/pvm_pack.C \
		src/node/dyn/hdyn/sdyn3/evolve/pvm_scatt.h \
		src/node/dyn/hdyn/sdyn3/evolve/rate3.C \
		src/node/dyn/hdyn/sdyn3/evolve/rate_stats.C \
		src/node/dyn/hdyn/sdyn3/evolve/scatter3_slave.C \
		src/node/dyn/hdyn/sdyn3/evolve/scatt_stats_1.C \
		src/node/dyn/hdyn/sdyn3/evolve/scatt_stats_2a.C \
		src/node/dyn/hdyn/sdyn3/evolve/scatt_stats_2.C \
		src/node/dyn/hdyn/sdyn3/evolve/scatt_stats.C \
		src/node/dyn/hdyn/sdyn3/evolve/sigma3.C \
		src/node/dyn/hdyn/sdyn3/evolve/SPECS \
		src/node/dyn/hdyn/sdyn3/evolve/template.C \
		src/node/dyn/hdyn/sdyn3/evolve/test_2.C \
		src/node/dyn/hdyn/sdyn3/evolve/test3.C \
		src/node/dyn/hdyn/sdyn3/xutil \
		src/node/dyn/hdyn/util/hdyn_sync.C \
		src/node/dyn/hdyn/util/kira_debug.h \
		src/node/dyn/hdyn/util/runf.C \
		src/node/dyn/hdyn/util/sync_system.C \
		src/node/dyn/hdyn/util/test.C \
		src/node/dyn/hdyn/util/wunf.C \
		src/node/dyn/hdyn/xutil/Demo.csh \
		src/node/dyn/hdyn/xutil/dummy.C \
		src/node/dyn/hdyn/xutil/dyn_util.h \
		src/node/dyn/hdyn/xutil/old_snap_to_image.C \
		src/node/dyn/hdyn/xutil/test.c \
		src/node/dyn/hdyn/xutil/write_png.C \
		src/node/dyn/hdyn/xutil/write_png.h \
		src/node/dyn/hdyn/xutil/xstarplot.h \
		src/node/dyn/init/addstoa.C \
		src/node/dyn/init/king.f \
		src/node/dyn/init/readgadget.C \
		src/node/dyn/init/split_particle.C \
		src/node/dyn/init/Testfile \
		src/node/dyn/kepler/save_kepler.C \
		src/node/dyn/tdyn/attach_new_node.C \
		src/node/dyn/tdyn/dyn_util.h \
		src/node/dyn/tdyn/inline.h \
		src/node/dyn/tdyn/Parti.todo \
		src/node/dyn/tdyn/Parti.txt \
		src/node/dyn/tdyn/pdyn_alloc.C \
		src/node/dyn/tdyn/print_local2.C \
		src/node/dyn/tdyn/print_local.C \
		src/node/dyn/tdyn/profile_new_interp.txt \
		src/node/dyn/tdyn/profile_old_interp.txt \
		src/node/dyn/tdyn/runf.C \
		src/node/dyn/tdyn/sizeof.C \
		src/node/dyn/tdyn/TO_DO \
		src/node/dyn/tdyn/update_node.C \
		src/node/dyn/tdyn/xstarplot2.C \
		src/node/dyn/tdyn/xstarplot.h \
		src/node/dyn/util/get_profile.C \
		src/node/dyn/util/merge_snaps.C \
		src/node/dyn/util/print_radii.C \
		src/node/util/makemultiples.C \
		src/node/util/mass.C \
		src/node/util/rmnode.C \
		src/node/util/story_test.C \
		src/node/util/test.C \
		src/node/util/Testfile \
		src/packages/autoconf-2.59.tar.gz \
		src/packages/automake-1.8.5.tar.gz \
		src/packages/banner-1.3.2.tar.gz \
		src/packages/clean-autotools \
		src/packages/gifsicle-1.40.tar.gz \
		src/packages/install-autotools \
		src/packages/install-packages \
		src/packages/libtool-1.5.6.tar.gz \
		src/packages/Makefile \
		src/packages/package-list \
		src/packages/shhmsg-1.4.1.tar.gz \
		src/packages/shhopt-1.1.7.tar.gz \
		src/star/dstar/evolve/binev.C \
		src/star/dstar/evolve/collexp.C \
		src/star/dstar/evolve/dbinev.C \
		src/star/dstar/evolve/nomoto.C \
		src/star/dstar/evolve/nSeBa.C \
		src/star/dstar/evolve/SeBa.C.GN \
		src/star/dstar/evolve/SeBa.C.SPZ \
		src/star/dstar/evolve/SeBa_GN.C \
		src/star/dstar/evolve/tripev.C \
		src/star/dstar/sdyn3/SeBa_tt.C \
		src/star/dstar/sdyn3/SeBa_tt.h \
		src/star/dstar/starclass/seba_counters.C \
		src/star/dstar/util/example.in \
		src/star/dstar/util/example.log \
		src/star/dstar/util/example.ps \
		src/star/dstar/util/makedouble.C \
		src/star/dstar/util/roche.C \
		src/star/dstar/util/roche.h \
		src/star/dstar/util/roche.txt \
		src/star/rdc \
		src/star/sstar/evolve/core_mass.C \
		src/star/sstar/evolve/evolve_star.C \
		src/star/sstar/evolve/SeBa.data \
		src/star/sstar/starclass/general_support.C \
		src/star/sstar/starclass/SeBa.data \
		src/star/sstar/util \
		src/star/util/BminV.tab \
		src/star/util/BolCor.tab \
		src/star/util/conv_dyn_to_star.C \
		src/star/util/conv_star_to_dyn.C \
		src/star/util/dummy.C \
		src/star/util/Longint \
		src/star/util/OortAB.C \
		src/star/util/rcore2_bk.C \
		src/star/util/rcore.C \
		src/star/util/RminI.tab \
		src/star/util/sint.C \
		src/star/util/sint.h \
		src/star/util/test.C \
		src/star/util/test_sint.C \
		src/star/util/UminB.tab \
		src/star/util/Version_1 \
		src/star/util/Version_2 \
		src/star/util/Version_3 \
		src/star/util/VminR.tab \
		src/star/xutil \
		src/std/mycat.C \
		src/std/set_up.C \
		src/std/setup_starlab.C \
		src/std/sofar_story.C \
		src/std/stats.C \
		src/std/test.C \
		src/std/xtest.C \
		TOOLS

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
