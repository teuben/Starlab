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
##	src/node/dyn/util/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libdyn-1.la
bin_PROGRAMS = compute_com compute_mcom compute_density compute_max_cod \
    compute_mean_cod dumbp dyn_io dyndiff dyndist energy extract_snap flatten \
    freeze freezebin lagrad lagradplot list_snap merge_low molecules \
    reflect_velocities scale plot_stars starplot snap2speck sys_stats set_com \
    to_com jiggle lagradplotsmooth get_std_center add_plummer add_power_law \
    add_tidal set_ignore_internal add_daughter_node getsnap print_time \
    radial_density radial_numdensity radial_vdisp radial_profile makemass \
    makesecondary

# this is for stpcpy(3) on the Pentium machines in Tokyo
libdyn_1_la_CPPFLAGS = -D_GNU_SOURCE
libdyn_1_la_SOURCES = compute_com.C compute_density.C compute_max_cod.C \
    compute_mcom.C compute_mean_cod.C dyn_di.C dyn_external.C dyn_io.C \
    dyn_stats.C dyn_story.C dyn_tt.C flatten.C lagrad.C merge_low.C scale.C \
    sys_stats.C to_com.C set_com.C plot_stars.C get_std_center.C refine_mass.C \
    refine_mass2.C add_plummer.C add_power_law.C radial_density.C \
    radial_vdisp.C radial_profile.C wrappers.C dyn_lexer.ll dyn_parser.yy \
    radial_numdensity.C

compute_com_CPPFLAGS = -DTOOLBOX
compute_mcom_CPPFLAGS = -DTOOLBOX
compute_density_CPPFLAGS = -DTOOLBOX
compute_max_cod_CPPFLAGS = -DTOOLBOX
compute_mean_cod_CPPFLAGS = -DTOOLBOX
dumbp_CPPFLAGS = -DTOOLBOX
dyn_io_CPPFLAGS = -DTOOLBOX
dyndiff_CPPFLAGS = -DTOOLBOX
dyndist_CPPFLAGS = -DTOOLBOX
energy_CPPFLAGS = -DTOOLBOX
extract_snap_CPPFLAGS = -DTOOLBOX
flatten_CPPFLAGS = -DTOOLBOX
freeze_CPPFLAGS = -DTOOLBOX
freezebin_CPPFLAGS = -DTOOLBOX
lagrad_CPPFLAGS = -DTOOLBOX
lagradplot_CPPFLAGS = -DTOOLBOX
list_snap_CPPFLAGS = -DTOOLBOX
merge_low_CPPFLAGS = -DTOOLBOX
molecules_CPPFLAGS = -DTOOLBOX
reflect_velocities_CPPFLAGS = -DTOOLBOX
scale_CPPFLAGS = -DTOOLBOX
plot_stars_CPPFLAGS = -DTOOLBOX
starplot_CPPFLAGS = -DTOOLBOX
snap2speck_CPPFLAGS = -DTOOLBOX
sys_stats_CPPFLAGS = -DTOOLBOX
set_com_CPPFLAGS = -DTOOLBOX
to_com_CPPFLAGS = -DTOOLBOX
jiggle_CPPFLAGS = -DTOOLBOX
lagradplotsmooth_CPPFLAGS = -DTOOLBOX
get_std_center_CPPFLAGS = -DTOOLBOX
add_plummer_CPPFLAGS = -DTOOLBOX
add_power_law_CPPFLAGS = -DTOOLBOX
add_tidal_CPPFLAGS = -DTOOLBOX
set_ignore_internal_CPPFLAGS = -DTOOLBOX
add_daughter_node_CPPFLAGS = -DTOOLBOX
getsnap_CPPFLAGS = -DTOOLBOX
print_time_CPPFLAGS = -DTOOLBOX
radial_density_CPPFLAGS = -DTOOLBOX
radial_numdensity_CPPFLAGS = -DTOOLBOX
radial_vdisp_CPPFLAGS = -DTOOLBOX
radial_profile_CPPFLAGS = -DTOOLBOX
makemass_CPPFLAGS = -DTOOLBOX
makesecondary_CPPFLAGS = -DTOOLBOX

compute_com_LDADD = ../libdyn.la
compute_mcom_LDADD = ../libdyn.la
compute_density_LDADD = ../libdyn.la
compute_max_cod_LDADD = ../libdyn.la
compute_mean_cod_LDADD = ../libdyn.la
dumbp_LDADD = ../libdyn.la
dyn_io_LDADD = ../libdyn.la
dyndiff_LDADD = ../libdyn.la
dyndist_LDADD = ../libdyn.la
energy_LDADD = ../libdyn.la
extract_snap_LDADD = ../libdyn.la
flatten_LDADD = ../libdyn.la
freeze_LDADD = ../libdyn.la
freezebin_LDADD = ../libdyn.la
lagrad_LDADD = ../libdyn.la
lagradplot_LDADD = ../libdyn.la
list_snap_LDADD = ../libdyn.la
merge_low_LDADD = ../libdyn.la
molecules_LDADD = ../libdyn.la
reflect_velocities_LDADD = ../libdyn.la
scale_LDADD = ../libdyn.la
plot_stars_LDADD = ../libdyn.la
starplot_LDADD = ../libdyn.la
snap2speck_LDADD = ../libdyn.la
sys_stats_LDADD = ../libdyn.la ../../../star/sstar/libsstar.la
set_com_LDADD = ../libdyn.la
to_com_LDADD = ../libdyn.la
jiggle_LDADD = ../libdyn.la
lagradplotsmooth_LDADD = ../libdyn.la
get_std_center_LDADD = ../libdyn.la
add_plummer_LDADD = ../libdyn.la
add_power_law_LDADD = ../libdyn.la
add_tidal_LDADD = ../libdyn.la
set_ignore_internal_LDADD = ../libdyn.la
add_daughter_node_LDADD = ../libdyn.la
getsnap_LDADD = ../libdyn.la
print_time_LDADD = ../libdyn.la
radial_density_LDADD = ../libdyn.la
radial_numdensity_LDADD = ../libdyn.la
radial_vdisp_LDADD = ../libdyn.la
radial_profile_LDADD = ../libdyn.la
makemass_LDADD = ../libdyn.la
makesecondary_LDADD = ../libdyn.la

compute_com_SOURCES = compute_com.C
compute_mcom_SOURCES = compute_mcom.C
compute_density_SOURCES = compute_density.C
compute_max_cod_SOURCES = compute_max_cod.C
compute_mean_cod_SOURCES = compute_mean_cod.C
dumbp_SOURCES = dumbp.C
dyn_io_SOURCES = dyn_io.C
dyndiff_SOURCES = dyndiff.C
dyndist_SOURCES = dyndist.C
energy_SOURCES = energy.C
extract_snap_SOURCES = extract_snap.C
flatten_SOURCES = flatten.C
freeze_SOURCES = freeze.C
freezebin_SOURCES = freezebin.C
lagrad_SOURCES = lagrad.C
lagradplot_SOURCES = lagradplot.C
list_snap_SOURCES = list_snap.C
merge_low_SOURCES = merge_low.C
molecules_SOURCES = molecules.C
reflect_velocities_SOURCES = reflect_velocities.C
scale_SOURCES = scale.C
plot_stars_SOURCES = plot_stars.C
starplot_SOURCES = starplot.C
snap2speck_SOURCES = snap2speck.C
sys_stats_SOURCES = sys_stats.C
set_com_SOURCES = set_com.C
to_com_SOURCES = to_com.C
jiggle_SOURCES = jiggle.C
lagradplotsmooth_SOURCES = lagradplotsmooth.C
get_std_center_SOURCES = get_std_center.C
add_plummer_SOURCES = add_plummer.C
add_power_law_SOURCES = add_power_law.C
add_tidal_SOURCES = add_tidal.C
set_ignore_internal_SOURCES = set_ignore_internal.C
add_daughter_node_SOURCES = add_daughter_node.C
getsnap_SOURCES = getsnap.C
print_time_SOURCES = print_time.C
radial_density_SOURCES = radial_density.C
radial_numdensity_SOURCES = radial_numdensity.C
radial_vdisp_SOURCES = radial_vdisp.C
radial_profile_SOURCES = radial_profile.C
makemass_SOURCES = makemass.C
makesecondary_SOURCES = makesecondary.C

BUILT_SOURCES = dyn_parser.h
AM_YFLAGS = -d
