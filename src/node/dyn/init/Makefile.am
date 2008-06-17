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
##	src/node/dyn/init/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libdyn-3.la
bin_PROGRAMS = makecube makedisk makeplummer makepyth makesphere makewrite \
    apply_plummer makebinary split_particles makeking makeheavystar \
    readp readp2 \
    readstoa makeplanetary makeblack_hole grow_black_hole makepowerlaw \
    cutrandomsample

if HAVE_F77
  bin_PROGRAMS += make_aniso_king
  make_aniso_king_CPPFLAGS = -DTOOLBOX
  make_aniso_king_LDADD = ../libdyn.la $(FLIBS)
  make_aniso_king_SOURCES = make_aniso_king.C aking.f
endif

libdyn_3_la_SOURCES = makesphere.C

makecube_CPPFLAGS = -DTOOLBOX
makedisk_CPPFLAGS = -DTOOLBOX
makeplummer_CPPFLAGS = -DTOOLBOX
makepyth_CPPFLAGS = -DTOOLBOX
makesphere_CPPFLAGS = -DTOOLBOX
makewrite_CPPFLAGS = -DTOOLBOX
apply_plummer_CPPFLAGS = -DTOOLBOX
makebinary_CPPFLAGS = -DTOOLBOX
split_particles_CPPFLAGS = -DTOOLBOX
makeking_CPPFLAGS = -DTOOLBOX
makeheavystar_CPPFLAGS = -DTOOLBOX
readp_CPPFLAGS = -DTOOLBOX
readp2_CPPFLAGS = -DTOOLBOX
readstoa_CPPFLAGS = -DTOOLBOX
makeplanetary_CPPFLAGS = -DTOOLBOX
makeblack_hole_CPPFLAGS = -DTOOLBOX
grow_black_hole_CPPFLAGS = -DTOOLBOX
makepowerlaw_CPPFLAGS = -DTOOLBOX
cutrandomsample_CPPFLAGS = -DTOOLBOX

makecube_LDADD = ../libdyn.la
makedisk_LDADD = ../libdyn.la
makeplummer_LDADD = ../libdyn.la
makepyth_LDADD = ../libdyn.la
makesphere_LDADD = ../libdyn.la
makewrite_LDADD = ../libdyn.la
apply_plummer_LDADD = ../libdyn.la
makebinary_LDADD = ../libdyn.la
split_particles_LDADD = ../libdyn.la
makeking_LDADD = ../libdyn.la
makeheavystar_LDADD = ../libdyn.la
readp_LDADD = ../libdyn.la
readp2_LDADD = ../libdyn.la
readstoa_LDADD = ../libdyn.la
makeplanetary_LDADD = ../libdyn.la
makeblack_hole_LDADD = ../libdyn.la
grow_black_hole_LDADD = ../libdyn.la
makepowerlaw_LDADD = ../libdyn.la
cutrandomsample_LDADD = ../libdyn.la

makecube_SOURCES = makecube.C
makedisk_SOURCES = makedisk.C
makeplummer_SOURCES = makeplummer.C
makepyth_SOURCES = makepyth.C
makesphere_SOURCES = makesphere.C
makewrite_SOURCES = makewrite.C
apply_plummer_SOURCES = apply_plummer.C
makebinary_SOURCES = makebinary.C
split_particles_SOURCES = split_particles.C
makeking_SOURCES = makeking.C
makeheavystar_SOURCES = makeheavystar.C
readp_SOURCES = readp.C
readp2_SOURCES = readp2.C
readstoa_SOURCES = readstoa.C
makeplanetary_SOURCES = makeplanetary.C
makeblack_hole_SOURCES = makeblack_hole.C
grow_black_hole_SOURCES = grow_black_hole.C
makepowerlaw_SOURCES = makepowerlaw.C
cutrandomsample_SOURCES = cutrandomsample.C
