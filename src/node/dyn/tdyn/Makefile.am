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
##	src/node/dyn/tdyn/Makefile.am
##..............................................................................
##	Sep 2000		Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

pkglib_LTLIBRARIES = libtdyn.la
bin_PROGRAMS = print_id scanbinary scanbundle tdyn_io worldlines xstarplot22

libtdyn_la_LIBADD = ../libdyn.la
libtdyn_la_SOURCES = tdyn_io.C tdyn_tt.C worldline.C worldline2.C

print_id_CPPFLAGS = -DTOOLBOX
scanbinary_CPPFLAGS = -DTOOLBOX
scanbundle_CPPFLAGS = -DTOOLBOX
tdyn_io_CPPFLAGS = -DTOOLBOX
worldlines_CPPFLAGS = -DTOOLBOX
xstarplot22_CPPFLAGS = -DTOOLBOX

print_id_LDADD = libtdyn.la
scanbinary_LDADD = libtdyn.la
scanbundle_LDADD = libtdyn.la
tdyn_io_LDADD = libtdyn.la
worldlines_LDADD = libtdyn.la
xstarplot22_LDADD = libtdyn.la ../../../gfx/liblux.la

print_id_SOURCES = print_id.C
scanbinary_SOURCES = scanbinary.C
scanbundle_SOURCES = scanbundle.C
tdyn_io_SOURCES = tdyn_io.C
worldlines_SOURCES = worldlines.C
xstarplot22_SOURCES = xstarplot22.C gfx_util.C dyn_util.C

.PHONY: libs
libs: libtdyn.la
