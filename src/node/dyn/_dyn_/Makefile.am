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
##	src/node/dyn/_dyn_/Makefile.am
##..............................................................................
##	Jun 1999		Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

pkglib_LTLIBRARIES = lib_dyn_.la
bin_PROGRAMS = flat_hermite _dyn_io no_xreal set_radius

lib_dyn__la_LIBADD = ../libdyn.la
lib_dyn__la_SOURCES = _dyn_ev.C _dyn_io.C _dyn_slow.C

flat_hermite_CPPFLAGS = -DTOOLBOX
_dyn_io_CPPFLAGS = -DTOOLBOX
no_xreal_CPPFLAGS = -DTOOLBOX
set_radius_CPPFLAGS = -DTOOLBOX

flat_hermite_LDADD = lib_dyn_.la
_dyn_io_LDADD = lib_dyn_.la
no_xreal_LDADD = lib_dyn_.la
set_radius_LDADD = lib_dyn_.la

flat_hermite_SOURCES = flat_hermite.C
_dyn_io_SOURCES = _dyn_io.C
no_xreal_SOURCES = no_xreal.C
set_radius_SOURCES = set_radius.C

PHONY: libs
libs: lib_dyn_.la
