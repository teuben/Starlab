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
##	src/gfx/util/Makefile.am
##..............................................................................
##	Feb 1996		Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libgfx-1.la
bin_PROGRAMS = data make_image

libgfx_1_la_SOURCES = colormaps.C write_sun.C sun_header.C write_gif.C \
    encode_gif.C write_png.C

# Seems to be needed so some configures (e.g. on MacOS) correctly
# include the X libraries:

AM_LDFLAGS = $(X_LIBS)

data_LDADD = ../libgfx.la
make_image_LDADD = ../libgfx.la

data_SOURCES = data.C
make_image_SOURCES = make_image.C
