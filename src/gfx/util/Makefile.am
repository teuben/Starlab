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

noinst_LTLIBRARIES = liblux-3.la
bin_PROGRAMS = data make_image

liblux_3_la_SOURCES = colormaps.C write_sun.C sun_header.C write_gif.C \
    encode_gif.C write_png.C

data_LDADD = ../liblux.la
make_image_LDADD = ../liblux.la

data_SOURCES = data.C
make_image_SOURCES = make_image.C
