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
##	src/gfx/lux/Makefile.am
##..............................................................................
##	Feb 1994		Biao Lu & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = liblux-1.la
bin_PROGRAMS = plot_data simple ximage

liblux_1_la_SOURCES = win.c draw.c draw1.c color.c dialog.c mcd.c interface.c \
    termio.c utility.c simple.c

simple_CPPFLAGS = -DTEST

plot_data_LDADD = ../liblux.la
simple_LDADD = ../liblux.la
ximage_LDADD = ../liblux.la
