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
##	src/star/xutil/Makefile.am
##..............................................................................
##	Mar 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

bin_PROGRAMS = xhrdplot

xhrdplot_CPPFLAGS = -DTOOLBOX
xhrdplot_LDADD = ../sstar/libsstar.la \
		 ../../node/dyn/hdyn/libhdyn.la \
		 ../../gfx/libgfx.la
xhrdplot_SOURCES = xhrdplot.C gfx_util.C dyn_util.C hrd_util.C
