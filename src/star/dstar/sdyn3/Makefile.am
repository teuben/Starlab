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
##	src/star/dstar/sdyn3/Makefile.am
##..............................................................................
##	Oct 1994		Simon F. Portegies Zwart
##	Sep 2004		Steve McMillan
##	Sep 2004		Annotated by Steve McMillan
##..............................................................................
##

# Executables to build:

bin_PROGRAMS = SeBa_scatter3

# Specific options and flags for each executable:

SeBa_scatter3_SOURCES  = SeBa_scatter3.C
SeBa_scatter3_CPPFLAGS = -DTOOLBOX
SeBa_scatter3_LDADD =	../libdstar.la \
			../../../node/dyn/hdyn/sdyn3/libsdyn3.la \
			../../../node/dyn/hdyn/libhdyn.la 

# Local library and sources (sub-libraries will be combined later):

noinst_LTLIBRARIES = libdstar-sdyn3.la
libdstar_sdyn3_la_SOURCES = SeBa_scatter3.C SeBa_tt.C
