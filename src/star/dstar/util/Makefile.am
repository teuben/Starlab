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
##	src/star/dstar/util/Makefile.am
##..............................................................................
##	Oct 1994		Simon F. Portegies Zwart
##	Jul 2004		Ernest Mamikonyan
##	Sep 2004		Annotated by Steve McMillan
##..............................................................................
##

# Executables to build:

bin_PROGRAMS = semi_to_period

# Specific options and flags for each executable:

semi_to_period_SOURCES  = semi_to_period.C
semi_to_period_CPPFLAGS = -DTOOLBOX
semi_to_period_LDADD =	../libdstar.la \
			../../../node/dyn/hdyn/libhdyn.la

# Local library and sources (sub-libraries will be combined later):

noinst_LTLIBRARIES = libdstar-util.la
libdstar_util_la_SOURCES = dstar_to_node.C semi_to_period.C
