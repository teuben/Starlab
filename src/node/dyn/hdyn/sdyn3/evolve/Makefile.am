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
##	src/node/dyn/hdyn/sdyn3/evolve/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

noinst_LTLIBRARIES = libsdyn3-2.la

libsdyn3_2_la_SOURCES = extend_or_end.C low_n3_evol.C scat_helper3.C \
    scat_init.C scatter3.C sdyn3_ev.C sig_helper3.C tree3_evolve.C \
    multiscatter3.C

bin_PROGRAMS = bound3 low_n3_evol rate3 scatter3 sigma3 tree3_evolve # pvm_bin

bound3_CPPFLAGS = -DTOOLBOX
low_n3_evol_CPPFLAGS = -DTOOLBOX
rate3_CPPFLAGS = -DTOOLBOX
scatter3_CPPFLAGS = -DTOOLBOX
sigma3_CPPFLAGS = -DTOOLBOX
tree3_evolve_CPPFLAGS = -DTOOLBOX

bound3_LDADD = ../libsdyn3.la
low_n3_evol_LDADD = ../libsdyn3.la
rate3_LDADD = ../libsdyn3.la
scatter3_LDADD = ../libsdyn3.la
sigma3_LDADD = ../libsdyn3.la
tree3_evolve_LDADD = ../libsdyn3.la

bound3_SOURCES = bound3.C
low_n3_evol_SOURCES = low_n3_evol.C
rate3_SOURCES = rate3.C
scatter3_SOURCES = scatter3.C
sigma3_SOURCES = sigma3.C
tree3_evolve_SOURCES = tree3_evolve.C
