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
##	src/node/dyn/hdyn/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jul 2004		Ernest Mamikonyan
##..............................................................................
##

SDYNS =
if USE_SDYN3
  SDYNS += sdyn3
endif
if USE_SDYN
  SDYNS += sdyn
endif

SUBDIRS = util evolve xutil stubs $(SDYNS)

pkglib_LTLIBRARIES = libhdyn.la

libhdyn_la_LDFLAGS = ../_dyn_/lib_dyn_.la
libhdyn_la_LIBADD = util/libhdyn-1.la evolve/libhdyn-2.la
libhdyn_la_SOURCES =

.PHONY: libs clibs
libs: clibs libhdyn.la
clibs:
	$(MAKE) $(AM_MAKEFLAGS) -C util libhdyn-1.la
	$(MAKE) $(AM_MAKEFLAGS) -C evolve libhdyn-2.la
