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
##	src/Makefile.am
##..............................................................................
##	Dec 1994		Piet Hut & Steve McMillan
##	Jun 2004		Ernest Mamikonyan
##..............................................................................
##


pkglib_LTLIBRARIES = libstd.la
bin_PROGRAMS = gethist pgetopt randinter runtime_help precision story \
    test_vector xtest xreal

libstd_la_LDFLAGS = -lm
libstd_la_SOURCES = randinter.C pgetopt.C gethist.C precision.C story.C \
    util_math.C util_io.C error.C time.C par_to_use.C printbits.C \
    runtime_help.C util_unformatted_io.C xreal.C

gethist_CPPFLAGS = -DTOOLBOX
pgetopt_CPPFLAGS = -DTOOLBOX
randinter_CPPFLAGS = -DTOOLBOX
runtime_help_CPPFLAGS = -DTOOLBOX
precision_CPPFLAGS = -DTOOLBOX
story_CPPFLAGS = -DTOOLBOX
test_vector_CPPFLAGS = -DTOOLBOX
xtest_CPPFLAGS = -DTOOLBOX
xreal_CPPFLAGS = -DTOOLBOX

gethist_LDADD = libstd.la
pgetopt_LDADD = libstd.la
randinter_LDADD = libstd.la
runtime_help_LDADD = libstd.la
precision_LDADD = libstd.la
story_LDADD = libstd.la
test_vector_LDADD = libstd.la
xtest_LDADD = libstd.la
xreal_LDADD = libstd.la

gethist_SOURCES = gethist.C
pgetopt_SOURCES = pgetopt.C
randinter_SOURCES = randinter.C
runtime_help_SOURCES = runtime_help.C
precision_SOURCES = precision.C
story_SOURCES = story.C
test_vector_SOURCES = test_vector.C
xtest_SOURCES = xtest.C
xreal_SOURCES = xreal.C

.PHONY: libs
libs: libstd.la
