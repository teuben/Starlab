# 
#       #=======================================================#     _\|/_
#      #  __  _____           ___                    ___       #       /|\
#     #  /      |      ^     |   \  |         ^     |   \     #           _\|/_
#    #   \__    |     / \    |___/  |        / \    |___/    #             /|\
#   #       \   |    /___\   |  \   |       /___\   |   \   #  _\|/_
#  #     ___/   |   /     \  |   \  |____  /     \  |___/  #    /|\
# #                                                       #             _\|/_
##=======================================================#               /|\
#
##
##  src/std/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC)
LIBFLAGS	= -L$(LIB) -lstd -lm
DEPEND_LIBS	= $(LIB)/libstd.a

EXECUTABLES =  gethist pgetopt randinter runtime_help precision \
	       story test_vector
OBJECT_FILES = randinter.o pgetopt.o gethist.o precision.o story.o \
		util_math.o util_io.o error.o time.o par_to_use.o \
		printbits.o runtime_help.o xreal.o
LIBRARY = libstd.a

EXTRA_CLEAN = $(LIB)/$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#randinter:	randinter.C
#	@make_cond_build $@ $(BIN)/$@ \
#		"$(LIBFLAGS) " $(STARLAB_CCC) $(CFLAGS) $(TOOLFLAGS)

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/std/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

