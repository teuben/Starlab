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
##  src/star/init/Makefile
##.............................................................................
##    version 1.1:	Dec 1993	Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC)
LIBFLAGS	= -L$(LIB) -lsstar -ldyn -lnode -lstd -lm
DEPEND_LIBS	= $(LIB)/libsstar.a $(LIB)/libdyn.a \
		  $(LIB)/libnode.a $(LIB)/libstd.a

# N.B. the -ldstar here is to satisfy references to make_state(double_star*)
#      and put_state(double_star*) forced in star_state...

EXECUTABLES = sstar_stats
OBJECT_FILES = sstar_to_dyn.o sstar_state.o stardyn_story.o addstar_dyn.o \
	       sstar_stats.o proto_star_dyn.o

LIBRARY = libsstar.a
LIB_DEP = .$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/star/init/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

