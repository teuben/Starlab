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
##  src/evolve/single/star/Makefile
##.............................................................................
##    version 2.1:	March 1994	Simon F. Portegies Zwart
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC) \
		  $(STARLAB_GFX_CFLAGS) $(STARLAB_EXTRAS)
LIBFLAGS	= -L$(LIB) -lsstar -lsdyn3 -lhdyn -ldyn -lnode \
	   	  -lstd -llux $(STARLAB_POSIX_LIB) \
	   	  -L$(STARLAB_X11LIB) -lX11 $(STARLAB_XLIB) -lm
DEPEND_LIBS	= $(LIB)/libstar.a $(LIB)/libsdyn3.a $(LIB)/libhdyn.a \
		  $(LIB)/libdyn.a $(LIB)/libnode.a $(LIB)/libnode.a

include $(STARLAB_PATH)/inc/Makefile.worker.inc

EXECUTABLES = xhrdplot

xhrdplot:	 xhrdplot.C gfx_util.o dyn_util.o hrd_util.o
	  $(STARLAB_CCC) $(CFLAGS) $(STARLAB_CCC_FLAGS) -o $(BIN)/$@ \
	  xhrdplot.C gfx_util.o dyn_util.o hrd_util.o $(LIBFLAGS)
	@touch xhrdplot


LIBRARY = 
LIB_DEP = 

#..............................................................................


#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/evolve/single/star/Makefile
##  +---------------+                     +------------------------------//
##========================= STAREV =====================================\\

