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
##  src/node/dyn/xutil/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) $(STARLAB_GFX_CFLAGS) \
		  $(STARLAB_EXTRAS)
LIBFLAGS	= -L$(LIB) -lsdyn3 -l_dyn_ -ldyn -lnode \
	   	  -lstd -llux $(STARLAB_POSIX_LIB) \
	   	  -L$(STARLAB_X11LIB) -lX11 $(STARLAB_XLIB) -lm

EXECUTABLES = xstarplot3

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................

xstarplot3:	xstarplot3.C gfx_util.o dyn_util.o
	$(STARLAB_CCC) $(CFLAGS) $(STARLAB_CCC_FLAGS) -o $(BIN)/$@ \
		xstarplot3.C gfx_util.o dyn_util.o $(LIBFLAGS) \
		-D_SRC_=\"$(PWD)/$@.C\"
	@touch xstarplot

new_xstarplot3:	new_xstarplot3.C gfx_util.o dyn_util.o
	$(STARLAB_CCC) $(CFLAGS) $(STARLAB_CCC_FLAGS) -o $(BIN)/$@ \
		new_xstarplot3.C gfx_util.o dyn_util.o $(LIBFLAGS) \
		-D_SRC_=\"$(PWD)/$@.C\"
	@touch new_xstarplot

#..............................................................................


#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/dyn/xutil/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
