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

# NOTE that this xstarplot will supercede the one in dyn/xutil if the
# hdyn option is selected.



#### Site-specific, for now: must add to configure...
#### Set these ONLY if the png and zlib libraries are installed.

#PNG_CC_OPTION = -DHAS_PNG
#PNG_LIB_OPTION = -lpng



CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) $(STARLAB_GFX_CFLAGS) \
		  $(STARLAB_EXTRAS) $(PNG_CC_OPTION)
LIBFLAGS	= -L$(LIB) -lhdyn -l_dyn_ -ldyn -lnode \
		  -lstd -llux $(STARLAB_POSIX_LIB) \
		  $(STARLAB_X11LIB) $(STARLAB_XLIB) $(PNG_LIB_OPTION) -lm

EXECUTABLES = xstarplot snap_to_image write_png
OBJECT_FILES = write_png.o

LIBRARY = libhdyn.a

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

# Special target (extra objects not in library):

xstarplot:	xstarplot.C gfx_util.o dyn_util.o
	$(STARLAB_CCC) $(CFLAGS) $(STARLAB_CCC_FLAGS) -o $(BIN)/$@ \
		xstarplot.C gfx_util.o dyn_util.o $(LIBFLAGS) \
		-D_SRC_=\"$(PWD)/$@.C\" \
		-D_COMPILE_DATE_=\"`date +%a_%b_%d_%Y,_%T_%Z`\"
	@touch xstarplot

#..............................................................................


#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/dyn/xutil/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
