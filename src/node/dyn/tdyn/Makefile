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
##  src/node/dyn/tdyn/Makefile
##.............................................................................
##    Sep 2000					Steve McMillan
#..............................................................................

OBJECT_FILES = tdyn_io.o worldline.o

LIBRARY = libtdyn.a

EXECUTABLES = read_tdyn sizeof tdyn_io test_tdyn worldline worldlines \
	      xstarplot2

NON_CCC_EXECUTABLES = 

DEPEND_LIBS	= $(LIB)/libtdyn.a $(LIB)/libdyn.a $(LIB)/libnode.a \
		  $(LIB)/libstd.a

CFLAGS    =	$(STARLAB_CCC_FLAGS) -I$(INC) -O0 $(STARLAB_GFX_CFLAGS)

CPPFLAGS = 

LIBFLAGS  =	-L$(LIB) -ltdyn -ldyn -lnode \
		  -lstd -llux $(STARLAB_POSIX_LIB) \
		  -L$(STARLAB_X11LIB) -lX11 $(STARLAB_XLIB) -lm

EXEFLAGS = 

# Special target (extra objects not in library):

xstarplot2:	xstarplot2.C gfx_util.o dyn_util.o
	$(STARLAB_CCC) $(CFLAGS) $(STARLAB_CCC_FLAGS) -o $(BIN)/$@ \
		xstarplot2.C gfx_util.o dyn_util.o $(LIBFLAGS) \
		-D_SRC_=\"$(PWD)/$@.C\" \
		-D_COMPILE_DATE_=\"`date +%a_%b_%d_%Y,_%T_%Z`\"
	@touch xstarplot2

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/dyn/_dyn_/evolve/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
