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
##  src/gfx/util/Makefile
##.............................................................................
##    version 1.1:	Feb 1996	Steve McMillan
##.............................................................................
##

##---------------------------------------------------------------------------

#### Site-specific, for now: should add to configure...
#### Set these ONLY if the png and zlib libraries are installed.
#### If the libraries don't exist, comment these lines out and the
#### PNG option in snap_to_image will default to GIF.

PNG_CC_OPTION = -DHAS_PNG
PNG_LIB_OPTION = -lpng

##---------------------------------------------------------------------------

CFLAGS		=  $(STARLAB_CCC_FLAGS) -I$(INC) $(PNG_CC_OPTION)
LIBFLAGS	= -L$(LIB) -llux -lstd  $(PNG_LIB_OPTION) -lm

DEPEND_LIBS	= $(LIB)/libstd.a 

LIBRARY = liblux.a
LIB_DEP = .$(LIBRARY)

OBJECT_FILES =	colormaps.o write_sun.o sun_header.o \
		write_gif.o encode_gif.o write_png.o

# make_header.o write_image.o

EXECUTABLES = make_image data

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................


#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/gfx/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

