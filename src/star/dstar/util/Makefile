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
##  src/evolve/triple/util/Makefile
##.............................................................................
##    version 1.0:	Okt 1994	Simon F. Portegies Zwart
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC) \
		  $(STARLAB_PGPLOTINC) \
		  $(STARLAB_EXTRAS)

#	          -fhandle-exceptions

LIBFLAGS	= -L$(LIB) -lsstar -ldstar -lhdyn -lsstar -ldyn -lnode -lstd\
	   	  $(STARLAB_X11LIB) $(STARLAB_XLIB) \
                  $(PGPLOT_PGPLOTLIB) -lm 
DEPEND_LIBS	= $(LIB)/libsstar.a  $(LIB)/libdstar.a  $(LIB)/libhdyn.a \
		  $(LIB)/libdyn.a $(LIB)/libnode.a $(LIB)/libstd.a \
		  $(STARLAB_PGPLOTINC)/libpgplot.a 

EXECUTABLES  = semi_to_period	# roche # *don't* include roche routinely...

OBJECT_FILES = dstar_to_node.o semi_to_period.o


LIBRARY = 
LIB_DEP = 

LIBRARY = libdstar.a
#LIB_DEP = .$(LIBRARY)

roche:	roche.o
	$(STARLAB_FORTRAN_COMPILE) $(CFLAGS) \
	    $(TOOLFLAGS) $(STARLAB_FORTRAN_CFLAGS) \
	    $@.o -o $(BIN)/$@ \
	    -D_SRC_=\"$(PWD)/$@.C\" \
	    -D_COMPILE_DATE_=\"`date +%a_%b_%d_%Y,_%T_%Z`\" \
	    $(LIBFLAGS) \
	    $(STARLAB_FORTRAN_LIBPATH) $(STARLAB_FORTRAN_LIB) 
	@touch $@

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         | src/dyn/evolve/triple/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STAREV =====================================\\

