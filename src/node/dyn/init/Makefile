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
##  src/node/dyn/init/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC)
LIBFLAGS	= -L$(LIB) -ldyn -lnode -lstd -lm
DEPEND_LIBS	= $(LIB)/libstd.a $(LIB)/libnode.a $(LIB)/libdyn.a

EXECUTABLES	= makecube makedisk makeplummer makepyth makesphere mkwrite \
		  addplummer mkbinary split_particles makeking \
	 	  mkheavystar mk_aniso_king readp readstoa \
	          mkplanetary mkblack_hole grow_black_hole

OBJECT_FILES	= mksphere.o
LIBRARY		= libdyn.a

aking.o:	aking.f
	$(STARLAB_FORTRAN_COMPILE) $(STARLAB_FORTRAN_FLAGS) -g3 -c aking.f

mk_aniso_king:	mk_aniso_king.C aking.o
	$(STARLAB_CCC) $(CFLAGS) \
	    $(TOOLFLAGS) $(STARLAB_FORTRAN_CFLAGS) \
	    $@.C aking.o -o $(BIN)/$@ \
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
##  |  the end of:  |         /|\         |  src/node/dyn/init/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
