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

#INC_GADGET	= "/home/spz/cpp/gadget"
#CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) -I$(INC_GADGET)

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC)
LIBFLAGS	= -L$(LIB) -ldyn -lsstar -ldyn -lnode -lstd -lm
DEPEND_LIBS	= $(LIB)/libstd.a $(LIB)/libnode.a $(LIB)/libdyn.a

# (The -lsstar here is needed to accommodate get_addstar...
#  Note the annoying repetition of libdyn.a.  The star library is not listed
#  in DEPEND_LIBS.)

EXECUTABLES	= makecube makedisk makeplummer makepyth makesphere makewrite \
		  apply_plummer makebinary split_particles makeking \
	 	  makeheavystar readp readstoa \
	          makeplanetary makeblack_hole grow_black_hole \
		  make_aniso_king makepowerlaw \
	          cutrandomsample 

#                 include Gadget read/write routines should include also the 
#		  libraries
#		  readgadget	

# (Steve 9/01: place mk_aniso_king last, as the build may break on new systems.)

OBJECT_FILES	= makesphere.o
LIBRARY		= libdyn.a

aking.o:	aking.f
	$(STARLAB_FORTRAN_COMPILE) $(STARLAB_FORTRAN_FLAGS) -g3 -c aking.f

make_aniso_king:	make_aniso_king.C aking.o
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
