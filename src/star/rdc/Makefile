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
##  src/star/rdc/Makefile
##.............................................................................
##    version 1.1:	Dec 1993	Piet Hut & Steve McMillan
##.............................................................................
##


CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC) -I$(STARLAB_CFITSIOINC)
LIBFLAGS	= -L$(LIB) -L$(STARLAB_CFITSIOLIB) -lrdc -ldstar -lsstar -lhdyn -l_dyn_ -ldyn -lnode -lstd  -lcfitsio -lm
DEPEND_LIBS	= $(LIB)/libsstar.a  $(LIB)/libdstar.a $(LIB)/librdc.a \
	          $(LIB)/libdyn.a $(LIB)/libnode.a   $(LIB)/libstd.a \
		  -L$(STARLAB_CFITSIOLIB)/libcfitsio.a

EXECUTABLES = rdc_SeBa rs_findtype rs_snapshot \
	      rs_population rs_scenarios rs_birthrate \
	      rs_countsn rs_reorder
#		makepovfile

# make_ccd	# removed by Steve, 8/01 -- no fitsio.h

OBJECT_FILES =  
#		SeBa_hist.o makepovfile.o

#OBJECT_FILES = sint.o
# N.B. sint.o suppressed until bugs in sint.C are repaired...

LIBRARY = librdc.a
LIB_DEP = .$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/star/rdc/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

