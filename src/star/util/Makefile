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

CFLAGS          = $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC)
LIBFLAGS        = -L$(LIB) -lsstar -ldstar -ldyn -lnode \
	          -lstd -lm

DEPEND_LIBS     = $(LIB)/libsint.a \
	          $(LIB)/libsstar.a  $(LIB)/libdstar.a  $(LIB)/libhdyn.a \
	          $(LIB)/libdyn.a $(LIB)/libnode.a $(LIB)/libstd.a

EXECUTABLES  = red_stellar_system rcore2 coagulation rearrange_mf \
	       red_star_cluster encountervelocity imf_evolve proj_lagr \
               lmt_to_ubvri mkdyn_from_input
OBJECT_FILES = proj_lagr.o kingfit.o lmt_to_ubvri.o

#OBJECT_FILES = sint.o
# N.B. sint.o suppressed until bugs in sint.C are repaired...

LIBRARY = libsstar.a
LIB_DEP = .$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/star/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

