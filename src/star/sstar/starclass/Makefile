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

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC) 
LIBFLAGS	= -L$(LIB) -lsstar -lstd -lm
DEPEND_LIBS	= $(LIB)/libsstar.a $(LIB)/libstd.a 
#	          $(LIB)/libdstar.a

EXECUTABLES = 
OBJECT_FILES = star.o single_star.o \
	       static_star.o SPZDCH_star.o proto_star.o \
	       brown_dwarf.o main_sequence.o hertzsprung_gap.o \
	       sub_giant.o horizontal_branch.o \
	       super_giant.o thorne_zytkow.o helium_star.o helium_giant.o \
	       hyper_giant.o white_dwarf.o neutron_star.o \
	       black_hole.o disintegrated.o \
	       star_support.o star_state.o stdfunc.o constants.o

LIBRARY = 
LIB_DEP = 

LIBRARY = libsstar.a
LIB_DEP = .$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/evolve/single/star/Makefile
##  +---------------+                     +------------------------------//
##========================= STAREV =====================================\\

