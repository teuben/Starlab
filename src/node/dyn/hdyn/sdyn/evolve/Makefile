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
##  src/node/dyn/hdyn/sdyn/evolve/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................

#CFLAGS		= -O3 $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC)

# Starlab libraries to link with (LIB is STARLAB_PATH/lib):

#LIBFLAGS	= -L$(LIB) -L$(STARLAB_LIB) -lstd -lm -lsdyn3
#CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) -DUSE_MPI

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC)
LIBFLAGS	=  -L. -lsc4 -L$(STARLAB_LIB) -lsdyn -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm 

#LIBFLAGS	=  -L. -lsc4 -L$(STARLAB_LIB) -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm -L$(STARLAB_LIB) -lsdyn -ldyn -lmpi
#LIBFLAGS	=  -L. -lsc4 -L$(STARLAB_LIB) -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm -L$(STARLAB_LIB) -lsdyn -ldyn -lmpi

DEPEND_LIBS	= ./libsc4.a \
		  $(STARLAB_LIB)/libsdyn.a \
		  $(STARLAB_LIB)/libhdyn.a \
	          $(STARLAB_LIB)/lib_dyn_.a \
		  $(STARLAB_LIB)/libdyn.a \
		  $(STARLAB_LIB)/libnode.a \
	          $(STARLAB_LIB)/libstd.a 

#EXECUTABLES = low_n_evolve tree_evolve scatter sigma \
#		scatter_MPI sigma_MPI
EXECUTABLES = low_n_evolve tree_evolve scatter sigma 

OBJECT_FILES = low_n_evolve.o sdyn_ev.o tree_evolve.o sdyn_ev.o \
	       scatter_exp.o scatter_hist.o \
	       multiscatter.o tree_evolve.o sig_helper.o get_sigma.o \
	       scatter_merge.o extend_or_end.o \
	       scatter_MPI.o multiscatter_MPI.o MPI_sigma.o \
	       specific_function.o	

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

BIN = .
STARLAB_LIB = $(STARLAB_PATH)/lib
LIB = .

LIBRARY = libsc4.a
LIB_DEP = .$(LIBRARY)

#..............................................................................
#..............................................................................

##===================================================================//
##  +---------------+      _\|/_      +------------------------------\\
##  |  the end of:  |       /|\       |  src/node/dyn/hdyn/sdyn/evolve/Makefile
##  +---------------+                 +------------------------------//
##======================= STARLAB ===================================\\

