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
##  src/node/dyn/util/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		=  $(STARLAB_CCC_FLAGS) -I$(INC) -g3

#LIBFLAGS	= -L$(LIB) -ldyn -lnode -lstd -lm

LIBFLAGS	= -L$(LIB) -ldyn -lsstar -lnode -lstd -lm

# (The -lsstar here is needed only to compile sys_stats.
#  Note that the library is not listed in DEPEND_LIBS.)

DEPEND_LIBS	= $(LIB)/libdyn.a $(LIB)/libnode.a $(LIB)/libstd.a 

EXECUTABLES =	compute_com compute_density compute_max_cod \
		compute_mean_cod dumbp dyn_io dyndiff dyndist energy \
		extract_snap flatten freeze lagrad lagradplot list_snap \
		merge_low molecules reflect_velocities scale \
		plot_stars starplot snap2speck \
		sys_stats to_cod to_com jiggle \
                lagradplotsmooth

OBJECT_FILES =	compute_com.o compute_density.o compute_max_cod.o \
		compute_mean_cod.o dyn_di.o dyn_io.o dyn_stats.o dyn_tt.o \
		flatten.o lagrad.o merge_low.o scale.o sys_stats.o to_com.o \
		plot_stars.o

LIBRARY = libdyn.a

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/dyn/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

