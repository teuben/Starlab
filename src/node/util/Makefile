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
##  src/node/util/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) $(STARLAB_HAS_GZIP)
LIBFLAGS	= -L$(LIB) -lnode -lstd -lm
DEPEND_LIBS	= $(LIB)/libnode.a $(LIB)/libstd.a

EXECUTABLES = node_io display_subtrees rmq snapprune mk_single_node \
		add_daughter_node pretty_print_tree print_normal \
		mknode mkmass mksecondary mass_dist m_bar \
	        setmass
OBJECT_FILES = node_io.o hydrobase_io.o starbase_io.o tree_io.o node_tt.o \
		rmtree.o print_normal.o starbase_init.o \
		mknode.o mkmass.o

LIBRARY = libnode.a

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

