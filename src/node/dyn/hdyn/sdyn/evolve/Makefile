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

# Starlab libraries to link with (LIB is STARLAB_PATH/lib):

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) 
LIBFLAGS	=  -L. -lscat -L$(STARLAB_LIB) -lsdyn -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm 

#LIBFLAGS	=  -L. -lsc4 -L$(STARLAB_LIB) -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm -L$(STARLAB_LIB) -lsdyn -ldyn -lmpi

DEPEND_LIBS	= ./libsc4.a \
		  $(STARLAB_LIB)/libsdyn.a \
		  $(STARLAB_LIB)/libhdyn.a \
	          $(STARLAB_LIB)/lib_dyn_.a \
		  $(STARLAB_LIB)/libdyn.a \
		  $(STARLAB_LIB)/libnode.a \
	          $(STARLAB_LIB)/libstd.a 

EXECUTABLES = low_n_evolve tree_evolve scatter sigma 

OBJECT_FILES = low_n_evolve.o sdyn_ev.o tree_evolve.o sdyn_ev.o \
	       scatter_exp.o scatter_hist.o \
	       multiscatter.o tree_evolve.o sig_helper.o \
	       scatter_merge.o extend_or_end.o \
	       scatter_MPI.o multiscatter_MPI.o MPI_sigma.o \
	       specific_function.o scatter.o

BIN = .
STARLAB_LIB = $(STARLAB_PATH)/lib
#LIB = .

LIBRARY = libscat.a
LIB_DEP = .$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

# Extra targets specific to MPI:

CFLAGS_MPI	= $(CFLAGS) -DHAS_MPI
LIBFLAGS_MPI	=  -L. -L$(STARLAB_LIB) -lsdyn -lhdyn -l_dyn_ -ldyn -lnode -lstd -lm 

.SUFFIXES:	$(SUFFIXES) .mpi .mpio

MPI_EXECUTABLES = scatter_MPI.mpi 

MPI_OBJECTS = low_n_evolve.mpio sdyn_ev.mpio \
	       scatter_exp.mpio scatter_hist.mpio \
	       tree_evolve.mpio \
	       scatter_merge.mpio extend_or_end.mpio \
	       scatter_MPI.mpio \
	       specific_function.mpio	\
	       MPI_sigma.mpio sig_helper.mpio multiscatter_MPI.mpio  

mpi_bin:	$(MPI_EXECUTABLES)
mpi_sure:	mpi_clean mpi_bin
mpi_lib:	$(MPI_OBJECTS)


.C.mpi:
	@if [ "X$(MPI_ROOT)" != "X" ]; then \
	    for o in $(MPI_OBJECTS) ; do \
		if [ ! -f $$o ]; then make $$o ; fi ; \
	    done ; \
	    echo $(STARLAB_MPICCC) $(CFLAGS_MPI) -I $(MPI_ROOT) $(TOOLFLAGS) \
		 $< $(MPI_OBJECTS) -o $(BIN)/$@ $(LIBFLAGS_MPI) \
		 $(STARLAB_MPI_LIBS) -lm ; \
	    $(STARLAB_MPICCC) $(CFLAGS_MPI) \
		-I $(MPI_ROOT) $(TOOLFLAGS) \
		 $< $(MPI_OBJECTS) -o $(BIN)/$@ $(LIBFLAGS_MPI) \
		 -D_SRC_=\"$(PWD)/$<\" \
		 -D_COMPILE_DATE_=\"Some_time_some_date\" \
		 $(STARLAB_MPI_LIBS) -lm ; \
	    touch $@ ; \
	fi

.C.mpio:
	$(STARLAB_MPICCC) $(CFLAGS) -DHAS_MPI -I$(MPI_ROOT) -c $< -o $@

mpi_clean:
	@/bin/rm -f *.mpi *.mpio

#..............................................................................
#..............................................................................

##===================================================================//
##  +---------------+      _\|/_      +------------------------------\\
##  |  the end of:  |       /|\       |  src/node/dyn/hdyn/sdyn/evolve/Makefile
##  +---------------+                 +------------------------------//
##======================= STARLAB ===================================\\

