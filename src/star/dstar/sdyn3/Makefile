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
##  src/node/dyn/hdyn/sdyn3/scatter3/Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##    Sep 1996					Steve McMillan: added PVM
##.............................................................................

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC_STAR) -I$(INC)
LIBFLAGS	= -L$(LIB) -ldstar -lsstar -lstrsct \
		  -lsdyn3 -l_dyn_ -ldyn -lnode \
		  -lstd -lm
DEPEND_LIBS	= $(LIB)/libdstar.a $(LIB)/libsstar.a \
		  $(LIB)/libsdyn3.a $(LIB)/lib_dyn_.a \
		  $(LIB)/libdyn.a $(LIB)/libnode.a \
		  $(lib)/libstrsct.a \
		  $(LIB)/libstd.a

# Note that most applications are NOT routinely built.

EXECUTABLES =	SeBa_scatter3

OBJECT_FILES =  SeBa_scatter.o SeBa_tt.o

LIBRARY = libstrsct.a

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................

# Extra targets specific to PVM:

.SUFFIXES:	$(SUFFIXES) .pvm .pvmo

PVM_EXECUTABLES = scatter3_slave.pvm sigma3.pvm rate3.pvm delta_e.pvm
PVM_OBJECTS =  multiscatter3_pvm.pvmo pvm_pack.pvmo

pvm_bin:	$(PVM_EXECUTABLES)

pvm_sure:	pvm_clean pvm_bin

.C.pvm:
	@if [ "X$(PVM_ROOT)" != "X" ]; then \
	    for o in $(PVM_OBJECTS) ; do \
		if [ ! -f $$o ]; then make $$o ; fi ; \
	    done ; \
	    echo $(STARLAB_CCC) $(CFLAGS) -I $(PVM_ROOT)/include $(TOOLFLAGS) \
		 $< $(PVM_OBJECTS) -o $(BIN)/$@ $(LIBFLAGS) \
		 $(STARLAB_PVM_LIBS) -lm ; \
	    $(STARLAB_CCC) $(CFLAGS) -I $(PVM_ROOT)/include $(TOOLFLAGS) \
		 $< $(PVM_OBJECTS) -o $(BIN)/$@ $(LIBFLAGS) \
		 -D_SRC_=\"$(PWD)/$<\" \
		 $(STARLAB_PVM_LIBS) -lm ; \
	    touch $@ ; \
	fi

.C.pvmo:
	$(STARLAB_CCC) $(CFLAGS) -I$(PVM_ROOT)/include -c $< -o $@

pvm_clean:
	@/bin/rm -f *.pvm *.pvmo

#..............................................................................
#..............................................................................

##===============================================================//
##  +---------------+    _\|/_    +------------------------------\\
##  |  the end of:  |     /|\     |  src/node/dyn/hdyn/sdyn3/scatter3/Makefile
##  +---------------+             +------------------------------//
##===================== STARLAB =================================\\

