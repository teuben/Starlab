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
##  src/node/dyn/hdyn/sdyn3/scatter3/interfaces/Makefile
##.............................................................................
##    Sep 1995					Steve McMillan
##.............................................................................

CFLAGS		= $(STARLAB_CCCFLAGS) -I$(INC) -O $(STARLAB_FORTRAN_CFLAGS)
LIBFLAGS	= -L$(LIB) -lsdyn3 -lhdyn -ldyn -lnode -lstd \
	   	  $(STARLAB_CCC_LIBPATH) $(STARLAB_CCC_LIB) -lm -lg

# Don't use EXECUTABLES, as those may be assumed to be C++!

NON_CCC_EXECUTABLES =	cscatter3 cscat3_template fscat3_template

OBJECT_FILES =  
LIBRARY = 

#..............................................................................

#.c.o:
#	$(STARLAB_C_COMPILE) -I$(INC) -c $<

cscatter3:	cscatter3.o c_interface.o
	$(STARLAB_C_COMPILE) -I$(INC) -o $@ \
		cscatter3.o c_interface.o $(LIBFLAGS)

cscat3_template:	cscat3_template.o c_interface.o
	$(STARLAB_C_COMPILE) -I$(INC) -o $@ \
		cscat3_template.o c_interface.o $(LIBFLAGS)

#..............................................................................

.f.o:
	$(STARLAB_FORTRAN_COMPILE) $(STARLAB_FORTRAN_FLAGS) -c $<

fscat3_template:	fscat3_template.o f_interface.o c_interface.o
	$(STARLAB_FORTRAN_COMPILE) $(STARLAB_FORTRAN_FLAGS) -o $@ \
		fscat3_template.o f_interface.o c_interface.o $(LIBFLAGS)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

##===============================================================//
##  +---------------+    _\|/_    +------------------------------\\
##  |  the end of:  |     /|\     |  .../sdyn3/scatter3/interfaces/Makefile
##  +---------------+             +------------------------------//
##===================== STARLAB =================================\\

