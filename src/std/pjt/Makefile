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
##  src/std/pjt/Makefile
##.............................................................................
##    Aug 1996					Peter Teuben
##.............................................................................
##

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC)
LIBFLAGS	= -L$(LIB) -lstd -lm
DEPEND_LIBS	= $(LIB)/libstd.a

EXECUTABLES = 
OBJECT_FILES = program.o getparam.o main.o \
    p_description.o p_keywords.o p_version.o p_examples.o p_usage.o

LIBRARY = libstd.a

EXTRA_CLEAN = $(LIB)/$(LIBRARY)

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................
#..............................................................................

# Local targets:

TAR = Makefile program.h program.C getparam.C main.C \
    p_description.C p_keywords.C p_version.C p_examples.C p_usage.C \
    test1.C nemo.sty latest.tex ui.inc README.pjt

tar:
	cp program.h ../../../inc/
	tar cvf pjt.tar $(TAR)

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/std/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

