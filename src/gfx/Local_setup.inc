
# Setup for all makefiles used in mcdraw construction.  Use this file
# to set up Starlab variables for use in non-Starlab environments.  The
# makefiles at lower levels should *not* make reference to any of the
# standard Starlab include files.

#========================================================================

# The following STARLAB_xxx variables must be set HERE if we are not
# operating within the Starlab environment.

# Note: This is an include file to be used by make.

#	STARLAB_PATH			root for libraries and binaries

#	STARLAB_CPP			C preprocessor
#	STARLAB_C_COMPILE		C compiler
#	STARLAB_SYSTEM_C_COMPILE	"System" C compiler
#	STARLAB_GFX_CFLAGS		Graphics-specific C-compiler flags

#	STARLAB_FORTRAN_COMPILE		FORTRAN compiler
#	STARLAB_FORTRAN_FLAGS		FORTRAN compiler flags
#	STARLAB_FORTRAN_CFLAGS		Flags for C compiler to use F77 objects

#	STARLAB_X11INC			Location of X header files
#	STARLAB_X11LIB			Location of X libraries
#	STARLAB_EXTRAS			Define any "extra" terms or options
#	STARLAB_POSIX_LIB		Link with POSIX library, if any

#	STARLAB_RANLIB			Define "ranlib" for this machine

# e.g. for Redhat 4.2 Linux, uncomment these lines:

#STARLAB_PATH = $(HOME)

#STARLAB_CPP = /usr/bin/cpp -C -Dlinux
#STARLAB_C_COMPILE = gcc
#STARLAB_SYSTEM_C_COMPILE = gcc
#STARLAB_GFX_CFLAGS =  

#STARLAB_FORTRAN_COMPILE = f77
#STARLAB_FORTRAN_FLAGS =  
#STARLAB_FORTRAN_CFLAGS = -DFORTRAN_TRAILING_UNDERSCORE

#STARLAB_X11INC = /usr/X11/include/X11
#STARLAB_X11LIB = /usr/X11/lib
#STARLAB_EXTRAS = -DHAS_USLEEP		# (also, -DHAS_NANOSLEEP)
#STARLAB_POSIX_LIB =  

#STARLAB_RANLIB = ranlib
##STARLAB_RANLIB = "echo > /dev/null"	# (use if ranlib is automatic in ar)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Location of binaries and libraries created by make:

BIN = $(STARLAB_PATH)/bin
#BIN = ..

LIB = $(STARLAB_PATH)/lib
#LIB = ..

# NOTE: The path is relative to libs and draw, so ".." refers to the
#	directory containing this setup file.

#========================================================================
#
# Probably not much reason to modify anything below this line...
#
#========================================================================

SHELL = /bin/sh

RM = /bin/rm
LS = /bin/ls
CP = /bin/cp

# Extension for new mcdraw:

EXT = 

# Compilers and flags:

FFLAGS   = -O $(STARLAB_FORTRAN_FLAGS)
CFLAGS   = -O0 -I$(STARLAB_X11INC) $(STARLAB_GFX_CFLAGS) $(STARLAB_FORTRAN_CFLAGS) $(STARLAB_EXTRAS) -O0

FCOMP = $(STARLAB_FORTRAN_COMPILE)

# Note:  gcc complains about the old mcdxsubs (not used now), and lux
# fails when compiled with gcc on hpux.

ANSICC = $(STARLAB_C_COMPILE)
CC = $(STARLAB_SYSTEM_C_COMPILE)

CPP = $(STARLAB_CPP)

#------------------------------------------------------------------------

# Set XLIB if we have X-windows:

XLIB  	= -L$(STARLAB_X11LIB) -lX11 $(STARLAB_POSIX_LIB) $(STARLAB_XLIB)
#XLIB	= 

#------------------------------------------------------------------------

# Set SUNCORELIB if we have SunCore (note that we must use static linking of
# the SunCore libraries for export, as run-time support no longer exists):

#SUNCORELIB = -Bstatic -L/usr/old -lcore77 -lcore -Bdynamic \
#	     -lsunwindow -lpixrect -lc -lm
SUNCORELIB = 

#------------------------------------------------------------------------

# Standard targets:

all:	lib bin
sure:	clean all

lib:		$(LIB_DEP)	
.$(LIBRARY):	$(OBJECT_FILES)
		@echo Building $(LIBRARY)
		ar ruv $(LIB)/$(LIBRARY) $(OBJECT_FILES)
		$(STARLAB_RANLIB) $(LIB)/$(LIBRARY)
		@touch .$(LIBRARY)

bin: $(EXECUTABLES)

tarlist:
	@make_list_files `$(LS) [A-Z]* *.[CcFfh] | grep -v "~"`

clean:
	@-$(RM) -f core $(EXECUTABLES) *.o .lib*.a *~ ...*.C .make.state

cleanbin:
	@-$(RM) -f $(EXECUTABLES)

cleanalt:
