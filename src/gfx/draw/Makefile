
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
##  src/gfx/mcdraw/draw/Makefile
##.............................................................................
##    version 1.2:	May 1994	Steve McMillan
##.............................................................................
##

#------------------------------------------------------------------------

MCDRAW_SOURCE = twodplot.f threedplot.f autocorrelation.f box.f \
		contexts.f decode.f errors.f expand.f \
		fitxy.f genfit.f \
		get_input.f get_token.f getline.f \
		getnext.f help.f histogram.f \
		integrate.f lsfit.f \
		nextcmd.f powerspectrum.f \
		readall.f readcols.f \
		readin.f rebin.f sarith.f \
		save_context.f smooth.f sort.f \
		strings.f svdfit.f swap.f \
		updhist.f varith.f wrout.f \
		xyplot.f xyzplot.f

MCDRAW_OBJECTS = $(MCDRAW_SOURCE:.f=.o)
MCDRAW_LIBRARIES   = -L$(LIB) -lmcpak -llux $(XLIB) $(SUNCORELIB) -lm

OBJECTS = mcdraw.o mcdraw_c.o ctox.o cfit.o cinter.o unix.o $(MCDRAW_OBJECTS)

#------------------------------------------------------------------------

LIBRARY =
EXECUTABLES = mcdraw

include ../Local_setup.inc

#------------------------------------------------------------------------

mcdraw:	$(OBJECTS) ../lux/.liblux.a ../libs/.libmcpak.a
	@echo Building mcdraw...
	@$(FCOMP) -o $(BIN)/$@$(EXT) $(OBJECTS) $(MCDRAW_LIBRARIES)
	@touch $@

.f:	../lux/.liblux.a ../libs/.libmcpak.a
	$(FCOMP) $(FFLAGS) $< -o $@ $(MCDRAW_LIBRARIES)

.c:	../lux/.liblux.a ../libs/.libmcpak.a
	$(ANSICC) $(CFLAGS) $< -o $@ $(MCDRAW_LIBRARIES)

# Don't optimize the driving routine (causes complaints on Sun,
# and errors on HP):

mcdraw.o:	mcdraw.f
	$(FCOMP) -c mcdraw.f

unix.o:	unix.F
	$(CPP) -P unix.F > unix.f
	$(FCOMP) $(FFLAGS) -c unix.f
	@ /bin/rm unix.f

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/gfx/mcdraw/libs/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

