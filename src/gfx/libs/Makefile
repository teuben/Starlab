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
##  src/gfx/mcdraw/libs/Makefile
##.............................................................................
##    version 1.2:	May 1994	Steve McMillan
##.............................................................................
##
#------------------------------------------------------------------------

# If we want to use the old X-windows interface, we must link with this:

#MCDX.o	= mcdxsubs.o
MCDX.o	= 

#------------------------------------------------------------------------

# If we have SunCore (note that we must use static linking of the
# SunCore libraries for export, as run-time support no longer exists):

#SUNCORE.o 	= mcyapp.o get_vs.o maps.o
SUNCORE.o	= 

#------------------------------------------------------------------------

# To satisfy X/SunCore/other references if the libraries don't exist:

#DUMMIES   	= ncardummies.o
DUMMIES		= sundummies.o ncardummies.o
#DUMMIES	= sundummies.o xdummies.o ncardummies.o

#------------------------------------------------------------------------

# FORTRAN routines:

MCPAK_F_OBJECTS = \
	aerrbars.o	alin.o		arrow.o		aspectratio.o \
	axes.o 		box.o		boxf.o \
	bounds.o	clrstr.o \
	cnomber.o	cnombr.o	compoz.o	contor.o \
	contor1.o	contor2.o	convec.o	csimbol.o \
	devicename.o	digital.o	displaytext.o	dplot.o \
	drline.o	eframe.o	enumbr.o	errbars.o \
	frinches.o	frlabfrm.o	frlgset.o	frlgxdr.o \
	frlgydr.o	frlndiv.o	frlnfnc.o	frlnset.o \
	frlnxdr.o	frlnydr.o	frnumbr.o	frplt.o \
	frspaces.o	frsymbl.o	frusers.o	getdevice.o \
	getp.o		getscaling.o	heart.o		hpsymout.o \
	initchars.o	jumpto.o	labels.o	lastp.o \
	lognum.o	lscale.o	mcinit.o	mcsymlims.o \
	minmax.o	mline.o		mlinec.o	modulo.o \
	ncolors.o	ngon.o		noclear.o \
	nomber.o	numbr.o		numsym.o	nxtvu.o \
	on.o		options.o	outp.o		pensize.o \
	pixel.o		plot.o		plotdetails.o	plotheight.o \
	plotin.o	plotsetup.o	plotwidth.o	plt3d.o \
	printid.o	printpage.o	psinit.o \
	ps.o		pushstr.o	reduce.o \
	save_context.o \
	setaxes.o	setfill.o	sethts.o \
	setmod.o	setorigin.o	setpat.o	setpens.o \
	setplf.o	setstars.o	settiks.o	setwts.o \
	shade.o		simbnds.o	simbol.o	simdraw.o \
	simsize.o	simwhe.o	spiral.o	square.o \
	strpos.o	symbl.o		typestring.o	uplot.o \
	unix.o		upoint.o	xyencode.o \
	$(DUMMIES)

#-----------------------------------------------------------------------

# C Routines:

MCPAK_C_OBJECTS = cutils.o

#-----------------------------------------------------------------------

MCPAK_OBJECTS = $(MCPAK_F_OBJECTS) $(MCPAK_C_OBJECTS) $(SUNCORE.o)
MCPAK_LIBRARIES   = -L$(LIB) -lmcpak -llux $(XLIB) $(SUNCORELIB) -lm

# Note: The routines in unix.F and cutils.c are used as interfaces
#	between FORTRAN and various non-standard UNIX utilities.

#-----------------------------------------------------------------------

LIBRARY = libmcpak.a
LIB_DEP = .$(LIBRARY)

OBJECT_FILES = $(MCPAK_OBJECTS)
EXECUTABLES = show_fonts testx

include ../Local_setup.inc

#-----------------------------------------------------------------------

show_fonts:	../lux/.liblux.a .libmcpak.a
	@echo Building show_fonts...
	$(FCOMP) $(FFLAGS) show_fonts.f -o $(BIN)/$@ $(MCPAK_LIBRARIES)
	@touch $@

.f:	../lux/.liblux.a .libmcpak.a
	$(FCOMP) $(FFLAGS) $< -o $(BIN)/$@ $(MCPAK_LIBRARIES)
	@touch $@

.c:	../lux/.liblux.a .libmcpak.a
	$(ANSICC) $(CFLAGS) $< -o $(BIN)/$@ $(MCPAK_LIBRARIES)
	@touch $@

#------------------------------------------------------------------------

# Special treatment for some FORTRAN/cpp files:

plotdetails.o:	plotdetails.F
	$(STARLAB_CPP) -P -UHAS_SUNCORE plotdetails.F plotdetails.f
	$(FCOMP) $(FFLAGS) -c plotdetails.f
	@/bin/rm plotdetails.f

psinit.o:	psinit.F
	$(STARLAB_CPP) -P psinit.F psinit.f
	$(FCOMP) $(FFLAGS) -c psinit.f
	@/bin/rm psinit.f

setfill.o:	setfill.F
	$(STARLAB_CPP) -P -UHAS_SUNCORE setfill.F setfill.f
	$(FCOMP) $(FFLAGS) -c setfill.f
	@/bin/rm setfill.f

unix.o:	unix.F
	$(STARLAB_CPP) -P unix.F unix.f
	$(FCOMP) $(FFLAGS) -c unix.f
	@/bin/rm unix.f

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/mcdraw/libs/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
