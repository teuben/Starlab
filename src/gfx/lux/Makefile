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
##  src/gfx/lux/Makefile
##.............................................................................
##    version 1.1:	Feb 1994	Biao Lu (Drexel U.)
##			Ported to Starlab by Steve McMillan
##.............................................................................
##

LIBRARY = liblux.a
LIB_DEP = .$(LIBRARY)

OBJECT_FILES = win.o draw.o draw1.o color.o dialog.o mcd.o interface.o termio.o utility.o simple.o
EXECUTABLES = plot_data ximage simple

include ../Local_setup.inc

win.o:	win.c
	$(ANSICC) $(CFLAGS) -c win.c

draw.o:	draw.c
	$(ANSICC) $(CFLAGS) -O0 -g -c draw.c # Optimization ==> Linux errors...

interface.o:	interface.c
	$(ANSICC) $(CFLAGS) -c interface.c

plot_data:       plot_data.c
	$(ANSICC) $(CFLAGS) plot_data.c -L$(LIB) -llux $(XLIB) -lm -o $(BIN)/$@
	@touch $@

simple:       simple.c
	$(ANSICC) $(CFLAGS) -DTEST simple.c -L$(LIB) -llux $(XLIB) -lm -o $(BIN)/$@
	@touch $@

ximage:       ximage.c
	$(ANSICC) $(CFLAGS) ximage.c -L$(LIB) -llux $(XLIB) -lm -o $(BIN)/$@
	@touch $@

ximage_compiled.o:	ximage_compiled.c
	$(ANSICC) $(CFLAGS) -c ximage_compiled.c

hydro2:       hydro2.c ximage_compiled.o
	$(ANSICC) $(CFLAGS) hydro2.c ximage_compiled.o -L$(LIB) -llux $(XLIB) -lm -o $(BIN)/$@
	@touch $@

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/gfx/lux/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

