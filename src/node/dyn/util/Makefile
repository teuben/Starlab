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

LEX		= flex
LFLAGS		= -f

YACC		= bison
YFLAGS		= 

EXTRA_CLEAN	= dyn_lexer.C dyn_parser.C dyn_parser.H

CFLAGS		= $(STARLAB_CCC_FLAGS) -I$(INC) -DSTARLAB_USE_FDYN

#-----------------------------------------------------------------------

# To omit the fget_dyn code (e.g. lex/bison/make not working properly),
# remove -DSTARLAB_USE_FDYN from CFLAGS and comment out the following line.

FDYN_OBJECTS = dyn_lexer.o dyn_parser.o

#-----------------------------------------------------------------------

#LIBFLAGS	= -L$(LIB) -ldyn -lnode -lstd -lm

LIBFLAGS	= -L$(LIB) -ldyn -lsstar -lnode -lstd -lm

# (The -lsstar here is needed to compile sys_stats and a few other functions.
#  The star library is not listed in DEPEND_LIBS.)

DEPEND_LIBS	= $(LIB)/libdyn.a $(LIB)/libnode.a $(LIB)/libstd.a 

EXECUTABLES =	compute_com compute_mcom compute_density compute_max_cod \
		compute_mean_cod dumbp dyn_io dyndiff dyndist energy \
		extract_snap flatten freeze freezebin lagrad lagradplot \
		list_snap merge_low molecules reflect_velocities scale \
		plot_stars starplot snap2speck \
		sys_stats set_com to_com jiggle \
                lagradplotsmooth get_std_center \
		add_plummer add_power_law add_tidal \
		set_ignore_internal add_daughter_node getsnap \
		print_time radial_density radial_vdisp radial_profile \
		makemass makesecondary

OBJECT_FILES =	compute_com.o compute_density.o compute_max_cod.o \
		compute_mcom.o compute_mean_cod.o \
		dyn_di.o dyn_external.o dyn_io.o dyn_stats.o \
		dyn_story.o dyn_tt.o \
		flatten.o lagrad.o merge_low.o scale.o sys_stats.o to_com.o \
		set_com.o plot_stars.o get_std_center.o \
		refine_mass.o refine_mass2.o \
		add_plummer.o add_power_law.o \
		radial_density.o radial_vdisp.o radial_profile.o \
		wrappers.o $(FDYN_OBJECTS)

LIBRARY = libdyn.a

#..............................................................................

include $(STARLAB_PATH)/inc/Makefile.worker.inc

#..............................................................................

dyn_lexer.C:	dyn_lexer.l dyn_parser.C
	$(LEX) $(LFLAGS) -o$@ dyn_lexer.l
#	$(LEX) $(LFLAGS) -o$@ $<		# better (needs gmake)

dyn_parser.C:	dyn_parser.yc
	$(YACC) $(YFLAGS) -do $@ dyn_parser.yc
#	$(YACC) $(YFLAGS) -do $@ $<		# better (needs gmake)

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  src/node/dyn/util/Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\

