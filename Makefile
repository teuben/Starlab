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
##  Makefile
##.............................................................................
##    Dec 1994					Piet Hut & Steve McMillan
##.............................................................................
##

SHELL = /bin/sh

RM = /bin/rm
MV = /bin/mv
CP = /bin/cp
DIFF = diff

MAKEWALK = makewalk

DIRS = etc src

#..............................................................................

# if make is invoked without argument, the following help message is displayed:

help:
	@echo This Makefile offers the following make options:
	@echo '  all:  executes    make libs  and  make bins'
	@echo '  libs: executes a  make lib  everywhere in src'
	@echo '  bins: executes a  make bin  everywhere in src'
	@echo '  alt:  executes a  make all  in alt'
	@echo '  usr: executes a  make all  in usr'
	@echo '  tar: makes tarfiles at all levels; the files will be named'
	@echo '       yymmdd.s_tar_lab_x with yymmdd being the current year,'
	@echo '       month, day, and x being the tar level'
	@echo '  tar0: tar level 0: makes a tarfile for the core of starlab'
	@echo '  tar1: tar level 1: same, plus the NBODY fortran files'
	@echo '  tar2: tar level 2: same, plus the /usr contributions'
	@echo '  clean: removes the tarfile tar.starlab, removes core file'
	@echo '         and emacs backup versions, here, in inc, demo,'
	@echo '         deletes all files in lib and bin (which means that all'
	@echo '         FOREIGN FILES, i.e. non-starlab files WILL BE LOST),'
	@echo '         and executes a  make clean  everywhere in src'
	@echo '  cleanbin: removes all files in bin and executes a'
	@echo '         make cleanbin  everywhere in src'

#..............................................................................

sure:	clean all

#..............................................................................

all:	sbins .libs bins
	@echo all done

#..............................................................................

sbins:
	@echo '==>' making sbin
	@cd sbin ; ${MAKE} $(STARLAB_MAKEFLAGS) all ; cd ..
	@if [ -d etc ] ; then true ; else mkdir etc ; fi
	@if [ -d etc/lab ] ; then true ; else mkdir etc/lab ; fi

#..............................................................................

.libs:
	@${MAKE} libs

libs:
	@if [ -d lib ] ; then true ; else mkdir lib ; fi
	@cd src ; ${MAKE} $(STARLAB_MAKEFLAGS) lib ; cd ..
	@touch .libs
	@echo libs done

lib:	libs

#..............................................................................

# The .libs dependency here should force all libs to be made in their
# entirety before any bins are built (parallel make).

bins:	.libs
	@if [ -d bin ] ; then true ; else mkdir bin ; fi
	@cd sbin ; ${MAKE} $(STARLAB_MAKEFLAGS) all ; cd ..
	@cd src ; ${MAKE} $(STARLAB_MAKEFLAGS) bin ; cd ..
	@echo bins done

bin:	bins

#..............................................................................

alt:
	@cd src/alt ; ${MAKE} $(STARLAB_MAKEFLAGS) all ; cd ../..
	@cd src/hydro/alt ; ${MAKE} $(STARLAB_MAKEFLAGS) all ; cd ../../..

#..............................................................................

usr:
	@cd usr ; ${MAKE} $(STARLAB_MAKEFLAGS) all ; cd ..

#..............................................................................

# Note: rather than simply copying all files in the directories below,
#	only files of specific type are copied in order to avoid the
#	inclusion of some files such as the zero-length placeholders
#	for the executables.

#	The nbody1 and nbody5 packages are themselves saved as tarfiles,
#	and the directories are pruned.

#	Options "tar1" and "tar2" may produce very large files, so
#	they are automatically compressed, as of 4/14/93.

#	Option tar0 is now compressed also (3/12/94).

#	As of 10/4/96, the "tar0" options check to see if a new tarfile
#	is actually needed before creating one.

COMPRESS = gzip
CEXT = gz

#..............................................................................

tar:	tar0

# The "tar0" option will save all files in . and below inc, src, and etc,
# that begin with a capital letter or start with a lowercase letter or "_"
# end in ".C", ".c", ".h", ".f", or ".F".

# Note: all "alt" directories, and the top-level "fp" directory, are
#       excluded, as is etc/lab and everything below it.  This is to
#	avoid overwriting work in progress.  However, etc/demo is
#	explicitly saved.  It is assumed that no source files are
#	kept in etc, so none are saved.

# The script check_make0 ensures that we only create a new tarfile when
# the old one really is out of date.

# Use hide_exe and show_exe in sbin to avoid carrying around executables.

tar0:
	@${MAKE} list_tar0 | grep -v 'make\[.\]' \
			| grep -v cxx_repository > tar0_tmp_file
	@check_make0 tar0_tmp_file
	@rm -f tar0_tmp_file

# *** Note special treatment of hdyn/evolve and hdyn/util Makefiles. ***

do_tar0:
	@$(CP) src/node/dyn/hdyn/evolve/Makefile local/Makefile.kira
	@$(CP) local/cshrc.starlab templates/cshrc.starlab.export
	@$(MV) src/node/dyn/hdyn/evolve/Makefile \
		src/node/dyn/hdyn/evolve/Makefile.export
	@$(MV) src/node/dyn/hdyn/util/Makefile \
		src/node/dyn/hdyn/util/Makefile.export
	@(cd sbin ; ${MAKE} hide_exe > /dev/null) >/dev/null
	@tarfile=`date '+%y%m%d%H'`.s_tar_lab_$(STARLAB_VERSION)_0.$(CEXT); \
		echo Tar file is $$tarfile; \
		${MAKE} list_tar0| grep -v 'make\[.\]' > tar0_tmp_file;	\
		sbin/do_big_tar tar0_tmp_file  | $(COMPRESS) -c > $$tarfile
	@(cd sbin ; ${MAKE} show_exe > /dev/null) >/dev/null
	@$(MV) src/node/dyn/hdyn/evolve/Makefile.export \
		src/node/dyn/hdyn/evolve/Makefile
	@$(MV) src/node/dyn/hdyn/util/Makefile.export \
		src/node/dyn/hdyn/util/Makefile

# Save:	(1) [A-Z]* make_log doc templates sbin
#	(2) inc, src:	[A-Z]*			(exclude alt)
#	(3) inc, src:	[_a-z]*.[CcFfh]		(exclude alt)
#	(5) etc:	[A-Z]*			(exclude lab)

list_tar0:
	@echo [A-Z]* doc templates sbin/[a-z]* sbin/[A-Z]* \
			| tr ' ' \\12
	@find inc src -name alt -prune -o \
		-name '*.o' -o \
		\( -name '[A-Z]*' -print \) -o \
		\( -name '[_a-z]*.[CcFfh]' -print \) -o \
		\( -name '*.cc' -print \) -o \( -name '*.cpp' -print \)
	@find etc     -name lab -prune -o -name '[A-Z]*'          -print

# Make tarfiles of node or star source only:

tar_node:
	@tar cf - `echo inc/*.h; find src/node -name '[_a-z]*.[CcFfh]' -print` \
		  | $(COMPRESS) -c > tar_node.gz 2> /dev/null

tar_star:
	@tar cf - `echo inc/star; find src/star -name '[_a-z]*.[CcFfh]' -print`\
		  | $(COMPRESS) -c > tar_star.gz 2> /dev/null

#..............................................................................

autosave:	autosave0

# autosave0: same as tar0, but with a standard file name:

autosave0:
	@check_autosave0 `${MAKE} list_tar0 | grep -v 'make\[.\]'`

do_autosave0:
	@(cd sbin ; ${MAKE} hide_exe)
	@(tar cf - `${MAKE} list_tar0 | grep -v 'make\[.\]'`		\
		| $(COMPRESS) -c > autosave$(STARLAB_VERSION).tar0.$(CEXT)) \
		2> /dev/null
	@(cd sbin ; ${MAKE} show_exe)

#..............................................................................

tar1:
	@(cd sbin ; ${MAKE} hide_exe)
	@if [ -d alt ] ; then (cd src/alt; $(MAKEWALK) tar nbody1 nbody5) ; fi
	@(tar cf - `${MAKE} list_tar1 | grep -v 'make\[.\]'` 		\
		| $(COMPRESS) -c > 					\
		  `date '+%y%m%d%H'`.s_tar_lab_$(STARLAB_VERSION)_1.$(CEXT)) \
		2> /dev/null
	@if [ -d alt ] ; then (cd src/alt; $(RM) -f *.tar > /dev/null) ; fi
	@(cd sbin ; ${MAKE} show_exe)

list_tar1:
	@echo [A-Z]* ${MAKE} templates doc sbin/[a-z]* sbin/[A-Z]*
	@find inc src -name 'nbody?' -prune -o -name '[A-Z]*' -print
	@find inc src -name 'nbody?' -prune -o -name '*.[CcFfh]' -print
	@find etc -name lab -prune -o -name '[A-Z]*' -print

#..............................................................................

tar2:
	@(cd sbin ; ${MAKE} hide_exe)
	@if [ -d alt ] ; then (cd src/alt; $(MAKEWALK) tar nbody1 nbody5); fi
	@(tar cf - `${MAKE} list_tar2 | grep -v 'make\[.\]'`		\
		| $(COMPRESS) -c > 					\
		  `date '+%y%m%d%H'`.s_tar_lab_$(STARLAB_VERSION)_2.$(CEXT)) \
		2> /dev/null
	@if [ -d alt ] ; then (cd src/alt; $(RM) -f *.tar > /dev/null) ; fi
	@(cd sbin ; ${MAKE} show_exe)

list_tar2:
	@echo [A-Z]* make_log templates doc sbin/[a-z]* sbin/[A-Z]*
	@find inc src usr -name 'nbody?' -prune -o -name '[A-Z]*' -print
	@find inc src usr -name 'nbody?' -prune -o -name '*.[CcFfh]' -print
	@find etc -name lab -prune -o -name '[A-Z]*' -print

#..............................................................................

# autosave2: same as tar2, but with a standard file name:

autosave2:
	@(cd sbin ; ${MAKE} hide_exe)
	@if [ -d alt ] ; then (cd src/alt; $(MAKEWALK) tar nbody1 nbody5); fi
	@(tar cf - `${MAKE} list_tar2 | grep -v 'make\[.\]'`		\
		| $(COMPRESS) -c > autosave$(STARLAB_VERSION).tar2.$(CEXT)) \
		2> /dev/null
	@if [ -d alt ] ; then (cd src/alt; $(RM) -f *.tar > /dev/null) ; fi
	@(cd sbin ; ${MAKE} show_exe)

#..............................................................................

# tarmin: minimal tar storage:

#	1. Basic top-level files.
#	2. Selected include files.
#	3. sbin utility files (no emacs backups)
#	4. src/std.
#	5. src/node, but not src/node/dyn and below.

tarmin:
	@(cd sbin ; ${MAKE} hide_exe)
	@(tar cf - Makefile templates README doc			    \
		`cat inc/min.inc`					    \
		`find sbin/* 	-name "*[~#]" -prune -o			    \
				\( -name "[a-z]*" -o -name "[A-Z]*" \)      \
				-print`					    \
		 src/Makefile						    \
		`find src/std	\( -name '*.[CcFfh]' -o -name '[A-Z]*' \)   \
				-print`					    \
		`find src/node 	-name dyn -prune -o			    \
				\( -name '*.[CcFfh]' -o -name '[A-Z]*' \)   \
				-print`					    \
		| $(COMPRESS) -c >					    \
		  `date '+%y%m%d%H'`.s_tar_lab_$(STARLAB_VERSION)_min.$(CEXT)) \
		2> /dev/null
	@(cd sbin ; ${MAKE} show_exe)

#..............................................................................

# tar_curr: Make a tarfile for the current "build" configuration.

tar_curr:
	@(cd sbin ; ${MAKE} hide_exe)
	@(cd templates ; $(MV) cshrc.starlab .cshrc.starlab; \
		         $(CP) ../local/cshrc.starlab .)
	@(tar cf - `${MAKE} tarlist`					 \
		| $(COMPRESS) -c > 					 \
		`date '+%y%m%d%H'`.s_tar_lab_$(STARLAB_VERSION).$(CEXT)) \
		2> /dev/null
	@(cd templates ; $(MV) .cshrc.starlab cshrc.starlab)
	@(cd sbin ; ${MAKE} show_exe)

# List files to tar, given current "build" preferences.

tarlist:
	@echo Makefile templates README doc etc/Makefile etc/demo
	@find sbin/* 	-name "*[~#]" -prune -o				\
			\( -name "[a-z]*" -o -name "[A-Z]*" \)		\
			-print
	@echo inc
	@$(MAKEWALK).quiet tarlist src

#..............................................................................
#  list all files you have newer than your CVS archived

cvsu:
	cvsu | grep ^M

#..............................................................................

clean:
	@-$(RM) -f .libs core *~ bin/* lib/* inc/*~
	@-$(RM) -f `find . -name .make.state -print`
	@$(MAKEWALK) clean sbin $(DIRS)

#..............................................................................

cleanalt:
	@echo '==>' descending to $(PWD)/src, to execute' ' "${MAKE} cleanalt"
	@cd src ; ${MAKE} cleanalt ; cd ..

#..............................................................................

cleanbin:
	-$(RM) -f bin/*
	@echo '==>' descending to $(PWD)/src, to execute' ' "${MAKE} cleanbin"
	@cd src ; ${MAKE} cleanbin ; cd ..

#..............................................................................
# configure and related stuff for new starlab4 installation

CONFIG_CLEAN = config.h config.cache config.log config.status 
CONFIG_EXTRA = cshrc.starlab config.h 

config: configure
	./configure


config_extra:	config_new

#			generic site is called 'linux', works on most machines
SITE = linux

config_new:
	-mkdir bin lib local	
	$(CP) cshrc.starlab local
	$(CP) config.h inc
	$(CP) templates/starlab_setup local
	-(cd src/node/dyn/hdyn/util; cp Makefile.$(SITE) Makefile)
	-(cd src/node/dyn/hdyn/evolve; cp Makefile.$(SITE) Makefile)

find_site:
	find . -name Makefile.$(SITE) -print

diff_new:
	-$(DIFF) cshrc.starlab local
	-$(DIFF) config.h inc
	-$(DIFF) templates/starlab_setup local

config.status: configure
	$(SHELL) ./config.status --recheck

configure: configure.in 
	autoconf

config_clean:
	rm -f $(CONFIG_CLEAN)


#		a CVS based distribution maker

STARLAB_VERSION = `cat VERSION`
DIST_DIR = starlab_$(STARLAB_VERSION)

version:	VERSION
	echo '#define STARLAB_VERSION  "'$(STARLAB_VERSION)'"' > inc/version.h
	

dist:
	rm -rf $(DIST_DIR)
	cvs -q export -D tomorrow -d $(DIST_DIR) starlab 2>&1 > /tmp/starlabdist.log
	tar -cf $(DIST_DIR).tar $(DIST_DIR)
	gzip $(DIST_DIR).tar
	rm -rf $(DIST_DIR)

#..............................................................................
#..............................................................................

##=======================================================================//
##  +---------------+        _\|/_        +------------------------------\\
##  |  the end of:  |         /|\         |  Makefile
##  +---------------+                     +------------------------------//
##========================= STARLAB =====================================\\
