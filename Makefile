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
##    Jul 2001					Steve McMillan & Peter Teuben
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

# If make is invoked without argument, the following help message is displayed:

help:
	@echo This Makefile offers the following make options:
	@echo '  all:  executes    make libs  and  make bins'
	@echo '  libs: executes a  make lib  everywhere in src'
	@echo '  bins: executes a  make bin  everywhere in src'
	@echo '  dist: creates a current CVS snapshot (net access needed)'
	@echo '  clean: removes the tarfile tar.starlab, removes core file'
	@echo '         and emacs backup versions, here, in inc, demo,'
	@echo '         deletes all files in lib and bin (which means that all'
	@echo '         FOREIGN FILES, i.e. non-starlab files WILL BE LOST),'
	@echo '         and executes a  make clean  everywhere in src'
	@echo '  cleanbin: removes all files in bin and executes a'
	@echo '         make cleanbin  everywhere in src'

#..............................................................................

install:
	@echo To install starlab:
	@echo '    ./configure'
	@echo '    make config_new'
	@echo '    source starlab_start   (or  . ./starlab_start.sh)'
	@echo '    make all'

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

packages:
	@cd src/packages; ./install; cd ../..

#..............................................................................

# Note from Steve:  All "tar" options have been removed as of July 2001
#                   with the long-awaited release of Starlab 4.0.

#..............................................................................

#  List all files you have newer than your CVS archived

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

CONFIG_CLEAN = config.h config.cache config.log config.status \
		cshrc.starlab starlab_start starlab_start.sh
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
	$(CP) templates/starlab_setup.sh local
	$(CP) templates/starlab.sh local
	@echo "###  Note: patching some Makefiles for SITE=$(SITE)"
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
