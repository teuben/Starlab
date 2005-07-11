LDFLAGS = -lm
EXECUTABLES =	check_dep earliest_time errecho getbytes ints \
		kiraformat kiraindex kiraindex2 \
		largest_int latest_time mytime readlink sqrt strip_starlab_path

# strings trn

all:	$(EXECUTABLES)
sure:	clean all

mytime:	routines.c mytime.c
	$(CC) -c -O routines.c -g3
	$(CC) -O0 mytime.c routines.o -lm -o mytime -g

exelist:
	@echo $(EXECUTABLES)

# Print a list of binary files explicitly not to be linked by mkshadow,
# in a format suitable for use by egrep.

skip:
	@for file in $(EXECUTABLES) routines.o ; do \
	    echo -n $$file'$$' ; \
	    if [ "$$file" != "routines.o" ]; then \
	        echo -n '|' ; else echo "" ; fi \
	done

maintainer-clean:
distclean:
clean:
	@/bin/rm -f $(EXECUTABLES)
	@/bin/rm -f *.o

libs:
	@echo No libraries in sbin

install:
	@echo Nothing to install in sbin

dist:
distdir:
	@echo Nothing extra to distribute in sbin

