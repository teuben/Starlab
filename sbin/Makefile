LDFLAGS = -lm
EXECUTABLES =	check_dep earliest_time ints kiraformat kiraindex largest_int \
		latest_time mytime readlink sqrt strip_starlab_path

# strings trn

all:	$(EXECUTABLES)
sure:	clean all

mytime:	routines.c mytime.c
	$(CC) -c -O routines.c -g3
	$(CC) -O0 mytime.c routines.o -lm -o mytime -g

distclean:
clean:
	@/bin/rm -f $(EXECUTABLES)
	@/bin/rm -f *.o

libs:
	@echo No libraries in sbin

install:
	@echo Nothing to install in sbin

