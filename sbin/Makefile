LDFLAGS = -lm
EXECUTABLES =	check_dep earliest_time kiraformat kiraindex largest_int \
		latest_time mytime sqrt strip_starlab_path

# strings trn

all:	$(EXECUTABLES)
sure:	clean all

.c:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.o:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.c.o:
	$(STARLAB_C_COMPILE) -c $<

mytime:	mytime.c
	$(STARLAB_C_COMPILE) -O0 mytime.c -lm -o mytime

hide_exe:
	@./hide $(EXECUTABLES)

show_exe:
	@./show $(EXECUTABLES)

clean:
	@/bin/rm -f $(EXECUTABLES)
