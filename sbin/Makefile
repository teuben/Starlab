LDFLAGS = -lm
EXECUTABLES = check_dep earliest_time kiraformat largest_int latest_time \
	      mytime sqrt

# strings trn

.c:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.o:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.c.o:
	$(STARLAB_C_COMPILE) -c $<

mytime:	mytime.c
	$(STARLAB_C_COMPILE) -O0 mytime.c -lm -o mytime

all:	$(EXECUTABLES)

hide_exe:
	@./hide $(EXECUTABLES)

show_exe:
	@./show $(EXECUTABLES)

clean:
	@/bin/rm -f $(EXECUTABLES)
