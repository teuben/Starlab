LDFLAGS = -lm
EXECUTABLES = check_dep earliest_time latest_time sqrt kiraformat # strings trn

.c:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.o:
	$(STARLAB_C_COMPILE) $< -o $* -lm

.c.o:
	$(STARLAB_C_COMPILE) -c $<


all:	$(EXECUTABLES)

hide_exe:
	@./hide $(EXECUTABLES)

show_exe:
	@./show $(EXECUTABLES)

clean:
	@/bin/rm -f $(EXECUTABLES)
