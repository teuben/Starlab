
Switching between force-evaluation schemes in kira
--------------------------------------------------

--> New structure, as of 6/04 (Steve).

All available kira code is now included in the compilation.  However,
some functionality may not actually exist on any given system.
Include GRAPE-4, GRAPE-6, treecode, etc. libraries by defining the
appropriate macros in Makefile.  (May be configurable?)  Use the
"stubs" library to resolve missing references.  So now only need a
single Makefile and compilation should be site-independent...

Switches between various options are included in the following
functions:

    int kira_calculate_top_level_acc_and_jerk()	  in kira_ev.C
    void kira_calculate_energies()		  in kira_energy.C
    void kira_calculate_densities()		  in kira_density.C

(We have also experimented with replacing the first of these by the
 static function pointer

    hdyn::kira_calculate_top_level_acc_and_jerk

 set up via kira_config.C, and parts of this code remain, commented
 out for now.  However, a uniform approach seems more useful.)

These functions simply make a run-time choice of which
algorithm/hardware to use, according to the static datum hdyn::config,
set by kira_config().  The configuration is set by default at startup
based on system configuration, but overridable from the command line.
We will distinguish between having a package (which can be determined
automatically) and using it (which can have a reasonable default, but
may be altered at startup -- or even during the run).

To add a new package:

* write functions conforming to the kira specs to compute top-level
  accs and jerks, energies, and densities

* add the appropriate new run-time options into the above files (may
  need additional flags, carried as static hdyn data)

* modify kira_config.C to reflect the new functionality

* provide the libraries and add them to the Makefiles in hdyn/evolve
  and hdyn/util

* provide stubs if the package is not distributed with Starlab
