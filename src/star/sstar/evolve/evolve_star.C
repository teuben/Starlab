//// evolve_star: evolve a cluster of single stars (see nstarev).
////              the single stars should be provided in the input stream.
////
//// Options:    -c    create a separate instear of taking it from input.
////             -d    time interval for output.
////             -M    Total mass of the stars.
////             -m    See -M option.
////             -n    number of output timesteps (timesteps are taken
////                   with constant time intervals) 
//++ Notes:
//++  Assymetry in supernova is taken care of by integrator (see starev).
//++
//++ Example of usage:      
//++  evolve_star -c -d 2 -m 10 -t 30
//++
//++ See also: addstar
//++           mkmass
//++           mknode
//++           nstarev
//++           starev


//// Latest version (SPZ:1.1) February 1993.

//   version 1:        1997   Steven McMillan
//   version 1.1: July 1998   Simon F. Portegies Zwart
//                            spz@grape.c.u-tokyo.ac.jp

// evolve_star:	Handy little program using Simon's stellar evolution
//		package to track the evolution of one or more stars.
//
//		NOTE: There is no particular reason to keep this
//		      program in the hdyn section of starlab...
//
//		      Dyn is necessary here because of flatten_node and
//		      system_time, and also because stellar_evolution()
//		      expects a dyn pointer.

#include "dyn.h"
#include <star/sstar_to_dyn.h>
#include <star/dstar_to_dyn.h>

main(int argc, char **argv)
{
    // Defaults:

    int n = 1;					// 1 star
    real m = 1.0;				// 1 solar mass
    real t = 0.0, dt = 10.0, t_end = 1000.0;	// 1 Gyr, 10 Myr steps

    bool create_system = false;			// expect system on stdin

    extern char *poptarg;
    int c;
    char* param_string = "cd:m:M:n:ct:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case 'c':	create_system = true;
			break;
	    case 'd':	dt = atof(poptarg);
			break;

	    case 'm':
	    case 'M':	m = atof(poptarg);	// Note: M, m, and n are used
			break;			// only if we are creating
	    case 'n':	n = atoi(poptarg);	// our own system ("-c").
			break;

	    case 't':	t_end = atof(poptarg);
			break;
	    default:
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}
    }

    dyn *b;					// root node

    if (create_system) {

	if (n < 1 || m <= 0.0) err_exit("Need N and M > 0");
	b = (dyn*) mknode_mass(n, m);

    } else {

	b = get_dyn();
	if (b == NULL) err_exit("Can't read input snapshot");

	n = 0;
	m = 0;
	b->flatten_node();
	for_all_daughters(dyn, b, bb) {
	    n++;
	    m += bb->get_mass();
	}

    }

    // Now total number of star nodes is n, total mass is m (solar),
    // in *all* cases.

    // Set stellar evolution mass unit to 1 solar mass, length unit
    // to 1 pc, time unit to 1 Myr:

    b->get_starbase()->set_stellar_evolution_scaling(1.0, 1.0, 1.0);

    // Start star(s) on the Main Sequence:

    addstar(b, 0.0, Main_Sequence);

    // End of initialization.  Now loop over evolution.

    while (t < t_end - 0.0001*dt) {

	int old_id_sum = 0;
	{for_all_daughters(dyn, b, bb) {
	    star* s = (star*) bb->get_starbase();
	    old_id_sum += s->get_element_type();
	}}

	// Evolve the system through one timestep.

	// NOTE that the stellar evolution routines expect
	// system_time and node->mass to be properly set.

	b->set_system_time(t += dt);
	stellar_evolution(b);

	int new_id_sum = 0;
	{for_all_daughters(dyn, b, bb) {

	    star* s = (star*) bb->get_starbase();
	    new_id_sum += s->get_element_type();

	    // Routinely update the mass.

	    bb->set_mass(s->get_total_mass());	// (Note: get_total_mass()
						//	  is in dyn units)
	}}

	// Change in identity of any star signals output from single_star_io.
	// Extra line feed here just cleans up output...

	if (new_id_sum != old_id_sum) cerr << endl;

	// Current verbose output on every star at every step:

	int count = 0;
	{for_all_daughters(dyn, b, bb) {

	    if (count++ == 0)
		fprintf(stderr, "time = %9.2f  ", t);
	    else
		fprintf(stderr, "                  ");

	    star* s = (star*) bb->get_starbase();
	    single_star* ss = (single_star*) s;

	    ss->dump(cerr, false);
	    
	    cerr << bb->format_label() << "  "
		 << type_string(s->get_element_type())
		 << "  M = " << s->get_total_mass()
		 << "  Mcore = " << s->get_core_mass()
		 << "  L = " << ss->get_luminosity()
		 << "  V = " << ss->magnitude()
		 << "  T = " << ss->temperature()
		 << "  R = " << ss->get_radius()
		 << endl;
	}}

    }
}
