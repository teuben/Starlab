
//// compute_mean_cod:  Determine the density center position and velocity
////                    for the input N-body system.  The density center is
////                    defined analogously to the center of mass, but instead
////                    of a mass weighting, the weighting used is proportional
////                    to the square of the density (the original suggestion
////                    in the paper below(*), to use a weighting factor
////                    linear in the density, may not converge well).
////
////                    Note: The computed density center is defined relative in
////                    absolute terms, and so includes the pos and vel of the
////                    parent node.
////
////                    Densities are not computed here -- run compute_density
////                    before invoking compute_mean_cod.
////
////                    Density center position and velocity are written to
////                    the dyn story of the top-level node; they are also
////                    optionally returned as function arguments in the
////                    library version.
////
//// Options:     -c    add a comment to the output snapshot [false]
////
//-----------------------------------------------------------------------------
//   version 1:  Nov 1994   Piet Hut
//   version 2:  Jul 1996   Steve McMillan & Jun Makino
//.............................................................................
//// (*) Stefano Casertano and Piet Hut: Astroph.J. 298,80 (1985).
////     To use their recipe, k >= 2 is required.
//.............................................................................
//   non-local function: 
//      compute_mean_cod
//.............................................................................
//   see also: density.C
//-----------------------------------------------------------------------------

#include "dyn.h"

#ifndef TOOLBOX

//-----------------------------------------------------------------------------
//  compute_mean_cod -- Returns the position and velocity of the density
//		        center of an N-body system.
//-----------------------------------------------------------------------------

#define MAX_COUNT 5

void compute_mean_cod(dyn *b, vec& pos, vec& vel)
{
    real total_weight = 0;
    bool print_message = true;
    int count = 0;

    pos = 0;
    vel = 0;
    
//  for_all_leaves(dyn, b, d)
    for_all_daughters(dyn, b, d) {

	real dens_time = getrq(d->get_dyn_story(), "density_time");

	if (print_message
	    && !twiddles(dens_time, (real) b->get_system_time(), 1.e-9)) {
	    warning("compute_mean_cod: using out-of-date densities.");
	    PRL(d->format_label());
	    int p = cerr.precision(HIGH_PRECISION);
	    PRL(b->get_system_time());
	    PRL(dens_time);
	    cerr.precision(p);
	    if (++count > MAX_COUNT) print_message = false;
	}

	real this_density = getrq(d->get_dyn_story(), "density");

	if (this_density > 0) {
	    real dens2 = this_density * this_density;	    // weight factor
	    total_weight += dens2;
	    pos += dens2 * d->get_pos();
	    vel += dens2 * d->get_vel();
	} else if (this_density <= -VERY_LARGE_NUMBER) {
	    warning("compute_mean_cod: density not set.");
	    PRL(d->format_label());
	    if (++count > MAX_COUNT) print_message = false;
	}
    }	

    if (total_weight > 0) {
	pos /= total_weight;
	vel /= total_weight;
    }

    // Include the parent quantities.

    pos += b->get_pos();
    vel += b->get_vel();

    putsq(b->get_dyn_story(), "density_center_type", "mean");
    putrq(b->get_dyn_story(), "density_center_time", b->get_system_time());
    putvq(b->get_dyn_story(), "density_center_pos", pos);
    putvq(b->get_dyn_story(), "density_center_vel", vel);
}

void compute_mean_cod(dyn *b)
{
    vec pos, vel;
    compute_mean_cod(b, pos, vel);
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use compute_mean_cod() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    dyn * b;
    bool  c_flag = FALSE;       // if TRUE: a comment given on command line

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
        }            

    if ((b = get_dyn()) == NULL)
       err_exit("compute_mean_cod: No N-body system on standard input");

    while (b) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        compute_mean_cod(b);

        put_dyn(b);
	rmtree(b);
	b = get_dyn();
    }
}

#endif


// endof: compute_mean_cod.C
