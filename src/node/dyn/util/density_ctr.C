/*
 *  density_ctr.C: coordinate transformation to the density center for
 *                 a Nbody system, defined as the position of the particle
 *                 with the highest local density
 *.............................................................................
 *    version 1:  May 1989   Piet Hut
 *    version 2:  Nov 1994   Piet Hut
 *.............................................................................
 *  non-local function: 
 *    density_center
 *.............................................................................
 *     Computes the position and velocity of the density center of a Nbody
 *  system, and shifts all positions and velocities so that the root node
 *  coincides with the density center (the shifted values are absorbed in the
 *  positions and velocities of the root node).
 *.............................................................................
 *  see also: density.C
 *.............................................................................
 */

#include "dyn.h"

/*-----------------------------------------------------------------------------
 *  density_center --  Computes the total mass of an Nbody system, and the
 *                     position and velocity of its density center.
 *                     The shift in position and velocity of the density
 *                     center is substracted from each individual body, so
 *                     that for each body the absolute position (density
 *                     center position + its own position) remains 
 *                     invariant; the same holds for the absolute velocity.
 *-----------------------------------------------------------------------------
 */
void  density_center(dyn *b)
    {
    real this_density;
    real max_density = 0;
    vector pos_shift;
    vector vel_shift;

    for_all_leaves(dyn, b, d)
	{
	this_density = getrq(d->get_dyn_story(), "density");

	if (max_density < this_density)
	    {
	    max_density = this_density;
	    pos_shift = - something_relative_to_root(d, &dyn::get_pos);
	    vel_shift = - something_relative_to_root(d, &dyn::get_vel);
	    }
	}

    for_all_daughters(dyn, b, dd)
	{
	dd->inc_pos(pos_shift);
	dd->inc_vel(vel_shift);
	}

    b->inc_pos(-pos_shift);
    b->inc_vel(-vel_shift);
    }

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver to use  density_center() as a tool
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
    {
    int  c;
    char  *comment;
    extern char *poptarg;
    dyn * b;
    bool  c_flag = FALSE;       /* if TRUE: a comment given on command line  */
    int  pgetopt(int, char **, char *);

    while ((c = pgetopt(argc, argv, "c:")) != -1)
	switch(c)
	    {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': err_exit("usage: density_center [-c]");
	    }            

    if ((b = get_dyn(cin)) == NULL)
       err_exit("density_center: No N-body system provided on standard input");

    while (b)
        {
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        density_center(b);

        put_dyn(cout, *b);
	rmtree(b);
	b = get_dyn(cin);
        }
    }

#endif

/*===========================================================================*/

/* endof: density_ctr.C */
