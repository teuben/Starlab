
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// set_com:  Modify positions and velocities to set the center-of-mass
////           position and velocity.  Write com_pos and com_vel to the root
////           dyn story.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -n    force interpretation of r and v in N-body units [no]
////              -r    specify center of mass position [(0,0,0)]
////              -v    specify center of mass velocity [(0,0,0)]
////
//// If an external field has been specified, the velocity is taken to be in
//// units of the circular orbit speed at the specified location.  Otherwise,
//// the velocity is taken as is.  Positions and velocities are in physical
//// units (parsecs and km/s, if relevant) if physical stellar parameters
//// are known, and N-body units otherwise (or if -n is specified).

//   version 1:  Aug/Sep 2001   Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void dyn::set_com(vector pos, vector vel)	// defaults = 0
{
    vector com_pos;
    vector com_vel;

    compute_com(this, com_pos, com_vel); 

    vector dpos = pos - com_pos;
    vector dvel = vel - com_vel;

    // Adjust only top-level nodes; don't touch the root node.

    for_all_daughters(dyn, this, bb) {
	bb->inc_pos(dpos);
	bb->inc_vel(dvel);
    }

    // Correct entries in the dyn story.

    putvq(get_dyn_story(), "com_pos", pos);
    putvq(get_dyn_story(), "com_vel", vel);

    // Note: Do not modify the "center" of any external field.
    // This function adjusts only the N-body center of mass.
}

#else

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    vector r = 0;
    vector v = 0;

    bool n_flag = false;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    char* param_string = "c:nr:::v:::";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
	    case 'n':	n_flag = true;
			break;

	    case 'r':	r = vector(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;
	    case 'v':	v = vector(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;

            case '?':	params_to_usage(cerr, argv[0], param_string);
	    		get_help();
	    		exit(1);
        }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	// Check if we have to reinterpret r and v in light of possible
	// external fields and physical parameters.

	real mass, length, time;
	bool phys = get_physical_scales(b, mass, length, time);

	// If phys, we are using physical units.  Mass, length, and time are
	// the physical equivalents of 1 N-body unit, in Msun, pc, and Myr.

	check_set_external(b, false);
	bool ext = (b->get_external_field() != 0);

	// If ext, we have an external field, assumed attractive.  Note
	// that this function doesn't make much sense if ext = false...
	// If ext, the external field has already been converted to
	// N-body units, regardless of phys.

	// Possibilities:
	//
	// no ext field, no physical units:	set r, v as is in N-body units
	// no ext field, physical units:	convert r, v to N-body and set
	// ext field, no physical units:	use N-body units, but use v/vc
	// ext field, physical units:		convert to N-body and use v/vc

	if (phys && !n_flag) {

	    cerr << "set_com:  converting input physical parameters"
		 << " to N-body units" << endl;
	    cerr << "          r = (" << r << ") pc,  v = ("
		 << v << ") v_circ" << endl;
	    cerr << "          N-body length scale = " << length << " pc"
		 << endl;

	    // Convert r and v to N-body units.

	    r /= length;
	    if (!ext) {

		// Input v is in km/s; length/time is the velocity unit
		// in pc/Myr = 3.086e13/3.156e13 km/s = 0.978 km/s.

		real vunit = 0.978*length/time;
		v /= vunit;
	    }

	} else

	    cerr << "set_com:  interpreting input parameters"
		 << " as N-body units" << endl;

	if (ext) {

	    // Convert v from units of vcirc into N-body units.

	    real vc = vcirc(b, r);

	    if (vc > 0) {
		v *= vc;
		cerr << "          circular orbit speed = " << vc;
		if (phys && !n_flag) cerr << " (N-body)";
		cerr << ",  period = " << 2*M_PI*abs(r)/vc
		     << endl;
	    }
	}

	// Now r and v are in N-body units, appropriately scaled.

        b->set_com(r, v);
	cerr << "          r = (" << r << "),  v = (" << v << ")"
	     << endl;

	put_dyn(b);
	rmtree(b);
    }
}

#endif
