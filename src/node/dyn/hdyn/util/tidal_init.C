
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// tidal_init.C:  initialization of tidal parameters.
//
// Externally visible functions:
//
//	real get_initial_mass
//	real get_initial_virial_radius
//	real get_initial_jacobi_radius
//	void set_tidal_params

#include "hdyn.h"

real get_initial_mass(hdyn* b,
		      bool verbose)	// default = false
{
    real initial_mass = getrq(b->get_log_story(), "initial_mass");

    if (initial_mass > 0) {

	// The input snapshot has the information needed.

	if (verbose)
	    cerr << "Read initial_mass = " << initial_mass
		 << " from input snapshot" << endl;

    } else {

	// Note: presently no option to set initial_mass from the command line.

	if (b->get_system_time() <= 0)

	    initial_mass = total_mass(b);

	else {

	    initial_mass = 1;
	    if (verbose)
		cerr << "initial_mass not known; adopting default = "
		     << initial_mass << endl;
	}

	putrq(b->get_log_story(), "initial_mass", initial_mass);
    }

    return initial_mass;
}

real get_initial_virial_radius(hdyn* b,
			       bool verbose,		// default = false
			       bool r_virial_set,	// default = false
			       real input_r_virial)	// default = 0
{
    real r_virial = getrq(b->get_log_story(), "initial_rvirial");

    if (r_virial > 0) {

	// The input snapshot has the information needed.  Snapshot
	// data OVERRIDES any command-line entry.

	if (verbose)
	    cerr << "Read initial r_virial = " << r_virial
		 << " from input snapshot" << endl;

    } else {

	// See if r_virial was set from the command-line, or use default.

	if (r_virial_set) {

	    if (verbose) {

		cerr << "Setting initial r_virial = " << r_virial
		     << " from command-line input" << endl;

		if (b->get_system_time() > 0)
		    cerr << endl
			 << "Warning: setting initial r_virial with time > 0"
			 << endl;
	    }

	    r_virial = input_r_virial;

	} else {

	    // r_virial is not set in the snapshot or command line.
	    // Assume r_virial = 1.  Probably better to compute
	    // r_virial directly in this case.  However, the "mk..."
	    // functions and "scale" already set this quantity, so there
	    // seems to be no strong reason to add an extra energy
	    // computation to the code...(?)

	    r_virial = 1;

#if 0
	    // Alternative: get r_virial from input data.

	    real kinetic, potential;
	    get_top_level_energies(b, 0.0, potential, kinetic);

	    r_virial = -0.5 * pow(total_mass(b), 2) / potential_energy;
#endif

	    if (verbose)
		cerr << "Adopting default initial r_virial = "
		     << r_virial << endl;
	}

	// Write the virial radius to the log story if t <= 0.

	if (b->get_system_time() <= 0)
	    putrq(b->get_log_story(), "initial_rvirial", r_virial);
    }

    return r_virial;
}

real get_initial_jacobi_radius(hdyn* b,
			       real r_virial,
			       bool verbose,		// default = false
			       bool r_jacobi_set,	// default = false
			       real input_r_jacobi)	// default = 0
{
    real initial_r_jacobi = -1;

    // Determine the initial Jacobi radius of the system.  This radius
    // may be specified on the command line, or it may be derived from
    // the information in the input snapshot.  Unlike r_virial, the
    // command-line takes precedence over initial data found in the
    // input snapshot, although it is superceded by any "kira" data
    // found in that file.  Scaling is such that, if a non-kira
    // version of initial_r_jacobi is found in the input snapshot, the
    // command-line value will scale that number.  A kira version is
    // used as is.  Otherwise, the command-line value scales the virial
    // radius (hence the two tests of r_jacobi_set below).

    // See if a kira_initial_jacobi_radius exists in the input file.
    // If it does, it TAKES PRECEDENCE over any other setting.

    if (find_qmatch(b->get_log_story(), "kira_initial_jacobi_radius")) {

	real kira_initial_jacobi_radius
	    = getrq(b->get_log_story(), "kira_initial_jacobi_radius");

	if (kira_initial_jacobi_radius > 0) {

	    initial_r_jacobi = kira_initial_jacobi_radius;

	    if (verbose) {
		cerr << "Using initial Jacobi radius ";
		PRL(kira_initial_jacobi_radius);
		cerr << "    from input snapshot" << endl;

		if (r_jacobi_set)
		    cerr << "Ignoring \"-J " << input_r_jacobi
			 << "\" found on command line"
			 << endl;
	    }

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "kira_initial_jacobi_radius "
		     << "from input snapshot"
		     << endl;
	    else
		err_exit(
	     "Error reading kira_initial_jacobi_radius from input snapshot");

	}
    }

    if (initial_r_jacobi < 0) {

	// See if we can determine the system's initial Jacobi radius from
	// the input snapshot.

	if (find_qmatch(b->get_log_story(), "initial_rtidal_over_rvirial")) {

	    // The input snapshot contains an initial "tidal" radius.  Most
	    // likely the initial model is a King profile and this is the
	    // King radius (misnamed for historical reasons -- it may or
	    // may not have anything to do with a tidal cutoff).

	    real r_jacobi_over_r_virial = getrq(b->get_log_story(),
						"initial_rtidal_over_rvirial");

	    if (r_jacobi_over_r_virial > 0) {

		initial_r_jacobi = r_jacobi_over_r_virial * r_virial;

		if (verbose) {
		    cerr << "Got r_jacobi_over_r_virial = "
			 << r_jacobi_over_r_virial
			 << " from input snapshot"
			 << endl
			 << "    setting ";
		    PRL(initial_r_jacobi);
		}

	    } else {

		if (verbose)
		    cerr << "Warning: error reading "
			 << "initial_rtidal_over_rvirial from input snapshot"
			 << endl;
		else
		    err_exit(
	      "Error reading initial_rtidal_over_rvirial from input snapshot");
	    }
	}

	// See if there was any information on the command line, and resolve
	// conflicts.

	if (r_jacobi_set) {

	    // A "-J" entry was specified on the command line.

	    if (initial_r_jacobi > 0) {

		// The Jacobi radius has been specified both in the input file
		// and on the command line.  If the model is an anisotropic
		// King model, we must use the snapshot data and ignore the
		// command line.  Otherwise, use the command line as a scaling
		// for the "known" value.

		if (find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {

		    if (verbose)
			cerr << "Ignoring command-line Jacobi radius (-J "
			     << input_r_jacobi << ") for anisotropic King model"
			     << endl;

		} else {

		    if (verbose)
			cerr << "Command-line Jacobi radius (-J "
			     << input_r_jacobi << ") scales initial value "
			     << initial_r_jacobi
			     << endl
			     << "    from input snapshot"
			     << endl;

		    initial_r_jacobi *= input_r_jacobi;

		}

	    } else {

		// Use the command-line data to scale the virial radius.

		if (verbose)
		    cerr << "Command-line Jacobi radius (-J "
			 << input_r_jacobi << ") scales initial virial radius "
			 << r_virial << endl;

		initial_r_jacobi = input_r_jacobi * r_virial;

	    }
	}

	// Save the Jacobi radius for restart.

	if (initial_r_jacobi > 0)
	    putrq(b->get_log_story(), "kira_initial_jacobi_radius",
		  initial_r_jacobi);

    }

    return initial_r_jacobi;
}

#define MATCH_TOL 0.001

local bool matches(real r, real v)
{
    return (abs(r/v - 1) < MATCH_TOL);
}

static char* tidal_type[5] = {"none",
			      "point-mass",
			      "isothermal",
			      "disk",
			      "custom"};

// Express Galactic parameters in "Stellar" units (see Mihalas 1968):
//
//	G		=  1
//	length unit	=  1 pc
//	velocity unit	=  1 km/s
//
// ==>	time unit 	=  0.978 Myr
//	mass unit	=  232 Msun

#define OORT_A	(0.0144)	// km/s/pc
#define OORT_B	(-0.012)	// km/s/pc
#define RHO_G	(0.11/232)	// (232 Msun)/pc^3

#define Rsun_pc 2.255e-8	// R_sun/1 parsec = 6.960e+10/3.086e+18;

#define OORT_ALPHA_RATIO \
	((4*M_PI*RHO_G + 2*(OORT_A*OORT_A - OORT_B*OORT_B)) \
			 / (-4*OORT_A*(OORT_A-OORT_B)))

local void tidal_msg(int i, real alpha3_over_alpha1)
{
    cerr << "Forcing tidal_field_type = " << i
	 << " for anisotropic King model "
	 << endl
	 << "    with alpha3_over_alpha1 = "
	 << alpha3_over_alpha1
	 << endl;
}

void set_tidal_params(hdyn* b,
		      bool verbose,
		      real initial_r_jacobi,
		      real initial_mass,
		      int& tidal_field_type)
{
    // Set the tidal parameters for the system.
    //
    // Procedure:
    //
    //		default tidal_field_type is 1 (point-mass)
    //		initial_mass and initial_r_jacobi set the value of alpha1
    //		tidal_field_type then sets the value of alpha3
    //
    // NOTE that tidal_field_type = 3, the field is intended to
    // model the Galactic field in the solar neighborhood, for
    // which alpha1 and alpha3 are actually determined directly by
    // the local Oort constants.  The resultant field may thus *not*
    // be consistent with the choice of radius used later in
    // establishing physical scales.

    // Tidal field type or anisotropic King model initial info in
    // input snapshot overrides command-line version, if any.

    int kira_tidal_field_type = -1;

    if (find_qmatch(b->get_log_story(), "kira_tidal_field_type")) {

	kira_tidal_field_type
	    = getiq(b->get_log_story(), "kira_tidal_field_type");

	if (kira_tidal_field_type > 0) {

	    if (verbose) {
		cerr << "Using tidal_field_type = " << kira_tidal_field_type
		     << " (" << tidal_type[kira_tidal_field_type]
		     << ") from input snapshot" << endl;

		if (tidal_field_type > 0)
		    cerr << "Ignoring \"-F " << tidal_field_type
			 << "\" found on command line"
			 << endl;
	    }

	    tidal_field_type = kira_tidal_field_type;

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "kira_tidal_field_type "
		     << "from input snapshot"
		     << endl;
	    else
		err_exit(
	     "Error reading kira_tidal_field_type from input snapshot");

	}

    } else if (find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {

	// Try to infer tidal_field_type from alpha3_over_alpha1.

	real alpha3_over_alpha1
	    = getrq(b->get_log_story(), "alpha3_over_alpha1");

	if (alpha3_over_alpha1 < 0
	    && alpha3_over_alpha1 > -VERY_LARGE_NUMBER) {

	    // See if the input ratio matches any of the standard types.
	    // If it does, use that type; if not, flag a warning and use
	    // the command-line type or the default.

	    if (matches(alpha3_over_alpha1, -1.0/3)) {
		if (tidal_field_type != 1) {
		    tidal_field_type = 1;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    } else if (matches(alpha3_over_alpha1, -1.0/2)) {
		if (tidal_field_type != 2) {
		    tidal_field_type = 2;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    }  else if (matches(alpha3_over_alpha1, OORT_ALPHA_RATIO)) {
		if (tidal_field_type != 3) {
		    tidal_field_type = 3;
		    if (verbose)
			tidal_msg(tidal_field_type, alpha3_over_alpha1);
		}
	    } else {
		cerr << "Warning: snapshot alpha3_over_alpha1 = "
		     << alpha3_over_alpha1
		     << " does not match any standard value"
		     << endl;
	    }

	} else {

	    if (verbose)
		cerr << "Warning: error reading "
		     << "alpha3_over_alpha1 "
		     << "from input snapshot"
		     << endl;

	}
    }

    if (tidal_field_type > 4)
	err_exit("Unknown tidal field type");

    if (verbose) {
	cerr << "Using ";

	if (tidal_field_type <= 1) {

	    cerr << tidal_type[1] << " ";
	    if (tidal_field_type <= 0) cerr << "(default) ";

	} else
	    cerr << tidal_type[tidal_field_type] << " ";

        cerr << "tidal field; "
	     << "initial Jacobi radius = " << initial_r_jacobi
	     << endl;
    }

    if (tidal_field_type <= 0) tidal_field_type = 1;

    // Set up the dynamical tidal parameters.

    real alpha1 = -initial_mass / pow(initial_r_jacobi, 3.0);	// (definition)
    real alpha3, omega;

    // Note that we don't use alpha3_over_alpha1, even if it is stored
    // in the snapshot.  We use the "standard" ratios, or rederive the
    // ratio from the stored Oort constants.

    if (tidal_field_type == 1) {

	// Point-mass field (see Binney & Tremaine, p. 452).

	alpha3 = -alpha1/3;
	omega = sqrt(alpha3);

    } else if (tidal_field_type == 2) {

	// Isothermal halo.

	alpha3 = -alpha1/2;
	omega = sqrt(alpha3);
	
    } else if (tidal_field_type == 3) {

	// Disk field.  Use the local Oort constants (defined above).
	// Conversion between Oort constants and alpha1/3 is taken
	// from Heggie & Ramamani (MNRAS 272, 317, 1995):
	//
	//	alpha1 = -4 A (A - B)
	//	alpha3 = 4 PI G RHO_G + 2 (A^2 - B^2)
	//
	// and recall that G = 1 by our above choice of units.

	alpha3 = alpha1 * OORT_ALPHA_RATIO;

	// Get omega from the definition of A and B:
	//
	//	A =  (1/2) (Vc/R - dVc/dR)
	//	B = -(1/2) (Vc/R + dVc/dR)
	//
	// ==>	omega = Vc/R = A - B,  dVc/dR = A + B
	//
	// so	alpha1 = -2 omega (omega - dVc/dR)
	//	       = -2 omega^2 (1 - (A + B)/(A - B))
	//	       =  4 omega^2 B / (A - B)

	omega = sqrt(0.25 * alpha1 * (OORT_A/OORT_B - 1));

    } else if (tidal_field_type == 4) {

	// Custom tidal field.  Require that the original physical
	// parameters be saved in the input snapshot.

	real Oort_A = 0;
	real Oort_B = 0;
	real rho_G = 0;

	if (find_qmatch(b->get_log_story(), "Oort_A_constant")) 
	    Oort_A = getrq(b->get_log_story(), "Oort_A_constant");
	else
	    err_exit("Oort A constant not specified");

	if (find_qmatch(b->get_log_story(), "Oort_B_constant")) 
	    Oort_B = getrq(b->get_log_story(), "Oort_B_constant");	
	else
	    err_exit("Oort B constant not specified");

	if (find_qmatch(b->get_log_story(), "local_mass_density")) 
	    rho_G = getrq(b->get_log_story(), "local_mass_density");
	else
	    err_exit("rho_G not specified");

	cerr << "Create custom tidal field from " << endl;
	PRI(4);PRC(Oort_A);PRC(Oort_B);PRL(rho_G);

	if (Oort_A != 0 && Oort_B != 0 && rho_G != 0) {

            // alpha1 = -4*Oort_A*(Oort_A-Oort_B);	// no!

	    real alpha3_over_alpha1 = 
		(4*M_PI*rho_G + 2*(pow(Oort_A, 2) - pow(Oort_B, 2)))
		    / (-4*Oort_A*(Oort_A - Oort_B));

	    alpha3 = alpha1 * alpha3_over_alpha1;
	    omega = sqrt(0.25*alpha1 * (Oort_A/Oort_B - 1));
	}
	else
	  err_exit("Insufficient information to reconstruct tidal field");
    }

    if (verbose) {
	PRI(4); PRC(alpha1); PRC(alpha3); PRL(omega);
    }

    b->set_tidal_field(tidal_field_type);
    b->set_omega(omega);
    b->set_alpha(alpha1, alpha3);

    // Save the field type information in the root log story, if necessary.

    if (kira_tidal_field_type <= 0)
	putiq(b->get_log_story(), "kira_tidal_field_type",
	      tidal_field_type);

}

void test_tidal_params(hdyn* b,
		       bool verbose,
		       real initial_r_jacobi,
		       real initial_r_virial,
		       real initial_mass)
{
    cerr << endl << "Comparing initial parameters for disk tidal field:"
	 << endl;

    fprintf(stderr, "    model r_virial = %.3f, r_jacobi = %.3f",
	    initial_r_virial, initial_r_jacobi);
    fprintf(stderr, ", ratio = %.3f\n", initial_r_jacobi/initial_r_virial);

    // Convert initial mass and virial radius using conversion factors.
    // Compute Jacobi radius from Oort constants.

    initial_mass
	= b->get_starbase()->conv_m_dyn_to_star(initial_mass);	    // Msun
    initial_r_virial
	= b->get_starbase()->conv_r_dyn_to_star(initial_r_virial)
	    * Rsun_pc;						    // pc

    initial_r_jacobi = pow((initial_mass/232)
			    / (4*OORT_A*(OORT_A-OORT_B)), 1.0/3);   // pc

    fprintf(stderr, "    real  r_virial = %.3f, r_jacobi = %.3f",
	    initial_r_virial, initial_r_jacobi);
    fprintf(stderr, ", ratio = %.3f\n", initial_r_jacobi/initial_r_virial);
}
