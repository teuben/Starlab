
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add a mass spectrum to an input snapshot.  Existing node masses
//// are overwritten.  Dyn version is copied from the node version.
////
//// Usage: makemass [OPTIONS] < input > output
////
//// Options:   
////           -C        force output data to be in 'col' format [no]
////           -e/E/x/X  exponent [-2.35 (Salpeter)]
////           -F/f      mass function option: 1) Power-law [default]
////                                           2) Miller & Scalo
////                                           3) Scalo
////                                           4) Kroupa
////                                           5) GdeMarchi
////                     Option -F requires one of the following strings:
////                     (Power_Law, Miller_Scalo, Scalo, Kroupa, GdeMarchi).
////                     Option -f requires the appropriate integer.
////           -i        (re)number stellar index from highest to lowest mass.
////           -l/L      lower mass limit [1]
////           -m/M      scale to specified total mass [don't scale]
////           -u/U      upper mass limit [1]
////           -s        random seed
////
//// Written by Steve McMillan and Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//++ Note: The conversion factor for scaling between dynamical and stellar masss
//++        is properly set in the output snapshot.

//		Steve McMillan, July 1996
//		Simon Portegies Zwart, Tokyo, December 1997

// need to repair bugs, see PJT comments

// Need to re-merge node and dyn versions!!  Steve, 7/04

#include "dyn.h"

#define  MAXIMUM_ZAMS_MASS 100
#define  SEED_STRING_LENGTH  256
static char  tmp_string[SEED_STRING_LENGTH];

typedef  struct
{
    dyn* str;
    real  mass;
} nm_pair, *nm_pair_ptr;

// Probably would be much better if these functions were contained in the
// node version and made global.  (Steve, 6/03)

//-----------------------------------------------------------------------------
//  compare_mass  --  compare the masses of two particles
//-----------------------------------------------------------------------------

local int compare_mass(const void * pi, const void * pj)
{
    if (((nm_pair_ptr) pi)->mass < ((nm_pair_ptr) pj)->mass)
        return(1);
    else if (((nm_pair_ptr)pi)->mass > ((nm_pair_ptr)pj)->mass)
        return(-1);
    else
        return(0);
}

local void makemass(dyn* b, mass_function mf,
		    real m_lower, real m_upper,
		    real exponent, real total_mass, bool renumber_stars) {

    real m, m_sum = 0;
    int n=0;
    for_all_daughters(dyn, b, bi) {
	m = get_random_stellar_mass(m_lower, m_upper, mf, exponent);
	n++;
	bi->set_mass(m);
	m_sum += bi->get_mass();
    }

    b->set_mass(m_sum);

    // Renumber the stars in order of mass.
    // Highest mass gets smallest number (strange choise, but).
    if (renumber_stars) {
      int istart = 1;
      bool M_flag = true;
      renumber(b, istart, M_flag);
    }

    if (total_mass > 0) {
	real m_factor = total_mass/m_sum;
	for_all_daughters(dyn, b, bi)
	    bi->set_mass(m_factor*bi->get_mass());
	b->set_mass(total_mass);
    }

    if(m_lower>=m_upper)
	sprintf(tmp_string,
		"       %s mass function, total mass = %8.2f",
		type_string(mf), m_sum);
//		type_string(Equal_Mass), m_sum);
    else
	sprintf(tmp_string,
		"       %s mass function, total mass = %8.2f Solar",
		type_string(mf), m_sum);
    b->log_comment(tmp_string);
    cerr << "mass function is " << type_string(mf) << ", ";
    PRC(m_lower); PRL(m_upper);

}

int main(int argc, char ** argv)
{
    bool C_flag = false;
    bool F_flag   = false;                      // input mf via string
    mass_function mf = mf_Power_Law;            // default = Power-law
    char *mfc = new char[64];
    real m_lower  = 1, m_upper = 1;		// default = equal masses
    bool x_flag   = false;
    real exponent = -2.35;			// Salpeter mass function
                                                // if no exponent is given
                                                // Miller & Scalo is used
    real m_total = -1;				// don't rescale

    bool renumber_stars = false;                // renumber stellar index
                                                // from high to low mass

    int random_seed = 0;
    char seedlog[64];

    bool lower_set = false, upper_set = false;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "CE:e:iF:f:L:l:M:m:s:U:u:X:x:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'C': C_flag = true;
		      break;
	    case 'F': F_flag = true;
		      strncpy(mfc, poptarg, 63);
	              break;
	    case 'f': mf = (mass_function)atoi(poptarg);
	              break;
	    case 'i': renumber_stars = true;
	              break;
	    case 'E':
	    case 'e':
	    case 'X':
	    case 'x': x_flag = true;
	              exponent = atof(poptarg);
		      break;
	    case 'L':
	    case 'l': m_lower = atof(poptarg);
		      lower_set = true;
		      break;
	    case 'M':
	    case 'm': m_total = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
	    case 'U':
	    case 'u': m_upper = atof(poptarg);
		      upper_set = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}

    if (lower_set && !upper_set) m_upper = m_lower;
    if (!lower_set && upper_set) m_lower = m_upper;

    if (m_lower <= 0 ||
	m_upper <= 0 ||
	m_lower > m_upper)
	err_exit("makemass: Illegal mass limits");

    if (F_flag)
      mf = extract_mass_function_type_string(mfc);
    delete [] mfc;

    dyn *b;
    b = get_dyn();

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    PRC(m_lower); PRL(m_upper);
    makemass(b, mf, m_lower, m_upper, exponent, m_total, renumber_stars);	
    real initial_mass = getrq(b->get_log_story(), "initial_mass");

    if (initial_mass > -VERY_LARGE_NUMBER)
        putrq(b->get_log_story(), "initial_mass", b->get_mass(),
	      HIGH_PRECISION);

    // Add stellar conversion data, so add_star will work.

    real m_sum = b->get_mass();
    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
    if(old_mtot<=0) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_sum,
							 old_r_vir,
							 old_t_vir);
    }

    if (C_flag) b->set_col_output(true);
    put_dyn(b);
    rmtree(b);
}
