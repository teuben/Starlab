
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

//// mkmass:  Add a mass spectrum to an input snapshot.  Existing
////          node masses are overwritten.
////
//// Options:   -e/E/x/X  exponent [-2.35 (Salpeter)]
////            -F/f      mass function option: 1) Power-law [default]
////                                            2) Miller & Scalo
////                                            3) Scalo
////                                            4) Kroupa
////            Option -F requires one of the following strings:
////                      (Power_Law, Miller_Scalo, Scalo, Kroupa)
////                   -f requires the appropriate interger.
////            -i        (re)number stellar index from highest to lowest mass.
////            -l/L      lower mass limit [1]
////            -m/M      scale to specified total mass [don't scale]
////            -u/U      upper mass limit [1]
////            -s        random seed
////
//++ Note: The conversion factor for scaling between dynamical and stellar masss
//++        is properly set in the output snapshot.

//		Steve McMillan, July 1996
//		Simon Portegies Zwart, Tokyo, December 1997

// need to repair bugs, see PJT comments

#include "node.h"

#define  MAXIMUM_ZAMS_MASS 100
#define  SEED_STRING_LENGTH  256
static char  tmp_string[SEED_STRING_LENGTH];

#ifndef TOOLBOX


//see: Eggleton, Fitchet & Tout 1989, ApJ 347, 998
local real mf_Miller_Scalo(real m_lower, real m_upper) {

    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.19*rnd
	    / (pow(1-rnd, 0.75) + 0.032*pow(1-rnd, 0.25));
    }
    while(m_lower>m || m>m_upper); 
    return m;
}

// see: de la Fuente Marcos, Aarseth, Kiseleva, Eggleton 
// 1997, in Docobo, Elipe, McAlister (eds.), Visual Double
// Stars: Formation, dynamics and Evolutionary Tracks, KAP: ASSL
// Series vol. 223,  165
local real mf_Scalo(real m_lower, real m_upper) {
    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.3*rnd
	    / pow(1-rnd, 0.55);
    }
    while(m_lower>m || m>m_upper); 
    return m;
}

local real Kroupa_Tout_Gilmore(real m_lower, real m_upper) {

    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.08 + (0.19*pow(rnd, 1.55) + 0.05*pow(rnd, 0.6))
	    /        pow(1-rnd, 0.58);
    }
    while(m_lower>m || m>m_upper); 
    return m;
}

real get_random_stellar_mass(real m_lower, real m_upper, 
			     mass_function mf, real exponent) { 

    real m;
    switch(mf) {
       case Equal_Mass:
	  if (m_lower==0) {
	     cerr << "get_random_stellar_mass:"<<endl;
	     cerr << "unambiguous choise of Equal_Mass."<<endl;
	     cerr << "Use -m option to set fixed stellar mass."<<endl;
	     exit(1);
	  }
	  m =  m_lower;
	       break;
       case mf_Power_Law:
	  m =  general_power_law(m_lower, m_upper, exponent);
	       break;
       case Miller_Scalo: 
	   m = mf_Miller_Scalo(m_lower, m_upper);
	       break;
       case Scalo:     
	  m = mf_Scalo(m_lower, m_upper);
	      break;
       case Kroupa: // Kroupa, Tout & Gilmore 1993, MNRAS 262, 545
	  m = Kroupa_Tout_Gilmore(m_lower, m_upper);
	      break;
       default:
          cerr << "WARNING: \n"
	       << "        real get_random_stellar_mass:\n"
	       << "        Mass function parameters not properly defined.\n";
	  exit(1);
    }

    return m;
}

char* type_string(mass_function mf) {
    
    local char  mf_name[SEED_STRING_LENGTH];	// permanent
    switch(mf) {

       case Equal_Mass:
            sprintf(mf_name, "Equal_Mass"); 
	    break;	    
       case mf_Power_Law:
            sprintf(mf_name, "Power_Law"); 
	    break;	    
       case Miller_Scalo:
            sprintf(mf_name, "Miller_Scalo"); 
	    break;	    
       case Scalo:
            sprintf(mf_name, "Scalo"); 
	    break;	    
       case Kroupa:
            sprintf(mf_name, "Kroupa"); 
	    break;
       default:
            sprintf(mf_name, "Unknown"); 
	    break;
    }
    return mf_name;
}

mass_function extract_mass_function_type_string(char* type_string) {

     mass_function type = Unknown_MF;

     if (!strcmp(type_string, "Equal_Mass"))
        type = Equal_Mass;
     else if (!strcmp(type_string, "Power_Law")) 
        type = mf_Power_Law;
     else if (!strcmp(type_string, "Miller_Scalo")) 
        type = Miller_Scalo;
     else if (!strcmp(type_string, "Scalo"))
        type = Scalo;
     else if (!strcmp(type_string, "Kroupa")) 
        type = Kroupa;
     else if (!strcmp(type_string, "Unknown")) 
        type = Unknown_MF;
     else {
	 cerr << "No proper mass function indicated in mkmass."<<endl;
	 exit(-1);
     }

     return type;
   }

real general_power_law(real lowerl, real upperl, real exponent) {

    real x, norm;

    // Normalizing factor.
    if (exponent == -1)
	norm = log(upperl/lowerl);
    else
	norm = pow(upperl/lowerl, 1+exponent) - 1;

    if (exponent == -1)
	x = lowerl*exp(norm*randinter(0,1));
    else
	x = lowerl*pow(norm*randinter(0,1) + 1, 1/(1+exponent));

    return x;
}

#else

typedef  struct
{
    node* str;
    real  mass;
} nm_pair, *nm_pair_ptr;

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

local void mkmass(node* b, mass_function mf,
		  real m_lower, real m_upper, 
		  real exponent, real total_mass, bool renumber) {
  
    real m, m_sum = 0;
    int n=0;
    for_all_daughters(node, b, bi) {
	m = get_random_stellar_mass(m_lower, m_upper, mf, exponent);
	n++;
	bi->set_mass(m);
	m_sum += bi->get_mass();
    }

    b->set_mass(m_sum);

    // Renumber the stars in order of mass.
    // Highest mass gets smallest number (strange choise, but).
    if (renumber) {
      nm_pair_ptr nm_table = new nm_pair[n];
      if (nm_table == NULL) {
	cerr << "mkmass: "
	     << "not enough memory left for nm_table\n";
	return;
      }

      int i=0;
      for_all_daughters(node, b, bi) {
	nm_table[i].str = bi;
	nm_table[i].mass = bi->get_mass();
	i++;
      }

      qsort((void *)nm_table, (size_t)n, sizeof(nm_pair), compare_mass);

      for (i=0; i<n; i++) {
	nm_table[i].str->set_index(i+1);  // Ok, lets number from 1 to n
      }
    }

    if (total_mass > 0) {
	real m_factor = total_mass/m_sum;
	for_all_daughters(node, b, bi)
	    bi->set_mass(m_factor*bi->get_mass());
	b->set_mass(total_mass);
    }

    if(m_lower>=m_upper)
	sprintf(tmp_string,
		"         %s mass function, total mass = %8.2f",
		type_string(Equal_Mass), m_sum); 
    else
	sprintf(tmp_string,
		"         %s mass function, total mass = %8.2f Solar",
		type_string(mf), m_sum); 
    b->log_comment(tmp_string);
    cerr << "Mass function is " << type_string(mf) <<endl;

}

void main(int argc, char ** argv)
{
    bool F_flag   = false;                        // Input mf via string
    mass_function mf = mf_Power_Law;             // Default = Power-law
    char *mfc = new char[64];
    real m_lower  = 1, m_upper = 1;		// Default = equal masses
    bool x_flag   = false;
    real exponent = -2.35;			// Salpeter mass function
                                                // If no exponent is given
                                                // Miller & Scalo is used
    real m_total = -1;				// Don't rescale

    bool renumber = false;                      // renumber stellar index
                                                // from high to low mass

    int random_seed = 0;
    char seedlog[64];

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "E:e:iF:f:L:l:M:m:s:U:u:X:x:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'F': F_flag = true;
		      mfc = poptarg;
	              break;
	    case 'f': mf = (mass_function)atoi(poptarg);
	              break;
	    case 'i': renumber = true;
	              break;
	    case 'E':
	    case 'e':
	    case 'X':
	    case 'x': x_flag = true;
	              exponent = atof(poptarg);
		      break;
	    case 'L':
	    case 'l': m_lower = atof(poptarg);
		      break;
	    case 'M':
	    case 'm': m_total = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
	    case 'U':
	    case 'u': m_upper = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}

    if (m_lower <= 0 ||
	m_upper <= 0 ||
	m_lower > m_upper)
	err_exit("mkmass: Illegal mass limits");

    if (F_flag)
      mf = extract_mass_function_type_string(mfc);
    delete mfc;

    node* b;
    b = get_node(cin);

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "         random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    mkmass(b, mf, m_lower, m_upper, exponent, m_total, renumber);	
    real initial_mass = getrq(b->get_dyn_story(), "initial_mass");

    if (initial_mass > -VERY_LARGE_NUMBER)
        putrq(b->get_dyn_story(), "initial_mass", m_total);

    real m_sum = b->get_mass();
    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
    if(old_mtot<=0) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_sum,
							 old_r_vir,
							 old_t_vir);
    }
    put_node(cout, *b);
    rmtree(b);

}
#endif

