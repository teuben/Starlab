////
////  OortAB.C: Utility program to calculate the Oort constants
////
////     input: -R # maximum distance from teh Galactic center [NA]
////            -r # Minimum distance from teh Galactic center [NA]
////            -N # number of output steps [1]
////                 output is equally devided in (rmax-rmin)/N
////            -M # Mass of the star cluster.
////
////            Simon Portieges Zwart, March 2000, 
////            Developed in Boston Univeristy: spz@komodo.bu.edu

#include "stdfunc.h"
#ifndef TOOLBOX

#else

//local real Oort_AB(const real X, int nfit, const real Xarray[],
//	           const real Yarray[]) {
//}

local real galaxy_mass_within_radius(real r) {

  real m_galaxy = 4.25E+6 * pow(r, 1.2); // Msun

  return m_galaxy;
}


local void print_Oort_constants(real r_gc, real m_cluster) {

  real m_galaxy = galaxy_mass_within_radius(r_gc);

  real r_tide = pow(m_cluster/(2*m_galaxy), ONE_THIRD) * r_gc;  // pc

  real v_c = 0.0657 * sqrt((m_galaxy+m_cluster)/r_gc);  // km/s
  real dvc_dr = 13.7 / pow(r_gc, 0.9);                  // km/s

  real Oort_A = (v_c/r_gc - dvc_dr)/2;                  // km/s/pc
  real Oort_B = -(v_c/r_gc + dvc_dr)/2;                 // km/s/pc

  real dr = 0.5;     // [pc]
  real rho_G = (galaxy_mass_within_radius(r_gc-dr)
		- galaxy_mass_within_radius(r_gc+dr))
             / ((4*PI/3) * (pow(r_gc-dr, 3) - pow(r_gc+dr, 3)));

  cerr << "    r_gc= "<<r_gc<<" m_cluster= "<<m_cluster
       << " m_galaxy= "<<m_galaxy
       << " A= " <<Oort_A<< " B= " <<Oort_B<<" rho_G= " << rho_G << endl;
  cerr << " r_tide= "<<r_tide<<" v_c= "<<v_c<<" dvc_dr= "<< dvc_dr << endl;
}


main(int argc, char **argv) {

  bool m_flag = false;
  bool n_flag = false;
  bool r_flag = false;

  real r_min = 30;     // [pc]
  real r_max = 30;     // [pc]
  real m_cluster = 20000;  // [Msun]
  int n = 1;

  extern char *poptarg;
  int c;
  char* param_string = "M:m:N:n:R:r:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'R': r_max = atof(poptarg);
		      r_flag = true;
		      break;
	    case 'r': r_min = atof(poptarg);
		      r_flag = true;
		      break;
	    case 'N': 
	    case 'n': n = atoi(poptarg);
		      n_flag = true;
		      break;
	    case 'M':
	    case 'm': m_cluster = atof(poptarg);
		      m_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

  if (r_min>r_max) {
    r_max = r_min;
    n=1;
  }

  cerr << "Calculate Oort constants. " << endl; 
   if (r_flag && m_flag) {

     real r_gc;
     real dr = (r_max-r_min)/n;

     for(int i=0; i<n; i++) { 
       r_gc = r_min + dr*i;
       print_Oort_constants(r_gc, m_cluster);
     }
   }
}
#endif
