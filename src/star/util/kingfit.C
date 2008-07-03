//  kingfit.C: Print out various diagnostic statistics on a system.

#include "stdfunc.h"
#include "sstar_to_dyn.h"
#ifndef TOOLBOX

local real kingfit(const real X, int nfit, const real Xarray[],
	           const real Yarray[]) {

  real parameter = -1;
  if (X>=Xarray[1] || X<=Xarray[nfit-1])
    return parameter;
  
  if (X>=Xarray[1])
    return Yarray[0];
  else if(X<=Xarray[nfit-1])
    return Yarray[nfit-1];
  else
    for (int i=1; i<nfit; i++) 
      if (X>=Xarray[i])
	return lineair_interpolation(X,
				     Xarray[i-1], Xarray[i],
				     Yarray[i-1], Yarray[i]);
  return parameter;
}

// Fits to King models.
// A total of 102400 stars (100 times 1024) is used per fitted parameter.
void print_fitted_king_model(const real fractional_radius,
			     const radius_fraction which_fraction) {
			     
  int nfit = 16;
  real WoKing[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  real CKing[] = {0., 0.2956, 0.5054, 0.6739, 0.8405, 1.03059,
		  1.2560, 1.5284, 1.8341, 2.1202, 2.3516, 2.5495,
		  2.7396, 2.9346, 3.1414, 3.3570};
  real RK_Rvir[] = {2.0, 1.29917, 0.87470, .66626, 0.52325, 0.40736,
                    0.30368, 0.20692, 0.12217, 0.06357, 0.03311,
		    0.01870, 0.01132, 0.00710, 0.00450, 0.00282}; 

  real Rc_Rvir[] = {1., 0.5078, 0.4711, 0.4274, 0.3718, 0.3136, 0.2479,
		    0.1744, 0.1087, 0.05839, 0.03240, 0.01895, 0.01335,
		    0.009796, 0.007413, 0.006709};
  
  real proj_r1r9[] = {0, 0.2370, 0.2282, 0.2146, 0.1945, 0.1716, 0.1397,
		      0.1050, 0.07199, 0.05437, 0.05082, 0.05855, 0.06782,
		      0.07320, 0.07627, 0.07587};
  
  real r1r9[] = {0, 0.3234, 0.3052, 0.2825, 0.2520, 0.2160, 0.1711,
		 0.1224, 0.08080, 0.06040, 0.05631, 0.06693, 0.08000,
		 0.08907, 0.09046, 0.08942};

  real Wo = -1;
  real C = -1;
  switch(which_fraction) {
     case rcore_rvirial:
          Wo = kingfit(fractional_radius, nfit, Rc_Rvir, WoKing);
	  C = kingfit(fractional_radius, nfit, Rc_Rvir, CKing);
          cerr <<"     (non projected rcore/rvir)\n";
       break;
     case r10_r90_light:
          Wo = kingfit(fractional_radius, nfit, r1r9, WoKing);
	  C = kingfit(fractional_radius, nfit, r1r9, CKing);
          cerr <<"     (non projected r at 10% / 90% light )\n";
       break;
     case proj_r10_r90_light:
          Wo = kingfit(fractional_radius, nfit, proj_r1r9, WoKing);
	  C = kingfit(fractional_radius, nfit, proj_r1r9, CKing);
          cerr <<"     (projected r at 10% / 90% light )\n";
       break;
     default:
       cerr << "    Unknown option in fitted_king_model()." << endl;
       return;
  }
  
  cerr << "    King model Wo = " << Wo
       << "   concentration c = " << C
       << endl;
  
}

#else

main(int argc, char **argv)
{
  real rc_rvir = 0;
  real r_lumi  = 0;
  int  nlagrad = 0;   // use the R10 over R90 % lagrangian radii.

  extern char *poptarg;
  int c;
  const char *param_string = "L:l:R:r:";

  if (nlagrad<0 || nlagrad>2) {
    cerr << " kingfit -n " << nlagrad << endl;
    cerr << "please specify n between 0 and 2." << endl;
    exit (1);
  }
    
    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'R':
	    case 'r': rc_rvir = atof(poptarg);
		      break;
	    case 'L':
	    case 'l': r_lumi = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

   if (rc_rvir)
     print_fitted_king_model(rc_rvir, rcore_rvirial);
   else
     print_fitted_king_model(r_lumi, proj_r10_r90_light);

}
#endif
