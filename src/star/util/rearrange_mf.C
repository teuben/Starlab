/*
 *  red_stellar system.C: reduce some useful information from the
 *  			 stellar system
 *.............................................................................
 *    version 1:  Jan 1997   Simon Portegies Zwart   email: spz@astro.uva.nl
 *.............................................................................
 *  non-local function: 
 *  all functions are local
 *.............................................................................
 *  Input options:
 *     -B:	Binning option of the Lagrangian radii.
 *     -C:	Cut off creterium, mass, luminosity, number of stars.
 *     -l:      Lower luminosity limit (see below).
 *     -n:	Number of Lagrangian radii bins.
 *		Might for some choises be forced.
 *     -O:      Output option, what information should be studied.
 *     -o:	Output the story with data written in root.
 *     -S:      Sort option, on what parameter must information be sorted.
 *.............................................................................
 *  Conserning the Lower luminosity limit option (-l, see above).
 *    Currently this option is only applied to the sorting functions.
 *    This means that the binning (using Lagrangian radii) is applied
 *    on all cluster members.
 *    whether or not this is realistic depends on the application.
 *    Consequently, the luminosity cut-off does not affect the
 *    Lagrangian radii.
 *.............................................................................
 */
#include "stardyn_util.h"

#ifdef TOOLBOX

#define SCRATCH_PAD_LINE_LENGTH 255

//-----------------------------------------------------------------------------
//  compare_parameters  --  compare parameters of two particles
//-----------------------------------------------------------------------------

local int compare_parameters(const void* pi, const void* pj)
{
  if ((*(const real *)pi) > (*(const real *)pj))
    return(1);
  else if ((*(const real *)pi) < (*(const real *)pj))
    return(-1);
  else
    return(0);
}

local void rearrange_mf(dyn* b) {

  int  n = b->n_daughters();

  int*  name = new int[n];
  real* mass = new real[n];

  int i=0;

  for_all_daughters(dyn, b, bi) {	

    name[i] = bi->get_index();
    mass[i] = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
    i++;
  }
  n=i-1;
  qsort((real*)mass, (size_t)n, sizeof(real), compare_parameters);

  for (i=0; i<n; i++)
     cerr << mass[i] << endl;
}

/*-----------------------------------------------------------------------------
 *  main  --  driver to use  compute_mass_radii() as a tool
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
  char  *comment;
  bool verbatim = false;
  bool  c_flag = FALSE;      /* if TRUE: a comment given on command line   */
    
  extern char *poptarg;
  int c;
  char* param_string = "c:";

  dyn *b;

  cerr.precision(STD_PRECISION);
  while (b = get_dyn(cin)) {
    
    if (c_flag == TRUE)
      b->log_comment(comment);

    b->log_history(argc, argv);

    rearrange_mf(b);

    delete b;
  }
}

#endif
