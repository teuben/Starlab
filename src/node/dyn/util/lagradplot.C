
//// lagradplot:  compute and plot Lagrangian radii for input N-body
////              system(s), based on the geometric center.
////
//// Options:     -c    add a comment to the output snapshot [false]
////              -t    use ten-percentiles rather than quartiles [quartiles]

//.............................................................................
//    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
//.............................................................................
//  non-local function: 
//    plot_mass_radii
//.............................................................................
//     ....
//  ....
//.............................................................................
//  see also: pnode.c
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

#define  MAX_NUMBER_OF_COLUMNS   79
#define  PMR_LENGTH_PER_COLUMN   0.1
#define  PMRP_LENGTH_PER_COLUMN  0.05

typedef  struct
    {
    real  radius;
    real  mass;
    } rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)
    {
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return(1);
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return(-1);
    else
        return(0);
    }

//-----------------------------------------------------------------------------
//  plot_mass_radii  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

void  plot_mass_radii(dyn * b)
    {
    int  i;
    int  n;
    real  cumulative_mass;
    real  half_mass;
    real  first_quartile_mass;
    real  third_quartile_mass;
    bool  within_first_quartile = TRUE;
    bool  within_half_mass = TRUE;
    bool  within_third_quartile = TRUE;
    story *st;
    rm_pair_ptr  rm_table;
    
//  quick fix to determine  n  for a flat tree:

    dyn * bi;
    for (n = 0, bi = b->get_oldest_daughter(); bi != NULL;
         bi = bi->get_younger_sister())
        n++;    

    rm_table = new rm_pair[n];
    if (rm_table == NULL)
        {
        cerr << "plot_mass_radii: not enough memory left for the rm_table[]\n";
	exit(1);
	}

    for (i = 0, bi=b->get_oldest_daughter(); bi != NULL;
	 i++, bi=bi->get_younger_sister())
        {
	(rm_table + i)->radius = abs(bi->get_pos());
	(rm_table + i)->mass = bi->get_mass();
        }	

    qsort((void *)rm_table, (size_t)n, sizeof(rm_pair), compare_radii);

    half_mass = 0.5 * b->get_mass();
    first_quartile_mass = 0.25 * b->get_mass();
    third_quartile_mass = 0.75 * b->get_mass();

    st = b->get_log_story();

    cumulative_mass = 0.0;
    i = 0;

    while (cumulative_mass < first_quartile_mass)
	cumulative_mass += (rm_table + i++)->mass;

    first_quartile_mass = (rm_table + i-1)->radius;

    while (cumulative_mass < half_mass)
	cumulative_mass += (rm_table + i++)->mass;

    half_mass = (rm_table + i-1)->radius;

    while (cumulative_mass < third_quartile_mass)
	cumulative_mass += (rm_table + i++)->mass;

    third_quartile_mass = (rm_table + i-1)->radius;

    for (i= 0; i < MAX_NUMBER_OF_COLUMNS; i++)
        {
	if (within_first_quartile)
	    {
	    if (PMR_LENGTH_PER_COLUMN * i < first_quartile_mass)
	        printf(" ");
	    else
	        {
		printf("*");
		within_first_quartile = FALSE;
	        }
	    }
	else if (within_half_mass)
	    {
	    if (PMR_LENGTH_PER_COLUMN * i < half_mass)
	        printf(" ");
	    else
	        {
		printf("*");
		within_half_mass = FALSE;
	        }
	    }
	else if (within_third_quartile)
	    {
	    if (PMR_LENGTH_PER_COLUMN * i < third_quartile_mass)
	        printf(" ");
	    else
	        {
		printf("*");
		within_third_quartile = FALSE;
	        }
	    }
        }
    printf("\n");
    }

//-----------------------------------------------------------------------------
//  plot_mass_radii_in_percentages  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

void  plot_mass_radii_in_percentages(dyn * b)
    {
    int  i, k;
    int  n;
    real  cumulative_mass;
    real  mass_percent[9];
    bool  within_mass_percent[9];
    rm_pair_ptr  rm_table;
    dyn * bi;
    
//  quick fix to determine  n  for a flat tree:

    for (n = 0, bi = b->get_oldest_daughter(); bi != NULL;
         bi = bi->get_younger_sister())
        n++;    

    for (i = 0; i < 9; i++)
        within_mass_percent[i] = TRUE;

    rm_table = new rm_pair[n];
    if (rm_table == NULL)
        {
        cerr << "plot_mass_radii_in_percentages: "
	     << "not enough memory left for the rm_table[]\n";
	exit(1);
	}

    for (i = 0, bi=b->get_oldest_daughter(); bi != NULL;
	 i++, bi=bi->get_younger_sister())
        {
	(rm_table + i)->radius = abs(bi->get_pos());
	(rm_table + i)->mass = bi->get_mass();
        }	

    qsort((void *)rm_table, (size_t)n, sizeof(rm_pair), compare_radii);

    for (i = 0; i < 9; i++)
        mass_percent[i] = ((1 + i) / 10.0) * b->get_mass();

    for (i = k = 0, cumulative_mass = 0; k < 9; k++)
        {
        while (cumulative_mass < mass_percent[k])
	    cumulative_mass += (rm_table + i++)->mass;
        mass_percent[k] = (rm_table + i-1)->radius;
	}

    for (i= 0; i < MAX_NUMBER_OF_COLUMNS; i++)
        {
	if (within_mass_percent[0])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[0])
	        printf(" ");
	    else
	        {
		within_mass_percent[0] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[1])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[2])
		    {
		    printf("X");
		    within_mass_percent[1] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[1] = FALSE;
		    within_mass_percent[2] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[1])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[1])
	        printf(" ");
	    else
	        {
		within_mass_percent[1] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[2])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[3])
		    {
		    printf("X");
		    within_mass_percent[2] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[2] = FALSE;
		    within_mass_percent[3] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[2])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[2])
	        printf(" ");
	    else
	        {
		within_mass_percent[2] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[3])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[4])
		    {
		    printf("X");
		    within_mass_percent[3] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[3] = FALSE;
		    within_mass_percent[4] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[3])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[3])
	        printf(" ");
	    else
	        {
		within_mass_percent[3] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[4])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[5])
		    {
		    printf("X");
		    within_mass_percent[4] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[4] = FALSE;
		    within_mass_percent[5] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[4])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[4])
	        printf(" ");
	    else
	        {
		within_mass_percent[4] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[5])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[6])
		    {
		    printf("X");
		    within_mass_percent[5] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[5] = FALSE;
		    within_mass_percent[6] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[5])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[5])
	        printf(" ");
	    else
	        {
		within_mass_percent[5] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[6])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[7])
		    {
		    printf("X");
		    within_mass_percent[6] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[6] = FALSE;
		    within_mass_percent[7] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[6])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[6])
	        printf(" ");
	    else
	        {
		within_mass_percent[6] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[7])
		    {
		    printf("|");
		    }
		else if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[8])
		    {
		    printf("X");
		    within_mass_percent[7] = FALSE;
		    }
		else
		    {
		    printf("@");
		    within_mass_percent[7] = FALSE;
		    within_mass_percent[8] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[7])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[7])
	        printf(" ");
	    else
	        {
		within_mass_percent[7] = FALSE;
		if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[8])
		    {
		    printf("|");
		    }
		else
		    {
		    printf("X");
		    within_mass_percent[8] = FALSE;
		    }
	        }
	    }
	else if (within_mass_percent[8])
	    {
	    if (PMRP_LENGTH_PER_COLUMN * i < mass_percent[8])
	        printf(" ");
	    else
	        {
		within_mass_percent[8] = FALSE;
	        printf("|");
	        }
	    }
        }
    printf("\n");
    }

//-----------------------------------------------------------------------------
//  main  --  driver to use  plot_mass_radii() as a tool
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    char  *comment;
    bool  c_flag = FALSE;      /* if TRUE: a comment given on command line   */
    bool  t_flag = FALSE;      /* if TRUE: percentiles rather than quartiles */

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "c:t";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 't': t_flag = TRUE;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            

    dyn *b;
    while (b = get_dyn(cin))
	{
	if (t_flag)
	    plot_mass_radii_in_percentages(b);
	else
            plot_mass_radii(b);

	delete b;
	}
    }

#endif

// endof: lagradplot.C
