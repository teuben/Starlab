/*
 *  tripev.C: evolves triple system
 *.............................................................................
 *    version 1:  Okt 1993   Simon F. Portegies Zwart 
 *.............................................................................
 */
#include "stdinc.h"
#include "double_star.h"
#include "single_star.h"
#include "main_sequence.h"
#include <stdlib.h>
#include <fstream.h>

#define FULL_OUTPUT

#  define  FALSE  0
#  define  TRUE   1

/*-----------------------------------------------------------------------------
 *  main  --
 *	usage:
 *		tripev -t # [options]  ,
 *
 *		where # is the initial age of the cluster.
 *	options:
 *            	The following options are allowed:
 *	cluster age:
 *		-t #	Where # stands for the initial age of the
 *			in Myear.
 *
 *		At present the running time of the integrator correspnds
 *		to the stellar age an a one by 10e6year basis.
 *		This however should be scaled to the cluster parameters.
 *-----------------------------------------------------------------------------
 */

/*-----------------------------------------------------------------------------
 *  binev  --
 *-----------------------------------------------------------------------------
 */
local void  binev(double_star * b, real start_time, real end_time, 
                  int n_steps) {

//		Setup star from input data.
    real ageint = (end_time - start_time)/n_steps;

    real time_done, until_time, dt;
    for(int j = 0; j<n_steps; j++) {
           b->evolve_element(start_time+ageint*(j+1));
           cout << b->get_binary_age() <<":";
           b->put_state();
#ifdef FULL_OUTPUT
    b->print_status();
#endif
    }

/*
cout << "result: " << type_string(b->get_bin_type()) 
                  << " (" << type_string(b->get_primary()->get_element_type())
                  << ", " << type_string(b->get_secondary()->get_element_type())
                  << ")         Random seed: " << seed << endl; 
cout <<"M = "<<b->get_primary()->get_total_mass()<<"\t"
     <<"R = "<<b->get_primary()->get_radius()<<endl;
cout <<"m = "<<b->get_secondary()->get_total_mass()<<"\t"
     <<"r = "<<b->get_secondary()->get_radius()<<endl;
cout <<"a = "<<b->get_semi()<<endl;
cout <<"e = "<<b->get_eccentricity()<<endl;
*/
}     

/* endof: addstar.C */

int main(int argc, char ** argv)
    {
/*
    main_sequence * m = new main_sequence();
    single_star * p = new single_star(m);
    main_sequence * n = new main_sequence();
    single_star * s = new single_star(n);
    double_star * b = new double_star(p, s);
*/

    star_state primary, secondary;
//    double_star * b = new double_star();
    double_init inner;
    double_init outer;
    int  c;
    bool  M_flag = FALSE;
    bool  Q_flag = FALSE;
    bool  c_flag = FALSE;
    bool  f_flag = FALSE;
    inner.start_time = outer.start_time = 0;           // default value;
    inner.end_time = outer.end_time = 23;
    inner.n_steps = outer.n_steps = 1;
    inner.mass_prim = 12;
    inner.q = 0.75;
    outer.mass_prim = (1+inner.q) * inner.mass_prim;
    outer.q = 0.2;
    inner.semi = 250;
    outer.semi = 1000;
    inner.eccentricity = 0.2;
    outer.eccentricity = 0.1;
    int   id = 1;
    char  *comment;
    int   seed = 0;
    int   srandinter(int);
    extern char *poptarg;
    int   pgetopt(int, char **, char *);
//    local void  binev(double_star *, real, real, int);
    
    if (argc <= 1)
       cerr <<"Default: binev -t 0 -T 23 -a 250 -e 0.25 -M 12 -m 10 -n 1 \n";
 
    while ((c = pgetopt(argc, argv, "t:n:M:Q:q:A:a:E:e:T:fs:c:")) != -1)
	switch(c)
	    {
            case 't': inner.start_time = 
                      outer.start_time = atof(poptarg);
                      break;
            case 'n': inner.n_steps = 
                      outer.n_steps = atoi(poptarg);
                      break;
            case 'M': M_flag = TRUE;
		      inner.mass_prim = 
                      outer.mass_prim = atof(poptarg);
                      break;
            case 'Q': Q_flag = TRUE;
  		      inner.q = atof(poptarg);
                      break;
            case 'q': outer.q = atof(poptarg);
                      break;
            case 'A': inner.semi   = atof(poptarg);
                      break;
            case 'a': outer.semi   = atof(poptarg);
                      break;
            case 'E': inner.eccentricity = atof(poptarg);
                      break;
            case 'e': outer.eccentricity = atof(poptarg);
                      break;
            case 'T': inner.end_time = 
                      outer.end_time = atof(poptarg);
                      break;
            case 'f': f_flag = TRUE;
                      break;
            case 's': seed = atoi(poptarg);
                      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': cerr <<"usage: binev [-t #] [-T #] "
                           <<             "[-M #] [-n #] "
                           <<             "[-Q #] [-q #] "
                           <<             "[-A #] [-a #] "
                           <<             "[-E #] [-e #] "
                           <<             "[-s #] [-c \"..\"] "
                           <<             "[-f < filename] \n";
		      exit(1);
	    }

    if (inner.end_time<=inner.start_time) {
       cerr << "time of final stage should exceed initial stage time" 
            << " (" << inner.start_time << ">" << inner.end_time << ")" << endl;
       exit(1);
    }

    if (M_flag || Q_flag)
       outer.mass_prim = (1+inner.q) * inner.mass_prim;

    double_star * t = triple_star(inner, outer);
/*
//    main_sequence * m = new main_sequence();
    single_star * p = new single_star();
    p->initialize(m_prim, id, t_start);
//    main_sequence * n = new main_sequence();
    single_star * s = new single_star();
    s->initialize(m_sec, id, t_start);
    double_star * b = new double_star(p, s, semi, ecc, t_start);
   
//    main_sequence * o = new main_sequence();
    single_star * d = new single_star();
    d->initialize(m_tert, id, t_start);
    double_star * t = new double_star(b, d, a2, e2, t_start);

    int random_seed = srandinter(seed);
    cout << "Random seed: " << random_seed << endl;
*/
/*
    int id = 1;
    if(!f_flag)
       b->initialize(m_prim, m_sec, semi, ecc, id, t_start);
    else
       b->initialize(id, t_start);
*/

//    primary.init_star_appeal((star*)b->get_primary());
//    secondary.init_star_appeal((star*)b->get_primary());
#ifdef FULL_OUTPUT
    t->put_element();
    t->print_status();
#endif

    binev(t, inner.start_time, inner.end_time, inner.n_steps);

/*
    if (primary.identity == t->get_primary()->get_identity()) {
       primary.make_star_appeal((star*)t->get_primary());
       secondary.make_star_appeal((star*)t->get_secondary());
       primary.put_star_appeal();
       cout <<", ";
       secondary.put_star_appeal();
    }
    else {
       secondary.make_star_appeal((star*)b->get_primary());
       primary.make_star_appeal((star*)b->get_secondary());
       secondary.put_star_appeal();
       cout <<", ";
       primary.put_star_appeal();
    }
*/

 }
