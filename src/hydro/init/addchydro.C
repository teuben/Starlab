/*
 *  addchydro.C: creates a chydro part for each body
 *.............................................................................
 *    version 1:  Jan 1993   Piet Hut
 *.............................................................................
 */
#include "starlab_vector.h"
#include "node.h"
#include "chydro.h"

#ifndef TOOLBOX

/*-----------------------------------------------------------------------------
 *  addchydro  --
 *-----------------------------------------------------------------------------
 */
void  addchydro(node * b, real R_eff, real r_core, real m_core)
    {
    node * bi;

    if (b->get_oldest_daughter() != NULL)
	for (bi=b->get_oldest_daughter(); bi != NULL;
	                                  bi=bi->get_younger_sister())
	    addchydro(bi, R_eff, r_core, m_core);
    else
	{
	real M_tot = b->get_mass();

	if (m_core > M_tot)
	    {
	    cerr << "addchydro: m_core = " << m_core << " > M_tot = "
		 << M_tot << endl;
	    exit(1);
	    }

	if (r_core > R_eff)
	    {
	    cerr << "addchydro: r_core = " << r_core << " > R_eff = "
		 << R_eff << endl;
	    exit(1);
	    }

	hydrobase * old_hydrobase = b->get_hydrobase();
	chydro * new_chydro = new chydro(R_eff, m_core, r_core);

	new_chydro->set_hydro_story(old_hydrobase->get_hydro_story());
	b->set_hydrobase((hydrobase *) new_chydro);
	old_hydrobase->set_hydro_story(NULL);     // if not, deleting the old 
	delete old_hydrobase;                     // hydrobase would delete its
	}                                          // story as well
    }

#else

/*-----------------------------------------------------------------------------
 *  main  --
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    int  c;
    bool  r_flag = FALSE;
    bool  R_flag = FALSE;
    bool  m_flag = FALSE;
    bool  c_flag = FALSE;
    real  R_eff = 0;           // default value;
    real  r_core = 0;           // default value;
    real  m_core = 0;           // default value;
    char  *comment;
    extern char *poptarg;

    while ((c = pgetopt(argc, argv, "R:r:m:c:",
		    "$Revision$", _SRC_)) != -1)
	switch(c)
	    {
	    case 'R': R_flag = TRUE;
		      R_eff = atof(poptarg);
		      break;
	    case 'r': r_flag = TRUE;
		      r_core = atof(poptarg);
		      break;
	    case 'm': m_flag = TRUE;
		      m_core = atof(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': cerr <<
		      "usage: addchydro [-R #] [-r #] [-m #] "
 		      << "[-c \"..\"]\n" <<
		      "  for effective radius, core radius, core mass\n";
		      exit(1);
	    }
    
    node *b;

    while (b = get_node()) {
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

        addchydro(b, R_eff, r_core, m_core);
	put_node(b);
	delete b;
    }
}

#endif

/* endof: addchydro.c */
