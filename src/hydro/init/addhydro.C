/*
 *  addhydro.C: creates a hydro part for each body
 *.............................................................................
 *    version 1:  Jan 1993   Piet Hut
 *.............................................................................
 */
#include "hydro.h"

#ifndef TOOLBOX

/*-----------------------------------------------------------------------------
 *  addhydro  -- for all particles
 *-----------------------------------------------------------------------------
 */
void  addhydro(node * b, real R_eff)
    {
    node * bi;

    if (b->get_oldest_daughter() != NULL)
	for (bi=b->get_oldest_daughter(); bi != NULL;
	                                  bi=bi->get_younger_sister())
	    addhydro(bi, R_eff);
    else
	{
	hydrobase * old_hydrobase = b->get_hydrobase();
	hydro * new_hydro = new hydro(R_eff);

	new_hydro->set_hydro_story(old_hydrobase->get_hydro_story());
	b->set_hydrobase((hydrobase *) new_hydro);
	old_hydrobase->set_hydro_story(NULL);     // if not, deleting the old 
	delete old_hydrobase;                     // hydrobase would delete its
	}                                         // story as well
    }

#else

/*-----------------------------------------------------------------------------
 *  main  --n
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    int  c;
    bool  R_flag = FALSE;
    bool  c_flag = FALSE;
    real  R_eff = 0;           // default value;
    char  *comment;
    extern char *poptarg;

    while ((c = pgetopt(argc, argv, "R:c:")) != -1)
	switch(c)
	    {
	    case 'R': R_flag = TRUE;
		      R_eff = atof(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': cerr <<
		      "usage: addhydro [-R #] [-c \"..\"]\n";
		      exit(1);
	    }            
    
    node *b;

    while (b = get_node()) {
        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	addhydro(b, R_eff);

	put_node(b);
	delete b;
    }
}

#endif

/* endof: addhydro.c */
