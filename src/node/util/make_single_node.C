
//// mk_single_node:  create a single node, to which other nodes can be
////                  added later, using `add_daughter_node'.
////
//// Options:    -c  comment            optional comment
////             -m  mass               optional mass [1]
////             -i                     label nodes
////

//   version 1:  Dec 1994   Piet Hut
//               1998-11-25 Peter Teuben        fixed help


#include "node.h"

/*===========================================================================*/

#ifdef TOOLBOX

/*-----------------------------------------------------------------------------
 *  main  --  driver to create directly a single node
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    bool  i_flag = FALSE;
    real m = 1;               // default mass: unity
    char  *comment;

    check_help();

    extern char *poptarg;
    int  c;
    char* param_string = "c:im:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'i': i_flag = TRUE;
		      break;
	    case 'm': m = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	}            
    
    node * root;

    root = new node();
    root->set_mass(m);

    if (c_flag == TRUE)
        root->log_comment(comment);
    root->log_history(argc, argv);

    if (i_flag)
        root->set_label(1);

    put_node(root);
    rmtree(root);
}

#endif

/*===========================================================================*/

/* endof: mk_single_node.C */

