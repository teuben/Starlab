
// Place temporary post-step debugging code here, to clean up kira.C.


    // if (fmod(steps, 100) == 0) print_recalculated_energies(b);


#if 0
	if (t >= 44.1875 && t <= 44.21875) {

	  int nstep = 100;
	  if (   (t >= 44.1903765  && t <= 44.1903767)
	      || (t >= 44.19358027 && t <= 44.1935804)
	      || (t >= 44.19804955 && t <= 44.19805098)) nstep = 1;

	    int p = cerr.precision(INT_PRECISION);
	    cerr << endl; PRL(t);
	    cerr.precision(p);

	    if (fmod(steps, nstep) == 0) {
	        cerr << "nodes (post): ";
		for (int ii = 0; ii < n_next; ii++) {
		    if (next_nodes[ii] && next_nodes[ii]->is_valid())
		      cerr << next_nodes[ii]->format_label() << " ";
		}
		cerr << endl; 	    
		print_recalculated_energies(b);
	    }

	    for (int ii = 0; ii < n_next; ii++) {
		if (next_nodes[ii] && next_nodes[ii]->is_valid()
		    && next_nodes[ii]->is_top_level_node()
		    && next_nodes[ii]->n_leaves() > 2)
		    pp3(next_nodes[ii]);
	    }

	}
	if (t > 44.1935804) exit(0);
#endif


#if 0
     for (int ii = 0; ii < n_next; ii++) {
 	if (!next_nodes[ii] || !next_nodes[ii]->is_valid()) {
 	    cerr << "next_node[" << ii << "] = " << next_nodes[ii]
 		 << " is now invalid" << endl << flush;
 	}
     }
#endif


#if 0
//  	if (b->get_system_time() >= 1.46875
//  	    && b->get_system_time() <= 1.5) {
	    int countii = 0;
  	    for (int ii = 0; ii < n_next; ii++) {
  		if (next_nodes[ii] && next_nodes[ii]->is_valid()
  		    && node_contains(next_nodes[ii]->get_top_level_node(),
				     "100393")) {
		    if (countii++ == 0) {
			int p = cerr.precision(HIGH_PRECISION);
			cerr << endl << "After step to time "
			     << b->get_system_time() << endl;
			cerr.precision(p);
		    }
		    PRL(next_nodes[ii]->format_label());
  		    pp3(next_nodes[ii]->get_top_level_node());
		}
  	    }
//  	}
#endif


#if 0
	if (t > 85.30960501 && t < 85.30960504) {
	    for (int ii = 0; ii < n_next; ii++) {
		hdyn *bb = next_nodes[ii];
		if (bb && bb->is_valid()
		    && node_contains(bb->get_top_level_node(), "829")) {
		    cerr << endl << "POST..." << endl;
		    pp3(bb->get_top_level_node());
		    break;
		}
	    }
	}
#endif


#if 0
	cerr << endl << "check_sync #2 at t = "
	     << b->get_system_time() << " (";
	xprint(b->get_system_time(), cerr, false);
	cerr << ")" << endl;
	int n_unp = 0;
	for (int ii = 0; ii < n_next; ii++) {
	    if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
		hdyn *bb = next_nodes[ii];
		if (bb->get_kepler()) n_unp++;
	    }
	}
	PRC(n_next); PRL(n_unp);
	for (int ii = 0; ii < min(2,n_next); ii++) {
	    if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
		hdyn *bb = next_nodes[ii];
		if (ii > 0) cerr << "  ";
		cerr << bb->format_label() << " "
		     << bb->get_timestep() << " ";
		xprint(bb->get_time(), cerr, false);
	    }
	}
	cerr << endl;
	    
	for (int ii = 0; ii < n_next; ii++) {
	    if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
		hdyn *bb = next_nodes[ii];

		if (bb->get_time() != b->get_system_time()
		    && bb->get_kepler() == NULL) {
		    cerr << "check_sync warning:  node "
			 << bb->format_label() << " not synchronized" << endl
			 << "                     system time = "
			 << b->get_system_time() << " (";
		    xprint(b->get_system_time(), cerr, false);
		    cerr << ")" << endl;
		    cerr << "                     node time   = "
			 << bb->get_time() << " (";
		    xprint(bb->get_time(), cerr, false);
		    cerr << ")" << endl;
		}
	    }
	}
	cerr << "...done" << endl;
#endif
