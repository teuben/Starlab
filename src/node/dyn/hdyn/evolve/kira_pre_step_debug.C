
// Place temporary pre-step debugging code here, to clean up kira.C.


#if 0
    if (t_prev == 252) {
	int ppre = cerr.precision(20);
	set_complete_system_dump(true);
	pp3(b);
	set_complete_system_dump(false);
	cerr.precision(ppre);
    }
#endif


#if 0
    if (t > 82.16960 && t < 82.1696256 && fmod(steps, 1) == 0) {

	cerr << "pre..." << endl;
	hdyn *top = NULL;
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid()) {
		cerr << n->format_label() << " ";
		if (n->name_is("445")) top = n->get_top_level_node();
	    }
	}
	cerr << endl;
	if (top) pp3(top);

	int p = cerr.precision(HIGH_PRECISION);
	PRL(t); 
	cerr.precision(p);
	print_recalculated_energies(b);

    }
#endif


#if 0
	if (t > 85.30960501 && t < 85.30960504) {
	    for (int ii = 0; ii < n_next; ii++) {
		hdyn *bb = next_nodes[ii];
		if (bb && node_contains(bb->get_top_level_node(), "829")) {
		    cerr << endl << "PRE..." << endl;
		    pp3(bb->get_top_level_node());
		    break;
		}
	    }
	}
#endif


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
	        cerr << "nodes (pre): ";
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
#endif


#if 0
 	for (int ii = 0; ii < n_next; ii++) {
 	    if (next_nodes[ii] && next_nodes[ii]->is_valid())
 		cerr << next_nodes[ii]->format_label() << " ";
 	}
 	cerr << endl; 	    
 	print_recalculated_energies(b);
 	pp3(next_nodes[0]->get_top_level_node());
#endif


#if 0	
 	real epot, ekin, etot;
 	int nunpert = 0;

  	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *bb = next_nodes[ii];
  	    if (bb && bb->is_valid()
 		&& bb->is_low_level_node()
 		&& bb->get_kepler()) nunpert++;
  	}

 	calculate_energies_with_external(b, epot, ekin, etot);
#endif
