
// Place temporary post-step debugging code here, to clean up kira.C.

#if 0
    if (t > 252 && t < 252.1) {
	cerr << endl;
	int ppost = cerr.precision(20);
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid() && n->name_is("3")) {
		PRC(n->format_label());
		PRL(n->get_next_time());
		pp3(n);
	    }
	}
	cerr.precision(ppost);
    }
#endif


#if 0
if (t > 255.29) {
    for (int ii = 0; ii < n_next; ii++) {
	if (next_nodes[ii] && next_nodes[ii]->is_valid())
	    cerr << next_nodes[ii]->format_label() << " ";
    }
    cerr << endl;
    int p = cerr.precision(20);
    PRL(t);
    hdyn *bbb = (hdyn*)node_with_name("(68,629)", b);
    if (bbb) cerr << "(68,629) t = " << bbb->get_time()
		  << " dt = " << bbb->get_timestep()<< endl;
    bbb = (hdyn*)node_with_name("674", b);
    if (bbb) cerr << "674 t = " << bbb->get_time()
		  << " dt = " << bbb->get_timestep()<< endl;
    for (int ii = 0; ii < n_next; ii++) {
	hdyn *n = next_nodes[ii];
	if (n && n->is_valid() && n->name_is("(68,629)")) {
	    PRC(n->format_label());
	    PRL(n->get_next_time());
	    pp3(n);
	}
    }
    cerr.precision(p);
}
#endif


#if 0
    if (t > 82.16960 && t < 82.1696256 && fmod(steps, 1) == 0) {

	cerr << "post..." << endl;
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
    if (t > 82.1421279 && t < 82.142131 && fmod(steps, 1) == 0) {

	cerr << "post..." << endl;
	hdyn *top = NULL;
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid()) {
		cerr << n->format_label() << " ";
		if (n->name_is("(925,6)")) top = n->get_top_level_node();
		if (n->name_is("925")) top = n->get_top_level_node();
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
    if (t > 82.1420 && t < 82.1423 && fmod(steps, 100) == 0)
	print_recalculated_energies(b);
#endif


#if 0
    fmod(steps, 100) == 0) print_recalculated_energies(b);
#endif


#if 0
    if (((t > 81.814 && t < 81.815) || (t > 81.829 && t < 81.832))
	 && fmod(steps, 1) == 0) {

	for (int ii = 0; ii < n_next; ii++) {
	    if (next_nodes[ii] && next_nodes[ii]->is_valid())
		cerr << next_nodes[ii]->format_label() << " ";
	}
	cerr << endl; 

	int p = cerr.precision(HIGH_PRECISION);
	PRC(t); 
	cerr.precision(p);
	print_recalculated_energies(b);

    }
#endif


#if 0
    if (   t > 81.500675 && t < 81.502364
	|| t > 81.508965 && t < 81.510330
	|| t > 81.535075 && t < 81.536905
	|| t > 81.558855 && t < 81.560755
	|| t > 81.566103 && t < 81.567278
	|| t > 81.574715 && t < 81.576000
	|| t > 81.588040 && t < 81.589356 )
	print_recalculated_energies(b);
#endif


#if 0
    fmod(steps, 100) == 0) print_recalculated_energies(b);
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
