
local void dump_node(hdyn* b)
{
    // Incomplete dump of hdyn data.

    cerr << endl;
    int p = cerr.precision(HIGH_PRECISION);

    PRL(b->format_label());
    PRL(b->get_time());
    PRL(b->get_timestep());
    PRL(b->get_unperturbed_timestep());
    PRL(b->get_pos());
    PRL(b->get_vel());
    PRL(b->get_acc());
    PRL(b->get_jerk());
    if (b->get_kepler())
	b->get_kepler()->print_all(cerr);

    cerr.precision(p);
}

local void dump_binaries(hdyn* b)
{
    for_all_daughters(hdyn, b, bb)
	if (bb->get_oldest_daughter())
	    for_all_nodes(hdyn, bb, bbb)
		dump_node(bbb);
}


// Sample debugging lines to precede integrate_list() in evolve_system.


#if 0
	for (int i = 0; i < n_next; i++) {
	    hdyn* n = next_nodes[i];
	    if (n->is_parent() || n->is_low_level_leaf()) {
		PRL(t);
		n->print_label(cerr);
		dump_node(n);
	    }
	}
#endif

#if 0
	for (int i = 0; i < n_next; i++) {
	    if (next_nodes[i]) {
		hdyn* bi = next_nodes[i];
		if (bi->get_top_level_node()->is_parent()) {
		    cerr << endl;
		    cerr << "evolve_system: integration list contains "
			 << bi->format_label() << endl;
		}
	    }
	}
#endif

#if 0
	cerr << endl << "time = " << b->get_system_time() << endl;

	for (int j = 0; j < n_next; j++) {
	    hdyn* n = next_nodes[j];
	    cerr << n->format_label() << "  ";
	}
	cerr << endl;
#endif

//------------------------------------------------------------------------

// Sample debugging lines to follow integrate_list() in evolve_system.

#if 0
	for (int j = 0; j < n_next; j++) {
	    hdyn* n = next_nodes[j];
	    if (n && n->is_valid()) {
		if (node_contains(n->get_top_level_node(), "13a")) {
		    if (n->is_top_level_node())
			count0++;
		    hdyn* p = n->get_parent();
		    if (p->is_top_level_node())
			count1++;
		    hdyn* pp = p->get_parent();
		    if (pp && pp->is_top_level_node())
			count2++;
		    if ((count0+count1+count2)%1000 == 0) {
			cerr << "counts:  " << count0 << "  " << count1
			     << "  " << count2 << "  " << count_ << "  ";
			PRC(steps); PRL(cpu_time());
		    }
		} else
		    count_++;
	    }
	}
#endif


#if 0
 	if (t > 18.999) {
 	    cerr << endl; PRC(t); PRL(n_next);
	    for (int j = 0; j < n_next; j++) {
		hdyn* n = next_nodes[j];
		if (n && n->is_valid()) {
		    cerr << n->format_label() << ":  ";
		    if (n->is_low_level_node())
			n->print_perturber_list(cerr);
		    pp3(n, cerr);
		}
	    }
	    print_recalculated_energies(b, eps * eps);
 	}
	if (t > 19.1) exit(0);
#endif


#if 0
	if (t > 0.46) {

	    cerr << endl, PRC(t);
	    cerr << "evolve_system: after integrate_list:";

	    for (int i = 0; i < n_next; i++) {
		hdyn* bi = next_nodes[i];
		cerr << endl << "    ";
		cerr << bi->format_label();
		if (bi->is_top_level_node())
		    cerr << " (t";
		else
		    cerr << " (l";
		if (bi->is_leaf())
		    cerr << "l";
		else
		    cerr << "n";
		if (bi->get_kepler())
		    cerr << "u";
		else
		    cerr << "p";
		cerr << "), nn = ";
		if (bi->get_nn())
		    cerr << bi->get_nn()->format_label();
		else
		    cerr << "(nil)";
		if (bi->get_top_level_node()->is_parent()) {
		    cerr << endl;
		    pp3(bi->get_top_level_node(), cerr);
		}
	    }
	    cerr << endl;
	    print_recalculated_energies(b, eps * eps);

	} else {

	    for (int i = 0; i < n_next; i++) {
		if (next_nodes[i]) {
		    hdyn* bi = next_nodes[i];
		    if (bi->get_top_level_node()->is_parent()) {
			cerr << endl;
			cerr << "evolve_system: pp3 triggered by "
			     << bi->format_label() << endl;
			pp3(bi->get_top_level_node(), cerr);
		    }
		}
	    }

	}
#endif


#if 0
	for (int i = 0; i < n_next; i++) {
	    if (next_nodes[i] != NULL) {
		char c[128];
		real dt = next_nodes[i]->get_timestep();
		int idt = (int)(log(1.0/dt)/log(2.0) + .1);
		sprintf(c, "time = %.14f, new timestep = %.14f = 2^-%d",
			t, dt, idt);
		next_nodes[i]->log_comment(c);
	    }
	}
#endif

#if 0
	if (t > 2.61) {
	    for_all_nodes(hdyn, b, bb) {
		if (bb != b && bb->get_oldest_daughter()
		    && !bb->get_kepler()) {
		    if (!bb->get_valid_perturbers()) {
			PRC(t); cerr << "  no valid perturbers for "
			    	     << bb->format_label() << endl;
		    } else if (bb->get_n_perturbers() > 0) {
			for (int ip = 0; ip < bb->get_n_perturbers(); ip++) {
			    hdyn* p = bb->get_perturber_list()[ip];
			    if (!p) {
				PRC(t); cerr << "  NULL perturber #" << ip
				    << " for " << bb->format_label() << endl;
			    } else if (!p->is_valid()) {
				PRC(t); cerr << "  invalid perturber #" << ip
				    << " for " << bb->format_label() << endl;
			    }
			}			
		    }
		}
	    }
	}

	if (t > 28.1679) {
	    cerr << endl;
	    PRL(t);
	    for (int j = 0; j < n_next; j++) {
		hdyn* n = next_nodes[j];
		PRC(j);
		if (n) pp3(n);
		else cerr << endl;
	    }
	}

	if (t > 2.67) {
	    bool print_energy = false;
	    for (int j = 0; j < n_next; j++) {
		hdyn* n = next_nodes[j];
		if (n && n->is_valid()
		    && n->is_top_level_node())
		    print_energy = true;
	    }
	    if (print_energy)
		print_recalculated_energies(b, eps * eps);
	}


	if (t > 1.3 && t < 1.375) {
	    for (int j = 0; j < n_next; j++) {
		hdyn* n = next_nodes[j];
		if (n && n->is_valid()
		    && n->name_is("(1376b,1376a)")) {
		    cerr << "node " << n->format_label() << "  ";
		    PRL(t);
		    dump_story(n->get_dyn_story());
		}
	    }
	}



	if (t > 566.246119)
	    pp3((hdyn*)node_with_name("(1840,2812)", b->get_root()));

	if (fmod(count, 100) == 0) {

	    // Primitive checksum...

	    real checkx = 0, checkv = 0, checka = 0, checkj = 0;
	    for_all_nodes(hdyn, b, bb) {
		checkx += abs(bb->get_pos());
		checkv += abs(bb->get_vel());
		checka += abs(bb->get_acc());
		checkj += abs(bb->get_jerk());
	    }
	    int p = cerr.precision(HIGH_PRECISION);
	    cerr << endl; PRL(t);
	    PRC(checkx); PRL(checkv);
	    PRC(checka); PRL(checkj);
	    cerr.precision(p);
	}

	if (t > 566.246127) exit(0);
		

//	if (b->get_system_time() > 566.8)
	    for (int j = 0; j < n_next; j++) {
		hdyn* n = next_nodes[j];
		if (n && n->is_valid()
		    && n->is_low_level_node()
		    && fmod(n->get_steps(), 100.0) == 0
		    && n->name_is("1139")
		    && n->get_perturbation_squared() < 0.0001) {

		    pp3(n);

		    real dist = abs(n->get_pos()
				    - n->get_younger_sister()->get_pos());
		    real dtff = sqrt(dist*dist*dist
				     / (2*n->get_parent()->get_mass()));
		    real dtv  = sqrt(square(n->get_pos())
				     / square(n->get_vel()));

		    real kepstep = eta*min(dtff, dtv);
		    real ratio = n->get_timestep()/kepstep;
		    PRC(kepstep); PRL(ratio);
		}
	    }
#endif
