
// Place temporary post-step debugging code here, to clean up kira.C.
// Current time is t, root node is b.  The integration list is
// next_nodes[], of length n_next.  The current time is t.  Real steps
// counts particle steps; real count counts block steps.  The variable
// tdbg may be used here to control debugging output.  It is
// initialized to -1.

#if 0
if (count > 163000) {
  b->get_kira_diag()->n_check_heartbeat = 1;
  int p = cerr.precision(20); PRC(t); cerr.precision(p); PRL(n_next);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && !n->get_kepler()) {
      PRC(n->format_label()); PRL(n->get_timestep());
    }
  }
 }
if (count > 165000) exit(0);
#endif

#if 0
if (count > 188000 && count < 205000) {
  b->get_kira_diag()->n_check_heartbeat = 1;
  int p = cerr.precision(20);
  PRC(t); PRL(n_next); cerr.precision(p);

  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && !n->get_kepler()) {
      PRC(n->format_label()); PRL(n->get_timestep());
#if 0
      if (n->name_is("(114,10114)")
	  || n->name_is("((1092,11092),222)")
	  || n->name_is("((114,10114),((1092,11092),222))")) {

	hdyn *b1 = (hdyn*)node_with_name("(114,10114)",
					 n->get_root());
	hdyn *b2 = (hdyn*)node_with_name("((1092,11092),222)",
					 n->get_root());
	real m1 = b1->get_mass();
	real m2 = b2->get_mass();

	cerr << "%%% ";
	vec force1 = m1*b1->get_acc();
	vec force2 = m2*b2->get_acc();
	int p = cerr.precision(HIGH_PRECISION); PRC(t); cerr.precision(p);
	PRC(count); PRC(force1); PRC(force2); PRL(force1+force2);

	real mu = m1*m2/(m1+m2);
	vec dx = b2->get_pos() - b1->get_pos();
	vec dv = b2->get_vel() - b1->get_vel();
	real dx1 = abs(dx);
	PRC(dx1); PRL(abs(dv));
	real E12 = 0.5*mu*square(dv)-m1*m2/dx1;
	PRC(E12/mu); PRL(E12);

	real dx3 = pow(dx1,3);
	real a1 = abs(b1->get_acc());
	real a2 = abs(b2->get_acc());
	PRC(a1); PRC(abs(m2*dx/dx3)); PRC(a2); PRL(abs(-m1*dx/dx3));
      }

      if (n->name_is("((114,10114),((1092,11092),222))"))
	pp3(n);
      else if (n->name_is("((1092,11092),222)")) {
	n->print_perturber_list();
	n->get_oldest_daughter()->print_pert();
      } else if (n->name_is("(1092,11092)"))
	n->print_pert();
      else if (n->name_is("114"))
	n->print_pert();
#endif
    }
  }
 }
if (count > 205000)
  b->get_kira_diag()->n_check_heartbeat = 1000;
#endif

#if 0
if (t > 0.03125) {
  int p = cerr.precision(20);
  PRC(t); PRC(get_effective_block(t)); PRL(n_next);
  cerr.precision(p);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && !n->get_kepler()) {
      PRC(n->format_label()); PRL(n->get_timestep());
      PRL(fmod(t,n->get_timestep()));
    }
  }
  if (t > 0.1) exit(1);
 }
#endif

#if 0
if (t >= 0) {
  real n2 = 1;
  int block = 1;
  while (fmod(t, n2) != 0) {block++; n2 /= 2;}
  if (block > 30 || t > 0.049) {
    PRC(block); PRL(n2);
    PRL(n_next);
    for (int ii = 0; ii < n_next; ii++) {
      hdyn *n = next_nodes[ii];
      if (n && n->is_valid()) {
	PRL(n->format_label());
	int p = cerr.precision(20);
	PRI(4); PRL(n->get_time());
	cerr.precision(p);
	PRI(4); PRL(n->get_timestep());
      }
    }
  }
 }
#endif

#if 0
if (t > 0.6) {
  hdyn *n = (hdyn*)node_with_name("233", b);
  if (n) {
    PRC(n->format_label());
    PRC(n->get_system_time());
    PRL(n->get_timestep());
  }
  n = (hdyn*)node_with_name("(233,10233)", b);
  if (n) {
    PRC(n->format_label());
    PRC(n->get_system_time());
    PRL(n->get_timestep());
  }
  real n2 = 1;
  int block = 1;
  while (fmod(t, n2) != 0) {block++; n2 /= 2;}
  PRC(block); PRL(n2);
}
#endif

#if 0
if (t > 0.045) {
  int p = cerr.precision(20);
  cerr << endl; cerr << "after "; PRC(t); PRL(t_prev);
  cerr.precision(p);
  //dump_node_list();
  hdyn *n = (hdyn*)node_with_name("(346,10346)", b);
  if (n) {
    PRC(n->format_label());
    PRL(n->get_timestep());
    int p = cerr.precision(20);
    PRC(n->get_time());
    PRL(n->get_next_time());
    cerr.precision(p);
  }
  n = (hdyn*)node_with_name("346", b);
  PRC(n->get_kepler()); PRL(n->get_fully_unperturbed());
  dump_node_list_for("(346,10346)");
  if (t > 0.051) {
    for (int ii = 0; ii < n_next; ii++) {
      hdyn *n = next_nodes[ii];
      if (n && n->is_valid()) {
	PRL(n->format_label());
	int p = cerr.precision(20);
	PRI(4); PRL(n->get_time());
	cerr.precision(p);
	PRI(4); PRL(n->get_timestep());
      }
    }
  }
}
#endif

#if 0
if (t < 400.5) {
  int p = cerr.precision(15);
  cerr << endl; cerr << "after "; PRL(t);
  cerr.precision(p);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && n->is_low_level_node()) {
      PRC(n->format_label()); PRL(n->get_timestep());
      PRL(n->get_valid_perturbers());
      PRL(n->get_parent()->get_valid_perturbers());
      PRL(n->get_perturbation_squared());
      PRL(n->get_kappa());
      hdyn *pnn = n->get_parent()->get_nn();
      PRL(pnn);
      if (pnn) PRL(pnn->format_label());
    }
  }
} else
  exit(1);
#endif

#if 0
if (t > 83.29) {
  int p = cerr.precision(15);
  cerr << endl; cerr << "after "; PRL(t);
  cerr.precision(p);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && (n->is_parent() || n->is_low_level_node())) {
      PRC(n->format_label()); PRL(n->get_timestep());
      if (n->get_kepler()) PRL(n->get_unperturbed_timestep());
    }
  }
}
#endif

#if 0
if (t > 16.4342) {
  bool first = true;
  int p = cerr.precision(15);
  cerr << endl; PRC(t); PRL(n_next);
  cerr.precision(p);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && (n->is_parent() || n->is_low_level_node())) {
      PRC(ii); PRC(n->format_label()); PRL(n->get_timestep());
      if (n->is_low_level_node() && n->get_timestep() < 1.e-9) {
	int prec = cerr.precision(HIGH_PRECISION);
	PRL(n->get_perturbation_squared());
	cerr.precision(prec);
	PRL(n->get_parent()->get_valid_perturbers());
	PRL(n->get_kepler());
	hdyn *p = n->get_top_level_node();
	hdyn *pnn = NULL;
	real r2nn = VERY_LARGE_NUMBER;
	for_all_daughters(hdyn, b, bb) {
	  if (bb != p) {
	    real r2 = square(bb->get_pos()-p->get_pos());
	    if (r2 < r2nn) {
	      pnn = bb;
	      r2nn = r2;
	    }
	  }
	}
	PRC(pnn->format_label()); PRC(pnn->get_mass()); PRL(sqrt(r2nn));
      }
    }
  }
}
#endif

#if 0
if (t > 1360.926581) {
  int p = cerr.precision(15);
  cerr << endl; PRL(t);
  cerr.precision(p);
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && (n->is_parent() || n->is_low_level_node())) {
      PRC(n->format_label()); PRL(n->get_timestep());
      if (n->is_low_level_node())
	PRC(n->get_posvel()); PRL(n->get_pos()*n->get_vel());
    }
  }
}
#endif

#if 0
if (t > 361.1566) {
  for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid() && n->get_oldest_daughter()) {
      hdyn *od = n->get_oldest_daughter();
      hdyn *yd = od->get_younger_sister();

      // Compute the CM acceleration the hard way, using the CM
      // approximation for other stars.

      // Compute the forces on od and yd, the hard way...

      vec xod = n->get_pos() + od->get_nopred_pos();
      vec xyd = n->get_pos() + yd->get_nopred_pos();
      vec acm = 0, aod = 0, ayd = 0;
      for_all_daughters(hdyn, od->get_root(), bb) {
	if (bb != n) {
	  vec x = bb->get_nopred_pos() - n->get_pos();
	  real r2 = x*x + n->get_eps2();
	  acm += bb->get_mass()*x/(r2*sqrt(r2));
	  x = bb->get_nopred_pos() - xod;
	  r2 = x*x + od->get_eps2();
	  aod += bb->get_mass()*x/(r2*sqrt(r2));
	  x = bb->get_nopred_pos() - xyd;
	  r2 = x*x + yd->get_eps2();
	  ayd += bb->get_mass()*x/(r2*sqrt(r2));
	}
      }
      vec x12 = yd->get_nopred_pos() - od->get_nopred_pos();
      real r2 = x12*x12 + n->get_eps2();
      vec a12 =  yd->get_mass()*x12/(r2*sqrt(r2));
      vec a21 = -od->get_mass()*x12/(r2*sqrt(r2));
      vec atop = (od->get_mass()*aod + yd->get_mass()*ayd)/n->get_mass();
      cerr << "comparing CM acc on " << n->format_label() << endl;
      PRL(acm);
      PRL(aod);
      PRL(ayd);
      PRL(a12);
      PRL(atop);
      PRL(n->get_acc());
      PRL(aod+a12);
      PRL(aod+a12-atop);
      PRL(od->get_acc());
      PRL(ayd+a21);
      PRL(ayd+a21-atop);
      PRL(yd->get_acc());
      vec od_abs_acc = hdyn_something_relative_to_root(od, &hdyn::get_acc);
      PRL(od_abs_acc);
    }
  }
}
#endif

#if 0
for (int ii = 0; ii < n_next; ii++) {
    hdyn *n = next_nodes[ii];
    if (n && n->is_valid()) {
	if (n->name_is("174") || n->name_is("881")) {
	  PRC(n->format_label()); PRL(n->get_timestep());
	}
    }
}
#endif

#if 0
//    if (t >= 11.868 && t <= 11.8695) {
    if (t >= 11.868 && t <= 11.870) {
	int p = cerr.precision(10);
	cerr << "post: "; PRC(t);
	cerr.precision(p);
	PRL(n_next);
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid()) {
		PRI(4); PRC(n->format_label());	PRL(n->get_timestep());
		if (n->is_parent()) pp3(n);
	    }
	}
	if (t > tdbg) {
	    print_recalculated_energies(b);
	    while (t > tdbg) tdbg += 1.e-6;
	}
    }
#endif

#if 0
    if (t > 1.10503 && t < 1.10608) {

	// Output part of the system for xstarplot purposes.

	int p = cout.precision(10);
	cout << ";; system_time = " << t << endl;
	int n = 0;
	int nmax = 10;
	vec pos, vel, pos0 = vec(0.1, 0.2, -0.1);
 	for_all_leaves(hdyn, b, bb) {
	    pos = hdyn_something_relative_to_root(bb,
						  &hdyn::get_nopred_pos);
	    vel = hdyn_something_relative_to_root(bb,
						  &hdyn::get_nopred_vel);
	    if (abs(pos-pos0) < 0.01 && ++n < 25)
		cout << bb->get_index() << " "
		     << bb->get_mass() << " "
		     << pos-pos0 << " " << vel << endl;
	}
	if (n < nmax) {
	    for (int i = 0; i < nmax-n; i++) {
		pos -= vec(10);
		vel -= vec(10);
		cout << "999 0.0001 " << pos-pos0 << " " << vel << endl;
	    }
	}
	cout.precision(p);
	cout << endl;
    }
#endif


#if 0
//    if (t > 1.10503 && t < 1.10608) {
    if (t > 1.1332 && t < 1.14) {
	int p = cerr.precision(10);
	PRC(t);
	cerr.precision(p);
	PRL(n_next);
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid()) {
		PRI(4); PRC(n->format_label());	PRL(n->get_timestep());
#if 1
		if (n->is_parent() && n->is_top_level_node()
		    && node_contains(n, 1267)) {
		    pp3(n);
#if 1
		    vec act = n->get_acc();

		    hdyn* bb = (hdyn*)node_with_name("1267", b);
		    hdyn *bskip1 = bb->get_top_level_node();
		    vec pos = hdyn_something_relative_to_root(bb,
							&hdyn::get_nopred_pos);
		    vec dpos = pos - n->get_pos();
		    real dr1 = abs(dpos);
		    vec acc = bb->get_mass()*dpos/pow(dr1,3);

		    bb = (hdyn*)node_with_name("267", b);
		    hdyn *bskip2 = bb->get_top_level_node();
		    pos = hdyn_something_relative_to_root(bb,
							&hdyn::get_nopred_pos);
		    dpos = pos - n->get_pos();
		    real dr2 = abs(dpos);
		    acc += bb->get_mass()*dpos/pow(dr2,3);

		    bb = bb->get_parent();
		    pos = hdyn_something_relative_to_root(bb,
							&hdyn::get_nopred_pos);
		    dpos = pos - n->get_pos();
		    real dr = abs(dpos);
		    vec acm = bb->get_mass()*dpos/pow(dr,3);

		    PRL(act);
		    PRC(acc); PRC(dr1); PRL(dr2);
		    PRC(acm); PRL(dr);

		    vec alt = 0;
		    for_all_daughters(hdyn, b, bb) {
			if (bb != n && bb != bskip1 && bb != bskip2) {
			    pos = hdyn_something_relative_to_root(bb,
							&hdyn::get_nopred_pos);
			    dpos = pos - n->get_pos();
			    real dr = abs(dpos);
			    alt += bb->get_mass()*dpos/pow(dr,3);
			}
		    }

		    PRL(alt);
#endif
		}
#endif
	    }
	}
	// print_recalculated_energies(b);
    }
//if (t > 1.14) exit(0);
#endif


#if 0
    int ppp = cerr.precision(HIGH_PRECISION);
    PRL(ds);
    pp3("1001");
    pp3("1002");
    cerr.precision(ppp);
    cerr << endl;
#endif


#if 0
    if ((t > 1250.854692 && t < 1250.854695)
	|| (t > 1250.854730 && t < 1250.854732)
	|| (t > 1250.854800 && t < 1250.854801)
	|| (t > 1250.8548319)) {
	for (int ii = 0; ii < n_next; ii++)
	if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
	    hdyn *n = next_nodes[ii];
	    hdyn *top = n->get_top_level_node();
	    if (node_contains(top, "97")) pp3(top);
	}
    }
#endif


#if 0
    if (t > 122.36 && t < 122.37) {
	for (int ii = 0; ii < n_next; ii++)
	if (next_nodes[ii] && next_nodes[ii]->is_valid()) {
	    hdyn *n = next_nodes[ii];
	    hdyn *top = n->get_top_level_node();
	    if (node_contains(top, "(94,10)")
		|| node_contains(top, "27"))
		pp3(n);
	}
	hdyn *n = (hdyn*)node_with_name("(94,10)", b);
	int p = cerr.precision(HIGH_PRECISION);
	PRC(n->get_time()); PRL(n->get_timestep());
	cerr.precision(p);
    }
#endif


#if 0
    if (t_prev == 252) {
	int ppost = cerr.precision(20);
	cerr << endl;
	for (int ii = 0; ii < n_next; ii++) {
	    hdyn *n = next_nodes[ii];
	    if (n && n->is_valid()) {
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
