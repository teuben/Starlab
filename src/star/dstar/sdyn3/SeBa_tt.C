#ifndef TOOLBOX

#include "SeBa_tt.h"

void merge_binary_components(dyn* bi) {

  //      print_result(target);
  //  cerr<<"merge_binary_components()"<<endl;
  dyn * root = bi->get_parent();
  dyn *pi = bi->get_oldest_daughter();
  pi->set_label(bi->get_index());
  pi->set_root(NULL);

  add_node(pi, root);
  detach_node_from_general_tree(bi);
}

void ionize_binary(dyn* target) {

  //  cerr<<"ionize_binary()"<<endl;
	
  dyn * cluster = target->get_parent();

  double_star * ds = dynamic_cast(double_star*, target->get_starbase());
  binary_type bin_type = ds->get_bin_type();

  dyn * root = target->get_parent();
  dyn *od = target->get_oldest_daughter();
  dyn *ys = target->get_oldest_daughter()->get_younger_sister();
  detach_node_from_general_tree(target);
  add_node(od, root);
  add_node(ys, root);

}

bool repair_tree(dyn * root) {

  //  cerr << "Repair tree"<<endl;
  bool repair = false;

  for_all_nodes(dyn, root, bi) {
    if (bi->is_parent() && !bi->is_root()) {

      double_star * ds = dynamic_cast(double_star*, bi->get_starbase());
      binary_type bin_type = ds->get_bin_type();
      
      //      cerr << "Binary type: " << type_string(bin_type) << endl;

      switch(ds->get_bin_type()) {
      case Merged:  
	merge_binary_components(bi);
	repair = true;
	break;
      case Disrupted:
	ionize_binary(bi);
	//	cerr << "Disrupt binary" << endl;
	repair = true;
	break;
	//      default:
	//	cerr << "No default in switch ds->get_element_type()"<<endl;
	// do nothing
      };
    }
    else {
      // single star does not require repair. 
      // Could and should eliminate the disintegrated stars

      single_star * ss = dynamic_cast(single_star*, bi->get_starbase());
      switch(ss->get_element_type()) {
      case Disintegrated:
	//	eliminate_single_star(bi);
	//	cerr << "Disintegrate single star" << endl;
	repair = true;
	break;
	//default:
	//	cerr << "No default in switch ss->get_element_type()"<<endl;
	// do nothing
      };
    }
  }
      
  return repair;
}

void update_binary_parameters(dyn *bi, 
			      final_state3 &final3) {

  double_star * ds = dynamic_cast(double_star*, bi->get_starbase());
  //  ds->dump(cerr);

  // update binary parameters
  ds->set_semi(((double_star*)ds)->get_semi()*final3.sma);
  ds->set_eccentricity(final3.ecc);
  //  ds->dump(cerr);
}

// Note that the velocity of the particularparticle and the binary are
// used. All stars (including the binary) are in equilibrium
// (equipartition of kinetic energy).
void preserve_initial_state(dyn* bi, dyn *si,
			    final_state3& final3) {
  //  cerr<<"preserve_initial_state"<<endl;

  double_star * ds = dynamic_cast(double_star*, bi->get_starbase());
  //  binary_type bin_type = ds->get_bin_type();
  
  //  ds->dump(cerr);
  ds->set_semi(ds->get_semi()*final3.sma);
  ds->set_eccentricity(final3.ecc);
  //  ds->dump(cerr);
 
  // Maybe we should still do something with the velocity and 
  // binary binding energy. 
}

// Exchange incoming 'bullet' star in binary, leaving the secondary intact; 
// replacing the primary star
void exchange_secondary(dyn *bi,
			dyn *si,
			final_state3& final3) {

  //      print_result(target);

  //  cerr<<"exchange_secondary"<<endl;

  dyn *root = bi->get_parent();
  dyn *pi = bi->get_oldest_daughter();
  starbase *sbpi = (starbase*)pi->get_starbase();
  starbase *sbsi = (starbase*)si->get_starbase();
  pi->set_starbase(sbsi);
  sbsi->set_node(pi);

  si->set_starbase(sbpi);
  sbpi->set_node(si);
 
  //  single_star * ss = dynamic_cast(single_star*, si->get_starbase());
  update_binary_parameters(bi, final3);
}

// Exchange incoming 'bullet' star in binary, leaving the primary intact; 
// replacing the secondary star
void exchange_primary(dyn *bi,
		      dyn *si,
                      final_state3& final3) {

  //      print_result(target);

  //  cerr<<"exchange_primary"<<endl;

  //  dyn *root = bi->get_parent();
  dyn *pi = bi->get_oldest_daughter()->get_younger_sister();
  starbase *sbpi = (starbase*)pi->get_starbase();
  starbase *sbsi = (starbase*)si->get_starbase();
  pi->set_starbase(sbsi);
  sbsi->set_node(pi);

  si->set_starbase(sbpi);
  sbpi->set_node(si);
 
  //  single_star * ss = dynamic_cast(single_star*, si->get_starbase());

  update_binary_parameters(bi, final3);
}


// merge two stars and put the result in the first argument (*pi) of the tree
void merge_two_stars(dyn *pi, dyn *si) {

  //  cerr << "merge two stars" << endl;

  starbase *spi = pi->get_starbase();
  starbase *ssi = si->get_starbase();
  starbase *sbf;
  if (spi->get_total_mass()>ssi->get_total_mass()) {
    sbf = spi->merge_elements(dynamic_cast(single_star*, ssi));
    //    sbf = spi;
  }
  else {
    sbf = ssi->merge_elements(dynamic_cast(single_star*, spi));
    //    sbf = ssi;
  }

  // insert the merger product in the tree
  pi->set_starbase(sbf);
  sbf->set_node(pi);
}

// Merge incoming 'bullet' star with the binarie's secondary, 
// leaving the primary intact and member of the 'new' binary. 
void merger_binary_primary(dyn* bi,
			   dyn *si,
			   final_state3& final3) {

  cout<<"merger_binary_primary"<<endl;

  dyn *pi = bi->get_oldest_daughter()->get_younger_sister();
  merge_two_stars(pi, si);

  // detach the single star
  detach_node_from_general_tree(si);
  update_binary_parameters(bi, final3);
}

// Merge incoming 'bullet' star with the binarie's primary
// leaving the secondary intact and member of the 'new' binary. 
void merger_binary_secondary(dyn* bi,
			     dyn *si,
			     final_state3& final3) {

  //  cout<<"merger_binary_secondary"<<endl;

  dyn *pi = bi->get_oldest_daughter();

  merge_two_stars(pi, si);
  update_binary_parameters(bi, final3);
  //  cerr << "Binary parameters updated"<<endl;
}

void merger_binary_bullet(dyn* bi,
			  dyn* si,
                          final_state3& final3) {

  //  cerr<<"start: merger_binary_bullet"<<endl;
  dyn *pip = bi->get_oldest_daughter();
  dyn *pis = pip->get_younger_sister();

  merge_two_stars(pip, pis);

  // insert the merger product in the tree
  pis->set_starbase(si->get_starbase());
  si->get_starbase()->set_node(pis);

  // detach the single star
  detach_node_from_general_tree(si);

  update_binary_parameters(bi, final3);
}

void merger_escape_primary(dyn *bi, dyn *si, 
			   final_state3 &final3) {

  //  cerr<<"merge_escape_primary"<<endl;

  dyn *pip = bi->get_oldest_daughter();
  dyn *pis = pip->get_younger_sister();

  ionize_binary(bi);
  merge_two_stars(pis, si);

  // detach the single star
  detach_node_from_general_tree(si);
}

void merger_escape_secondary(dyn *bi, dyn *si, 
			     final_state3 &final3) {

  //  cerr<<"merge_escape_secondary"<<endl;

  dyn *pip = bi->get_oldest_daughter();

  ionize_binary(bi);
  merge_two_stars(pip, si);

  // detach the single star
  detach_node_from_general_tree(si);
}

void merger_escape_bullet(dyn *bi, dyn *si, 
			  final_state3 &final3) {

  //  cerr<<"merge_escape_bullet"<<endl;

  dyn *pip = bi->get_oldest_daughter();
  dyn *pis = pip->get_younger_sister();

  ionize_binary(bi);
  merge_two_stars(pip, pis);

  // detach the single star
  detach_node_from_general_tree(pis);
}

void merger_triple(dyn *bi, dyn *si, 
		   final_state3 &final3) {

  dyn *pip = bi->get_oldest_daughter();
  dyn *pis = pip->get_younger_sister();

  ionize_binary(bi);
  merge_two_stars(pip, si);
  detach_node_from_general_tree(si);
  merge_two_stars(pip, pis);
  detach_node_from_general_tree(pis);
}

void stopped_scatter3() {

  cerr << "void stopped_scatter3() from scattter3" << endl;
  cerr << "Discission: don't adjust binary or single star."<<endl;
  cerr << "do perform an extra collision."<<endl;

  //      perform_scatter(target, bullet);
}

void unknown_final_scatter3() {
  cerr << "void unknown_final_scatter3() from scattter3" << endl;
  cerr << "Discission: don't adjust binary or single star."<<endl;
  cerr << "do perform an extra collision."<<endl;

  //      perform_scatter(target, bullet);
}

void error_scatter3() {
  cerr << "void error_scatter3() from scattter3" << endl;
  cerr << "Discission: don't adjust binary or single star."<<endl;
  cerr << "do perform an extra collision."<<endl;

  //      perform_scatter(target, bullet);
}

// Process the result of the interaction. 
// The result is in final_state3::final.
// Routines for processing the data are in the header-file scatter_support.h
final_descriptor3 form_scatter_product(dyn * bi, 
				       dyn * si,
				       initial_state3& init3,
				       intermediate_state3& inter3,
				       final_state3 & final3) { 

  //		Assume the binary is bound (which is tue).
  final_descriptor3 fate = preservation;
  //		Have a look at our single star.
  //		And save it for later use.
  //      star_state sa;
  //      star_state str = present->get_pre_str();
  //      real mss = str.mass;
  //      real rss = str.radius;
  //      real vss = str.velocity;
  //      stellar_type type = str.type;
  //      sa=str;

  // for test purposes:
  //  final3.descriptor = triple_merger;
  switch (final3.descriptor) {
  case preservation      : 
    preserve_initial_state(bi, si, final3);
    //              sa.velocity = vss;
    break;
  case exchange_1        :
    //              sa.make_star_appeal((star*)target->get_primary());
    exchange_primary(bi, si, final3);
    //              sa.velocity = vss;
    break;
  case exchange_2        : 
    //              sa.make_star_appeal((star*)target->get_secondary());
    exchange_secondary(bi, si, final3);
    //              sa.velocity = vss;
    break;
  case merger_binary_1   : 
    //              sa.velocity = -1;
    merger_binary_primary(bi, si, final3);
    break;
  case merger_binary_2   : 
    //              sa.velocity = -2;
    merger_binary_secondary(bi, si, final3);
    break;
  case merger_binary_3   : 
    //              sa.velocity = -3;
    merger_binary_bullet(bi, si, final3);
    break;
  case ionization:	
    ionize_binary(bi); //, mss, vss, type, init3, final3);
    //              sa.velocity = vss;
    //              fate = ionized;
    break;
  case merger_escape_1   : 	// ...
    merger_escape_primary(bi, si, final3);
    //              sa.velocity = 0;
    //              fate = coalided;
    break;
  case merger_escape_2   : 	// ...
    merger_escape_secondary(bi, si, final3);
    //              sa.velocity = 0;
    //              fate = coalided;
    break;
  case merger_escape_3   : 	// ...
    merger_escape_bullet(bi, si, final3);
    //              sa.velocity = vss;
    //              fate = coalided;
    break;
    //             merger_escape(target, mss, vss, type, init3, final3);
    //              sa.velocity = vss;
    //              fate = coalided;
    //              break;
  case triple_merger     : 
    merger_triple(bi, si, final3);
    //              sa.clean();
    //              fate = coalided;
    break;
  case error             : 
    error_scatter3();
    //		Whenever an error occurs in acatter3 the ineraction must
    //		be performed anyway. So just try it again.
    //		Of course, this can give trouble by converging to
    //		a solution in the 3-body model.
    //		perform_scatter(target, present, init3, inter3,
    //                                final3, v_esc);
    break;
  case stopped           : 
    stopped_scatter3();
    break;
  case unknown_final     : 
    unknown_final_scatter3();
    break;
  default                : 
    cout << "Unexpected switch " << endl;
    cout <<"in: form_scatter_product()"<<endl;
  }
  //      
  //         present->set_post_str(sa);
  //         present->update_present(target, inter3, final3);
  //         if(target->get_velocity()>=v_esc)
  //            fate = escaped;
  //         present->set_fate(fate);
  return fate;
}

#endif // TOOLBOX
