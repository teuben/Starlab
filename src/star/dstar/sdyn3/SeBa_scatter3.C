//
// SeBa_scatter3.C:	Evolve a binary with a single scatter event.
//
//// Input is a single binary and a single star.
//// at some time specified time the star and binary encounter eachother.
//// This encounter is resolved with the STARLAB scatter3 package.
//// usage:  SeBa_scatter3
//// options      -d #      Impact parameter in Rsun [random]
////              -v #      velocity at infinity in km/s [random]
////              -t #      time of encounter in Myr [not specified]
////              -T #      end time 
////
////              the initial single star and binary are read in via 
////              standard IO.
////              A suitable input snapshot can be created using the 
////              following STARLAB command line:
//// example for creating initial binary and single star:
//// example: % makenode -n 2 | \
////            makemass -u 13 -l 13 | \
////            makesecondary -q -u .1 -l 0.1 -f 0.5 |
////            add_star -R 1 | \
////            makebinary -f 2 -o 1 -l 1000 -u 1000 -e 0 > infile
////
//// creates a zero-age binary with a 13Msun primay and a 1.3Msun secondary 
//// star together with a 13Msun tertiary.
//// The binary is cirular with an orbital separation of 1000Rsun.
////
//// SeBa_scatter3 is subsequently invoked with:
//// example: % SeBa_scatter3 -t 10 -T 100 -v 10 -d 10 < infile
////

#include "SeBa_tt.h"

#ifdef TOOLBOX

#define REPORT_ADD_DOUBLE false
#define  SEED_STRING_LENGTH  60

//		Relative velocity accordig to 
//		Verbunt, F., \& Meylan, G., 1988, A\&A 203, 297
real relative_velocity(real vt, real vb) {

  return sqrt(pow(vt, 2.) + pow(vb, 2.));
}

real relative_velocity(vec vt, vec vb) {

  return relative_velocity(abs(vt),abs(vb));
}

real relative_velocity(dyn* target, dyn* bullet) {

  return relative_velocity(target->get_vel(), bullet->get_vel());
}

real random_impact_parameter(const real m1,
			     const real m2,
			     const real d_max, 
			     const real v_inf) {
      

//	Maximum  impact parameter at this velocity: v_inf.
//	Apoparently I am not convinced that this is the proper 
//	way to do it.
        real sigma_max = cross_section(m1, m2, d_max, v_inf);
        real random    = randinter(0, 1);
        real sigma_rnd = random*sigma_max;

      return sqrt(sigma_rnd/PI);


}

// Here's all the hairy stuff concerning the scattering event.
// external routine scatter3 is called.
// scatter3 performs the tree body encounter and sets-up the system 
// from scratch.
// all it needs to know is the masses of the stars, their radii
// and their relative velocity.
local void perform_scatter(dyn *bi,
			   dyn *si,
                           initial_state3& initial3,
                           intermediate_state3& inter3,
                           final_state3& final3,
			   real dt_snap) { 

  double_star *binary = (double_star*)bi->get_starbase();
  single_star *single = (single_star*)si->get_starbase();
  

//	Setup the effective size and mass of the single star.
  real mss = single->get_total_mass();
  real rss = single->get_effective_radius();
  real vss = single->get_velocity();
  PRC(mss);PRC(rss);PRL(vss);

//	Setup the effective size and mass of the binary.
     real m_bin = binary->get_total_mass();
     real v_bin = binary->get_velocity();	
     real r_bin = binary->get_effective_radius();
     real e_bin = binary->get_eccentricity();

//		Scale sizes and masses the way scatter3 wants it.
//		If we are nice to scatter3, scatter3 is noce to us.
//      real v_crit = critical_velocity(binary, mss);
     real v_crit = 1;

      real r_scale = 1./binary->get_semi();
      real m_scale = 1./binary->get_total_mass();
      real v_scale = 1./v_crit;

      real v_rel  = relative_velocity(v_bin, vss);
      //      real v_rel_esc = sqrt(2.)*v_esc;
//      real d_max = 2*(r_bin + rss);
//      real d_max = maximum_closest_approach(m_bin, r_bin, e_bin,
//					    mss, rss, v_rel);
      real d_max = 3/r_scale;
      

      initial3.r1  = r_scale*binary->get_primary()->get_effective_radius();
      initial3.r2  = r_scale*binary->get_secondary()->get_effective_radius();
      initial3.r3  = r_scale*rss;
      initial3.m2  = m_scale*binary->get_secondary()->get_total_mass();
      initial3.m3  = m_scale*mss;
      initial3.sma = r_scale*binary->get_semi();
      real v_esc = 10*v_rel;
      real vinf_rnd  = random_focussed_maxwellian_velocity(m_bin, mss, d_max,
                                                           v_rel, v_esc);
      //      initial3.v_inf = v_scale*vinf_rnd;
      //      initial3.rho = r_scale*sqrt(1./r_scale);
      //      initial3.rho = r_scale*random_impact_parameter(m_bin, mss, 
      //      initial3.rho = r_scale*random_impact_parameter(m_bin, mss, 
      //      						     d_max, vinf_rnd);
      real rho_max = r_scale*sqrt(cross_section(m_bin, mss, d_max, vinf_rnd)/PI);

      initial3.ecc = binary->get_eccentricity();

//		Error output for testing cases.

cerr<<"scatter3 -x "<<initial3.r1<<" -y "<<initial3.r2<<" -z "<<initial3.r3
    <<" -m "<<initial3.m2<<" -M "<<initial3.m3<<" -e "
    <<initial3.ecc<<" -r "<<initial3.rho<<" -v "<<initial3.v_inf
    <<" -s "<<get_rand_seed() <<endl;

      unsigned rand_seed = get_rand_seed();

//		Be sure to randomize the conditions of the
//		scatetring event.
      real cpu_time_check = 3600;
      real dt_out = VERY_LARGE_NUMBER;
      randomize_angles(initial3.phase);
      scatter3(initial3, inter3, final3,
	       cpu_time_check, dt_out, dt_snap);

//		Still nice to see what happens on the screen.
//      cout << "result:  " << state_string(inter3.descriptor) << " "
//           << state_string(final3.descriptor) << endl;


//		Print it all to the standard output.
      print_initial(cout, initial3);
      print_intermediate(cout, inter3);
      print_final(cout, final3);
//      form_scatter_product(binary, present, initial3, inter3, final3);

//		Keep track of errors made by scatter3.
      if (final3.descriptor==error) {
         ofstream error("error.dat", ios::app|ios::out);
         if(!error) cerr << "error: couldn't crate file error.dat" <<endl;
         error<<"scatter3 -x "<<initial3.r1<<" -y "<<initial3.r2
              <<" -z "<<initial3.r3 <<" -m "<<initial3.m2
              <<" -M "<<initial3.m3<<" -e " <<initial3.ecc
              <<" -r "<<initial3.rho<<" -v "<<initial3.v_inf
              <<" -s "<<rand_seed <<endl;
         error.close();
      }
      else {
         ofstream scttr3("scttr3.dat", ios::app|ios::out);
         if(!scttr3) cerr << "error: couldn't crate file scttr3.dat" <<endl;
         scttr3 << binary->get_primary()->get_element_type() << " " 
                << binary->get_secondary()->get_element_type() << " " 
                << single->get_element_type() << endl;
         scttr3 << m_scale << " " << r_scale << " " <<initial3.r1<<" "
                <<initial3.r2 <<" "<<initial3.r3 <<" "<<initial3.m2 <<" "
                <<initial3.m3<<" " <<initial3.ecc <<" "<<initial3.rho<<" "
                <<initial3.v_inf <<" "<<rand_seed <<endl;
         scttr3 << inter3.n_osc << " " << inter3.r_min_min << " " 
                << inter3.descriptor << endl;
         scttr3 << final3.sma << " " << final3.ecc << " " 
                << final3.descriptor << endl;
         if (initial3.rho>=0.9*rho_max) {
            if (initial3.rho>=rho_max) {
               cerr << "ERROR: rho = " << initial3.rho << " > rho_max = "
	            << rho_max << endl;
	    }
         ofstream outer90("outer90_enc.dat", ios::app|ios::out);
         if(!outer90) cerr << "error: couldn't crate file outer90.dat" <<endl;
         outer90 << binary->get_primary()->get_element_type() << " " 
                << binary->get_secondary()->get_element_type() << " " 
                << single->get_element_type() << endl;
         outer90 << m_scale << " " << r_scale << " " <<initial3.r1<<" "
                <<initial3.r2 <<" "<<initial3.r3 <<" "<<initial3.m2 <<" "
                <<initial3.m3<<" " <<initial3.ecc <<" "<<initial3.rho<<" "
                <<initial3.v_inf <<" "<<rand_seed <<endl;
         outer90 << inter3.n_osc << " " << inter3.r_min_min << " " 
                << inter3.descriptor << endl;
         outer90 << final3.sma << " " << final3.ecc << " " 
                << final3.descriptor << endl;
         }
/*
//		Additional output, which we do not want to know.
         scttr3 << "init:" << endl;
         scttr3 << initial3.system[0].mass<<" "<<initial3.system[0].vel << endl;
         scttr3 << initial3.system[1].mass<<" "<<initial3.system[1].vel << endl;
         scttr3 << initial3.system[2].mass<<" "<<initial3.system[2].vel << endl;
         scttr3 << "final:" << endl;
         scttr3 << final3.system[0].mass << " " << final3.system[0].vel << endl; 
         scttr3 << final3.system[1].mass << " " << final3.system[1].vel << endl;
         scttr3 << final3.system[2].mass << " " << final3.system[2].vel << endl;
*/
         scttr3.close();
      }

}

real critical_velocity(const real a, const real m1, const real m2, 
		       const real m3) {

  real GM_R = cnsts.physics(G)*cnsts.parameters(Msun)
	    / cnsts.parameters(Rsun);
  real v_crit = sqrt(GM_R*m1*m2*(m1+m2+m3) / (a*m3*(m1+m2))); //Km_s;

  return v_crit/cnsts.physics(km_per_s);
}

real critical_velocity(dyn* bi, dyn* si) {

      real GM_R = cnsts.physics(G)*cnsts.parameters(Msun)
	        / cnsts.parameters(Rsun);
      real a  = bi->get_starbase()->get_semi();
      real m1 = ((double_star*)bi->get_starbase())->get_primary()->get_total_mass();
      real m2 = ((double_star*)bi->get_starbase())->get_secondary()->get_total_mass();

      real m3  = si->get_starbase()->get_total_mass();

      return critical_velocity(a, m1, m2, m3);
}

dyn * get_binary(dyn *root, int index=-1) {

  for_all_nodes(dyn, root, bi) {
    if (bi->is_parent() && !bi->is_root()) {
      if(index<0 || bi->get_index()==index)
	return bi;
    }
  }
  
  cerr << "No binary found in star_cluster::get_binary()" << endl;
  exit(-1);
}

dyn * get_single(dyn *root, int index=-1) {

  for_all_leaves(dyn, root, bi) {
    if (!(bi->is_parent() && !bi->is_root())) {
      if(index<0 || bi->get_index()==index)
	return bi;
    }
  }
  
  cerr << "No single star found in star_cluster::get_single()" << endl;
  exit(-1);
}

// Finding the apprpriate elements that for an encounter is done
// with a probability lineair to their cross section.
// First cose from all particles.
// secondly chose from the single (excluding possibly first chosen
// element) stars only.
// Finally perform the encounter.
// So, actually this routine does it all.
local final_descriptor3 collide_elements(dyn * root,
					 real v_inf, 
					 real impact_parameter,
					 real dt_snap) {

      initial_state3 initial3;
      make_standard_init(initial3);

      dyn* bi = get_binary(root);
      dyn* si = get_single(root);

      real v_crit = critical_velocity(bi, si);
      initial3.v_inf = v_inf/v_crit;
      initial3.rho = impact_parameter/bi->get_starbase()->get_semi();

      intermediate_state3 inter3;
      final_state3 final3;

      perform_scatter(bi, si, initial3, inter3, final3, dt_snap);
      return form_scatter_product(bi, si,
                                  initial3, inter3, final3); 
      
   }

local void do_the_adding(dyn *bi, real stellar_time,
			 binary_type type,
			 real sma, real ecc)
{
  dyn* od = bi->get_oldest_daughter();

  if (stellar_time < 0) {
    cerr << "No dyn time, implement hdyn" << endl;
    exit(-1);
  }

  if (REPORT_ADD_DOUBLE) {
    int p = cerr.precision(HIGH_PRECISION);
    cerr << "Before adddouble: at " << stellar_time << " to: ";
    bi->pretty_print_tree(cerr);
    od->get_kepler()->print_all(cerr);
    cerr.precision(p);
  }

  // addstar(bi, stellar_time); // Should not be needed.

  int id;
  if((od->is_low_level_node() &&
      od->get_younger_sister()->is_low_level_node()) &&
     (od->get_elder_sister() == NULL) &&
     bi->n_leaves()==2 ) {

    story * old_story = bi->get_starbase()->get_star_story();
    bi->get_starbase()->set_star_story(NULL);

    id = bi->get_index();

    double_star* new_double
      = new_double_star(bi, sma, ecc, stellar_time, id, type);
    // synchronize_binary_components(dynamic_cast(double_star*,
    // 						  bi->get_starbase()));

    if (REPORT_ADD_DOUBLE) {
      cerr<<"Adddouble: id, t, a, e: "
	  <<new_double->get_identity()<<" "
	  <<new_double->get_current_time() << " "
	  <<new_double->get_semi()<<" "
	  <<new_double->get_eccentricity()<<endl;
      put_state(make_state(new_double));
    }

    // Give the new binary the old star_story.

    new_double->set_star_story(old_story);

    if (REPORT_ADD_DOUBLE) {
      int p = cerr.precision(HIGH_PRECISION);
      cerr<<"After adddouble: at "<<stellar_time<< " to: ";
      new_double->dump(cerr);
      od->get_kepler()->print_all(cerr);
      put_state(make_state(new_double));
      cerr << endl;
      cerr.precision(p);
    }
  }
}

local void adddouble(dyn *b,
		     real dyn_time, 
		     binary_type type) {

  if (REPORT_ADD_DOUBLE) {
    int p = cerr.precision(HIGH_PRECISION);
    cerr<<"adddouble: "<<b<<" at T="<<dyn_time<<endl;
    cerr.precision(p);
  }

  real stellar_time = b->get_starbase()->conv_t_dyn_to_star(dyn_time);
  addstar(b, stellar_time);
  // switch off hdyn to include SeBa supernova treatment.
  b->get_starbase()->set_use_hdyn(false);

  real ecc, sma;
  real binary_age=0;
  int id;

  for_all_nodes(dyn, b, bi) {
    if (bi->is_parent() && !bi->is_root()) {

      story * old_story = b->get_starbase()->get_star_story();
      b->get_starbase()->set_star_story(NULL);
      binary_age = Starlab::max(bi->get_oldest_daughter()->get_starbase()
				->get_current_time(),
				bi->get_oldest_daughter()->get_binary_sister()
				->get_starbase()->get_current_time());

      id = bi->get_index();

      if (!bi->get_oldest_daughter()->get_kepler()) {
	new_child_kepler(bi, dyn_time, MAX_CIRC_ECC);
	//cerr << "Kepler added" << endl;
      }
      else {
	//cerr << "Kepler already present...?" << endl;
	bi->get_oldest_daughter()->get_kepler()->	  // Probably
	  set_circular_binary_limit(MAX_CIRC_ECC);  // unnecessary
	dyn_to_child_kepler(bi, dyn_time);
      }
		
      bi->get_oldest_daughter()->get_kepler()->print_elements();

      sma = b->get_starbase()->
	conv_r_dyn_to_star(bi->get_oldest_daughter()->
			   get_kepler()->get_semi_major_axis());
      ecc = bi->get_oldest_daughter()->get_kepler()->
	get_eccentricity();
      PRC(sma);PRL(ecc);

      do_the_adding(bi, binary_age, type, sma, ecc);

    }
  }
}

local bool  evolve_binary_until_next_time(dyn * bi,
					  real end_time,
			  bool stop_at_merger_or_disruption = false,
			  bool stop_at_remnant_formation = false) { 

  
  real dt; 

  if (!bi->is_root() &&
      bi->get_parent()->is_root()) {

    // only check once for double_star address (contrary to single star).
    double_star* ds = dynamic_cast(double_star*, 
				   bi->get_starbase());
    ds->dump("SeBa.data", true);
    
    real time = ds->get_current_time();
    do {

      dt = ds->get_evolve_timestep() + cnsts.safety(minimum_timestep);
      time = Starlab::min(time+dt, end_time);

      ds->evolve_element(time);

      //      cerr << "Binary type (t="<<time<<"): " << type_string(ds->get_bin_type()) << endl;

      if (stop_at_merger_or_disruption &&
	  (ds->get_bin_type() == Merged || 
	   ds->get_bin_type() == Disrupted))
	return false;

      if (stop_at_remnant_formation &&
	 (ds->get_primary()->remnant() || ds->get_secondary()->remnant()))
	return false;

    }
    while (time<end_time);

    ds->dump("SeBa.data", true);
  }

  return true;
}

#define EPSILON 1.e-10

local void evolve_star_until_next_time(dyn* bi, const real end_time) {

  single_star* ss = dynamic_cast(single_star*, bi->get_starbase());

  //  real current_time = ss->get_current_time();
  real dt;
  real time = ss->get_current_time();
  do {

    // have to check for address change due to stellar evolution
    ss = dynamic_cast(single_star*, bi->get_starbase());

    dt = ss->get_evolve_timestep();
    time = Starlab::min(time+dt, end_time);

    ss->evolve_element(time);
    
    bool verbose = false;
    if(verbose) {
      star_state ss(dynamic_cast(star*, bi->get_starbase()));
      put_state(ss, cerr);
    //             print_star(bi->get_starbase(), cerr);
      int p = cerr.precision(HIGH_PRECISION);
      bi->get_starbase()->dump(cerr, false);
      cerr.precision(p);
    }

  }
  while (time<end_time);

  //  bi->get_starbase()->evolve_element(out_time);
  //  print_star(bi->get_starbase(), cerr);
}

local void  local_evolve_the_stellar_system(dyn* b, real time) {

  dyn *bi;

  for_all_daughters(dyn, b, bi) {
    if (bi->is_parent() && !bi->is_root()) {
      evolve_binary_until_next_time(bi, time);
    }
    else {
      evolve_star_until_next_time(bi, time);
    }
  }
}

void evolve_element(dyn *root, real t_end) {

  local_evolve_the_stellar_system(root, t_end);
  //  put_dyn(cerr, *root);
  root->set_system_time(t_end);
  repair_tree(root);
}

void evolve_stars_and_binaries(dyn * root, real t_end, real dt_out) {

  do {
    real time = Starlab::min(t_end, 
			     dynamic_cast(real, 
					  root->get_system_time()+dt_out));
    evolve_element(root, time);
  }
  while(root->get_system_time()<t_end);
}

bool binary_still_exists(dyn *root) {

  bool binary_ok = false;
  for_all_nodes(dyn, root, bi) {
    if (bi->is_parent() && !bi->is_root()) 
      binary_ok = true;
  }

  return binary_ok;
}

/*-----------------------------------------------------------------------------
 *  main  --
 *      usage:
 *              scatterev many options.  ,
 *
 *              where # is the initial age of the cluster.
 *      options:
 *              The following options are allowed:
 *      cluster age:
 *              -t #    Where # stands for the initial age of the
 *                      in Myear.
 *
 *              At present the running time of the integrator correspnds
 *              to the stellar age an a one by 10e6year basis.
 *              This however should be scaled to the cluster parameters.
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
  int  c;
  bool c_flag = false;

  real  t_start = 0;     
  real  t_end  = 12000;
  real  t_scatter = -1;   // 
  real  v_enc = -1;            // [km/s]
  real  impact_parameter = -1; // [Rsun]
  real  dt_snap = VERY_LARGE_NUMBER;

  int   input_seed = 0;
  char  *comment;

  //  int srandinter(int);
  char  seedlog[SEED_STRING_LENGTH];
  extern char *poptarg;
  int  pgetopt(int, char **, char *);

  int n_output = 1;
  while ((c = pgetopt(argc, argv, "c:D:d:n:s:T:t:v:")) != -1)
    switch(c)
      {
      case 'c': c_flag = TRUE;
	comment = poptarg;
	break;
      case 'D': dt_snap = atof(poptarg);
	break;
      case 'd': impact_parameter = (profile)atoi(poptarg);
	break;
      case 'n': n_output = atoi(poptarg);
	break;
      case 's': input_seed = atoi(poptarg);
	break;
      case 'T': t_end = atof(poptarg);
	break;
      case 't': t_scatter = atof(poptarg);
	break;
      case 'v': v_enc = atof(poptarg);
	break;
      case '?': cerr <<
		  "usage: addstar [-T \"..\"] [-c \"..\"] "
		     << "for relative age\n";
	exit(1);
      }

  dyn *root = get_dyn(cin);
  adddouble(root, t_start, Detached);
  // Switch stellar evolution on
  root->set_use_sstar(true);
  // Switch internal (N-body) dynamics off
  root->get_starbase()->set_use_hdyn(false);

  char paramlog[255];
  sprintf(paramlog, 
	  "         alpha  = %3.1f\n         lambda = %3.1f\n         beta   = %3.1f\n         gamma  = %4.2f",
	  cnsts.parameters(common_envelope_efficiency),
	  cnsts.parameters(envelope_binding_energy),
	  cnsts.parameters(specific_angular_momentum_loss),
	  cnsts.parameters(dynamic_mass_transfer_gamma));

  root->log_comment(paramlog);
  int actual_seed = srandinter(input_seed);
  cerr << "random number generator seed = " << actual_seed << endl;
  root->log_history(argc, argv);
  sprintf(seedlog, "       random number generator seed = %d",actual_seed);
  root->log_comment(seedlog);
  if (c_flag == TRUE)
    root->log_comment(comment);

  bool mass_spectrum = true;
  vec center = 0;
  bool verbose = true;
  dstar_stats(root, mass_spectrum, center, verbose);

  real dt_out = t_end/n_output;
  if(t_scatter>0 && t_scatter<t_end) {
    evolve_stars_and_binaries(root, t_scatter, dt_out);
    dstar_stats(root, mass_spectrum, center, verbose);
  
    if(binary_still_exists(root)) {
      final_descriptor3 fate = preservation;
      fate = collide_elements(root, v_enc, impact_parameter, dt_snap);
    }
  }

  evolve_stars_and_binaries(root, t_end, dt_out);
  dstar_stats(root, mass_spectrum, center, verbose);

  put_dyn(root, cerr);
}


#endif  // TOOLBOX


/* endof: mkstarfield.C */





