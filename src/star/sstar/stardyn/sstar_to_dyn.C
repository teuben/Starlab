//
// sstar_to_dyn.C 
//
// Communication between the dynamics and the single stellar evolution.

#include "sstar_to_dyn.h"
#include "star.h"
#include "dyn.h"
//#define MASS_UPDATE_LIMIT 0.01

#define DEBUG false

bool has_sstar(dyn * bi) {

  if(bi->is_leaf() && bi->get_starbase()->get_element_type()!=NAS) {
    return true;
  }else{
    return false;
  }
}

vec conv_v_star_to_dyn(vec& v, real rf, real tf) {

//              Internal velocity is km/s
//              Internal size is solar raii
//              Internal time is Myr.
//              km/s to Rsun/Myr

      real to_Rsun_Myr = cnsts.physics(km_per_s) * cnsts.physics(Myear)
	               / cnsts.parameters(solar_radius);
      real to_dyn      = rf/tf;
      //cerr<<"scaling: "<<to_dyn<<endl;
      //cerr <<"kick: " << v << " dyn_V: " << to_Rsun_Myr * to_dyn * v << endl;
      return to_Rsun_Myr * to_dyn * v;
}

vec anomalous_velocity(dyn* b) {

    //cerr << "in anomalous_velocity..." << endl;
    //PRL(b->get_starbase()->get_anomal_velocity());

     vec anomal_velocity = b->get_starbase()->get_anomal_velocity();
     vec zero_vector=0;
     b->get_starbase()->set_anomal_velocity(zero_vector);

     vec new_velocity = conv_v_star_to_dyn(anomal_velocity,
                       b->get_starbase()->conv_r_star_to_dyn(1),
                       b->get_starbase()->conv_t_star_to_dyn(1));

    //PRL(b->get_starbase()->get_anomal_velocity());

     return new_velocity;
}

real get_effective_radius(dyn *b) {
     return b->get_starbase()->conv_r_star_to_dyn( 
            b->get_starbase()->get_effective_radius());
}

real get_total_mass(dyn *b) {
      return b->get_starbase()->conv_m_star_to_dyn(
             b->get_starbase()->get_total_mass());
}

stellar_type get_element_type(dyn *b) {
      return (stellar_type)b->get_starbase()->get_element_type();
}

#define EPSILON 1.e-10

local void evolve_the_stars(node* bi, const real end_time) {

  if(DEBUG) 
      cerr << "evolve_the_stars: " <<bi->format_label() <<" "<<end_time;
  
    real current_time = ((star*)bi->get_starbase())->get_current_time();
    real time_step    =  bi->get_starbase()->get_evolve_timestep();
    
    // PRC(end_time); PRC(current_time); PRC(time_step);

    // SPZ 20 Aug 2004: Just trying this, as it did not seem to work properly...

    while (end_time>current_time+time_step) {

      //  while (end_time>current_time) {
      //      cerr << "running loop at " << time_step << " for star "
      //           << bi->get_index() << endl;

       bi->get_starbase()->evolve_element(current_time+time_step);
       bi->get_starbase()->evolve_element(
             Starlab::min(current_time+time_step+EPSILON, end_time));
       current_time = ((star*)bi->get_starbase())->get_current_time();
       time_step    =  bi->get_starbase()->get_evolve_timestep();

      if(DEBUG) 
	  cerr <<" -dt="<<current_time<<" "<<time_step<<endl;
  }

  // bi->get_starbase()->evolve_element(end_time);

  // bi->get_starbase()->get_seba_counters()->step_sstar++;

    if(DEBUG) 
	cerr<<" and leave"<<endl;
}

#undef EPSILON

bool stellar_evolution(dyn *b)
{
    bool update_all_masses = false;
    
    real stellar_evolution_time = b->get_starbase()->
		conv_t_dyn_to_star(b->get_system_time());

    // Special treatment of black-hole "radii" implemented by Steve (1/05).
    // For efficiency in the integrator, it is useful to assign a negative
    // radius to each black hole.  The value returned by get_radius() will
    // always be positive, but member functions that "need to know" can
    // test the value directly and take appropriate action.
    
    for_all_leaves(dyn, b, bi) {
	if (! bi->is_low_level_node()
	    || !((star*)(bi->get_starbase()))->is_binary_component()) {

	  starbase *sb = bi->get_starbase();

	  if (DEBUG) ((star*)sb)->dump(cerr);

//	  cerr << "Time = " << bi->get_system_time() << " "
//	       << bi->format_label() << endl;
//	  PRC(get_total_mass(bi) - bi->get_mass());

	  real old_dyn_mass = bi->get_mass();
	  evolve_the_stars(bi, stellar_evolution_time);
//	  sb->evolve_element(stellar_evolution_time);

	  real new_dyn_mass_from_star = get_total_mass(bi);
	  real new_radius =  b->get_starbase()->
	      conv_r_star_to_dyn(sb->get_effective_radius());
	  if (old_dyn_mass <= 0 || new_dyn_mass_from_star <= 0) {
	      PRC(old_dyn_mass);PRC(new_dyn_mass_from_star);PRL(new_radius);
	      ((star*)sb)->dump(cerr);
	  }

	  // Special treatment of black hole dynamical radius.

	  if (sb->get_element_type() == Black_Hole) new_radius *= -1;
	  bi->set_radius(new_radius);

//	  cerr << bi->get_system_time() << " "
//	       << bi->format_label() << endl;
//	  PRC(get_total_mass(bi) - bi->get_mass());

	  if (DEBUG) ((star*)sb)->dump(cerr);
 
	  if (abs(new_dyn_mass_from_star-old_dyn_mass)/old_dyn_mass
	      >=cnsts.star_to_dyn(stellar_mass_update_limit)) {
	      update_all_masses = true;
	  }

//	  if (abs(new_dyn_mass_from_star-old_dyn_mass)/old_dyn_mass
//		>=MASS_UPDATE_LIMIT)
//	      update_all_masses = true;
      }
    }

    return update_all_masses;
}

real sudden_mass_loss(dyn* stellar_node) {

   real dm_sun =  stellar_node->get_starbase()->sudden_mass_loss();
   real dm_dyn = stellar_node->get_starbase()->conv_m_star_to_dyn(dm_sun);

   if (dm_sun!=0 && dm_dyn!=0)
     cerr << "sudden mass loss form single star "
	  << stellar_node->format_label()
	  << " dm = " << dm_dyn << " (" << dm_sun << " [Msun]" << endl;

   return dm_dyn;
}


// Normally the primary (in mass) coalesces with the secondary.
// However, this can be rather dangerous if the primary is a giant
// and the total mass of the merger product is likely to exceed that
// for a Wolf-Rayet star.
// Here we check whether or not this exception occurs.
// Normally completely unimportant. But for 30-Doradus it can very well
// happen. (SPZ:02/1998).

bool merge_with_primary(star* primary, star *secondary) {

  bool merge_with_primary = true;
  
  real m_conserved = primary->get_total_mass() 
      + secondary->get_total_mass();

  if (primary->get_total_mass() < secondary->get_total_mass()) {

    merge_with_primary = false;

    if (m_conserved >= cnsts.parameters(massive_star_mass_limit) &&
	!(secondary->get_element_type()==Main_Sequence ||
	  secondary->remnant())) {

      merge_with_primary = true;   // Merge with primary anyway

      if (!(primary->get_element_type()!=Main_Sequence ||
	    primary->remnant()))
      
	cerr << "Merge secondary with primary, although secondary is a "
	     << type_string(secondary->get_element_type())
	     << " and the total mass exceeds that for a Wolf-Rayet."
	     << endl;
    }

  } else {

    if (m_conserved >= cnsts.parameters(massive_star_mass_limit) &&
	!(primary->get_element_type()==Main_Sequence ||
	  primary->remnant())) {
      
      if (secondary->get_element_type()==Main_Sequence ||
	  secondary->remnant())
	
	  merge_with_primary = false;   // Merge with secondary.
      
      else
	  cerr << "Merge primary with secondary, although primary is a "
	       << type_string(primary->get_element_type())
	       << " and the total mass exceeds that for a Wolf-Rayet."
	       << endl;
    }
  }

  return merge_with_primary;
}

