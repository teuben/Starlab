//
// horizontal_branch.C
//

#include "horizontal_branch.h"
#include "sub_giant.h"
#include "hertzsprung_gap.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.
//

horizontal_branch::horizontal_branch(sub_giant & g) : single_star(g) {

    delete &g;

// (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;

    adjust_next_update_age();

    instantaneous_element();
    update();

    post_constructor();
}

horizontal_branch::horizontal_branch(hertzsprung_gap & h) : single_star(h) {

    delete &h;

// (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;

    adjust_next_update_age();

    instantaneous_element();
    update();

    post_constructor();
      
}

#if 0
void horizontal_branch::adjust_initial_star() {

  if(relative_age<=0) {
    real t_ms = main_sequence_time();
    real t_giant = t_ms + hertzsprung_gap_time(t_ms)
                 + base_giant_branch_time(t_ms);
    relative_age = max(t_giant, relative_age);
  }
}
#endif


void horizontal_branch::instantaneous_element() {

    luminosity = min(helium_giant_luminosity(), 
		     maximum_luminosity());

    radius = (0.25*pow(luminosity, 0.4)
           + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);

    real t_ms = main_sequence_time();
    real t_giant = t_ms + hertzsprung_gap_time(t_ms)
                 + base_giant_branch_time(t_ms);
    real t_he = helium_giant_time(t_ms);
  
    real ff = (relative_age - t_giant)/t_he;
    if (radius<=25)
      radius *= 1.0 - 0.1*ff;
    else
      radius *= pow(25.0/radius, ff); // Giant

    // (SPZ+GN:  1 Aug 2000)
    // coming from previous type the effective readius should 
    // remain the same.
    //    effective_radius = radius;
}

// evolve a horizontal_branch star upto time argument according to
// the model discribed by Eggleton et al. 1989.
void horizontal_branch::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {
         luminosity = min(helium_giant_luminosity(), 
                          maximum_luminosity());

	 evolve_core_mass(dt);
	 
         radius = (0.25*pow(luminosity, 0.4)
                + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);

         real t_ms = main_sequence_time();
         real t_giant = t_ms + hertzsprung_gap_time(t_ms)
                             + base_giant_branch_time(t_ms);
         real t_he = helium_giant_time(t_ms);

         real ff = (relative_age - t_giant)/t_he;
         if (radius<=25)
            radius *= 1.0 - 0.1*ff;
         else
            radius *= pow(25.0/radius, ff); // Giant
      }
      else {
	// Core mass remains unchanged: second dredge-up 
	// accounted for in Super_Giant (SPZ+GN: 27 Jul 2000)
	star_transformation_story(Super_Giant);
	new super_giant(*this);
         return;
      }

      update();
      stellar_wind(dt);
   }

// Post-evolution stellar update.
// Makes sure age and radius are updated.
void horizontal_branch::update() {

  // New core mass determination occurs in ::evolve_element.
  // (SPZ+GN:09/1998)
  // real m_tot = get_total_mass();
  // core_mass = helium_core_mass();
  // envelope_mass = m_tot - core_mass;

  core_radius = helium_core_radius();

// (GN+SPZ Apr 28 1999)
// effective_radius can not be larger than radius for horizontal branch stars
  effective_radius = radius;

  detect_spectral_features();

}

star* horizontal_branch::reduce_mass(const real mdot) {

      if (envelope_mass<=mdot) {
         envelope_mass = 0;
         star_transformation_story(Helium_Star);
         return dynamic_cast(star*, new helium_star(*this));
      }

      envelope_mass -= mdot;
      return this;
   }

star* horizontal_branch::subtrac_mass_from_donor(const real dt,
						 real& mdot) {

      real mdot_temp = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_temp);

      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;
         star_transformation_story(Helium_Star);
         return dynamic_cast(star*, new helium_star(*this));
      }

      envelope_mass -= mdot;
      return this;
   }

void horizontal_branch::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {

      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_ms = main_sequence_time();
      real t_giant_old = t_ms + hertzsprung_gap_time(t_ms) 
                              + base_giant_branch_time(t_ms);
      real t_he_old = helium_giant_time(t_ms);

      t_ms = main_sequence_time(m_rel_new);
      real t_giant_new = t_ms + hertzsprung_gap_time(m_rel_new, t_ms) 
                              + base_giant_branch_time(m_rel_new, t_ms);
      real t_he_new = helium_giant_time(m_rel_new, t_ms); 

      real dtime = relative_age - t_giant_old;

// (GN+SPZ May  4 1999) update last_update_age
      last_update_age = t_giant_new;
      relative_age = t_giant_new
                   + dtime*(t_he_new/t_he_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new);

      relative_age = max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep));

      // next_update_age should not be reset here
      // next_update_age = t_giant_new + t_he_new;

   }

real horizontal_branch::zeta_adiabatic() {
      return 15;	
   }

real horizontal_branch::zeta_thermal() {
      return 15;
   }

void horizontal_branch::adjust_next_update_age() {
     
      real t_ms = main_sequence_time();
      real t_giant = t_ms + hertzsprung_gap_time(t_ms)
                          + base_giant_branch_time(t_ms);
      real t_he = helium_giant_time(t_ms);
      next_update_age = t_giant + t_he;
   }

real horizontal_branch::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

void horizontal_branch::evolve_core_mass(const real dt) {

  real X = cnsts.parameters(hydrogen_fraction);
  // (SPZ+GN: 27 Jul 2000) Factor 0.5 introduced to model
  // that the luminosity is partially generated by core burning.
  // see Nelemans YPZV 2000
  real delta_mcore = 0.5 * luminosity*dt 
                   / (X * cnsts.parameters(energy_to_mass_in_internal_units));

  //(SPZ+GN: 31 Jul 2000)
  // removed the safety, nothing against type change by core growth.
  //  delta_mcore = min(delta_mcore, 
  //		    envelope_mass-cnsts.safety(minimum_mass_step));
  delta_mcore = min(delta_mcore, envelope_mass);

  if (delta_mcore<0) {
      PRC(luminosity);PRC(dt);PRL(delta_mcore);
      cerr << "Core mass increas is negative in horizontal_branch::evolve_core_mass()"<<endl;
      delta_mcore = 0;
      exit(-1);
//      return;
  }

  envelope_mass -= delta_mcore; 
  core_mass += delta_mcore;
      
}

// wind_constant is fraction of envelope lost in nuclear lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 1 Oct 1998)
void horizontal_branch::update_wind_constant() {
  
// (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_constant = meader_fit_dm;
    else {// constant
      real final_mass = 30;
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // 5% of initial envelope

    wind_constant = 0.05*(relative_mass - final_core_mass());
  }
  
  wind_constant = max(wind_constant, 0.0);

}
