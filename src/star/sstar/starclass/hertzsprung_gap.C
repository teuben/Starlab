#include "hertzsprung_gap.h"
#include "main_sequence.h"

// A low mass (M<25M_\odot) main_sequence star will transform
// into a hertzsprung_gap star after the exhaustion of hydrogen 
// in the core.
// This constructor 
// copies the old main_sequence star into the newly created 
// hertzsprung_gap star and destroys the old main_sequence star 
// and finally evolves the hertzsprung_gap in order to determine its
// appearence.
//
// ANSI C++ first creates the base class before the derived classes are
// created. 

     hertzsprung_gap::hertzsprung_gap(main_sequence & m) : single_star(m) {

      delete &m; 	

      real m_tot    = get_total_mass();
      core_mass     = TAMS_helium_core_mass();
      envelope_mass = m_tot - core_mass;
      core_radius   = helium_core_radius();

//      PRC(current_time);PRC(relative_age);PRL(next_update_age);
// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

      adjust_next_update_age();

      instantaneous_element();
      update();
      
      post_constructor();
   }

void hertzsprung_gap::instantaneous_element() {

      real alpha, beta, gamma, delta;

      real log_mass = log10(relative_mass);

      if (relative_mass>1.334) {
         alpha = 0.1509 + 0.1709*log_mass;
         beta  = 0.06656 - 0.4805*log_mass;
         gamma = 0.007395 + 0.5083*log_mass;
         delta = (0.7388*pow(relative_mass, 1.679)
               - 1.968*pow(relative_mass, 2.887))
               / (1.0 - 1.821*pow(relative_mass, 2.337));
       }
       else {
          alpha = 0.08353 + 0.0565*log_mass;
          beta  = 0.01291 + 0.2226*log_mass;
          gamma = 0.1151 + 0.06267*log_mass;
          delta = pow(relative_mass, 1.25)
                * (0.1148 + 0.8604*relative_mass*relative_mass)
                / (0.04651 + relative_mass*relative_mass);
       }

      real l_bgb  = base_giant_branch_luminosity();
      real rt = delta*pow(10., (alpha + beta + gamma));
	  
      luminosity = l_bgb;
      // (SPZ+GN:  1 Aug 2000)
      // coming from the main sequence the effective readius should 
      // remain the same.
//      effective_radius = radius     = rt;
      radius     = rt;

   }

#if 0
void hertzsprung_gap::adjust_initial_star() {

  if(relative_age<=0)
    relative_age = Starlab::max(main_sequence_time(), relative_age);
}
#endif

// evolve a hertzsprung_gap star upto time argument according to
// the model described by Eggleton et al. 1989.
void hertzsprung_gap::evolve_element(const real end_time) {
//cerr<<"void hertzsprung_gap::evolve_element(T="<<end_time<<")"<<endl;
//cerr<<luminosity<<" "<<radius<<endl;

      real alpha, beta, gamma, delta;

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      real log_mass = log10(relative_mass);

      if (relative_mass>1.334) {
         alpha = 0.1509 + 0.1709*log_mass;
         beta  = 0.06656 - 0.4805*log_mass;
         gamma = 0.007395 + 0.5083*log_mass;
         delta = (0.7388*pow(relative_mass, 1.679)
               - 1.968*pow(relative_mass, 2.887))
               / (1.0 - 1.821*pow(relative_mass, 2.337));
       }
       else {
          alpha = 0.08353 + 0.0565*log_mass;
          beta  = 0.01291 + 0.2226*log_mass;
          gamma = 0.1151 + 0.06267*log_mass;
          delta = pow(relative_mass, 1.25)
                * (0.1148 + 0.8604*relative_mass*relative_mass)
                / (0.04651 + relative_mass*relative_mass);
       }

       real t_ms   = main_sequence_time();
      
       if (relative_age<=next_update_age) {

          if(relative_age<t_ms) relative_age = t_ms; // safety
          real t_hg   = hertzsprung_gap_time(t_ms);
          real l_g    = giant_luminosity();
          real l_bgb  = base_giant_branch_luminosity();

          real rt = delta*pow(10., (alpha + beta + gamma));
          real rg = (0.25*pow(l_g, 0.4) + 0.8*pow(l_g, 0.67))
                  / pow(relative_mass, 0.27);
          real ff = (relative_age - t_ms)/t_hg;
          luminosity = l_bgb*pow(l_g/l_bgb, ff);
          radius = rt*pow(rg/rt, ff);

	  evolve_core_mass(dt);
       }
       else {

	 // Stars with M > 20.97184...Msun do not ascend the sub_giant branch
	 // but become directly a horizontal branch star.
	 // See Tout et al. 1997, implemented (SPZ+GN:23 Sep 1998)
	 if (base_giant_branch_time(t_ms)) {
          star_transformation_story(Sub_Giant);
          new sub_giant(*this);
          return;
	 }
	 else {
          star_transformation_story(Horizontal_Branch);
          new horizontal_branch(*this);
          return;
        }
       }

       update();
       stellar_wind(dt);
}

star* hertzsprung_gap::reduce_mass(const real mdot) {

      if (envelope_mass<=mdot) {
	envelope_mass = 0;

	// (SPZ+GN: 27 Jul 2000)
	// non degenerate core < helium_dwarf_mass_limit always(!) become
	// white dwarfs
	if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
	    relative_mass < cnsts.parameters(
			    upper_ZAMS_mass_for_degenerate_core)) {
	  star_transformation_story(Helium_Dwarf);
	  return dynamic_cast(star*, new white_dwarf(*this));
	} 
	else {
	  star_transformation_story(Helium_Star);
	  return dynamic_cast(star*, new helium_star(*this));
	}
      }

      envelope_mass -= mdot;
      return this;
   }

star* hertzsprung_gap::subtrac_mass_from_donor(const real dt, real& mdot) {
//cerr<<"real hertzsprung_gap::subtrac_mass_from_donor( dt= "<<dt<<")"<<endl;

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();

      mdot = mass_ratio_mdot_limit(mdot);

      if (envelope_mass<=mdot) {
//	  cerr << "Transform hertzsprung_gap star into Heliun stars"<<endl;
         mdot = envelope_mass;
         envelope_mass = 0;

	// (SPZ+GN: 27 Jul 2000)
	// non degenerate core < helium_dwarf_mass_limit always(!) become
	// white dwarfs
	if (get_total_mass() < cnsts.parameters(helium_dwarf_mass_limit) &&
	    relative_mass < cnsts.parameters(
			    upper_ZAMS_mass_for_degenerate_core)) {
	  star_transformation_story(Helium_Dwarf);
	  return dynamic_cast(star*, new white_dwarf(*this));
	} 
	else {
	  star_transformation_story(Helium_Star);
	  return dynamic_cast(star*, new helium_star(*this));
	}
      }

// (GN+SPZ Apr 29 1999) 
      adjust_donor_radius(mdot);

      envelope_mass -= mdot;
      return this;
   }

//  relative age adjusted lineairly with accreted mass.
void hertzsprung_gap::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {
//cerr <<"void hertzsprung_gap::adjust_accretor_age("<<mdot<<")"<<endl;

      real m_rel_new;
      real m_tot_new = get_total_mass() + mdot;
      if (m_tot_new>relative_mass)
         m_rel_new = m_tot_new;
      else m_rel_new = relative_mass;

      real t_ms_old = main_sequence_time();
      real t_hg_old = hertzsprung_gap_time(t_ms_old);

      real t_ms_new = main_sequence_time(m_rel_new);
      real t_hg_new = hertzsprung_gap_time(m_rel_new, t_ms_old);

      real dtime = relative_age - t_ms_old;

// (GN+SPZ May  4 1999) update last_update_age
      last_update_age = t_ms_new;

      relative_age = t_ms_new 
                   + dtime*(t_hg_new/t_hg_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new); 

      relative_age = Starlab::max(relative_age, 
			 last_update_age + cnsts.safety(minimum_timestep)); 
      
      // next_update_age should not be reset here
      // next_update_age = t_ms_new + t_hg_new;
   }

// Adiabatic response function for hertzsprung_gap star.
// Polynomial fit to Hjellming and Webbink 1987 ApJ, 318, 804
real hertzsprung_gap::zeta_adiabatic() {
//cerr<<"real hertzsprung_gap::zeta_adiabatic(): " << endl;

#if 0
      real z;

      real x = core_mass/get_total_mass();
      real A = -0.220823;
      real B = -2.84699;
      real C = 32.0344;
      real D = -75.6863;
      real E = 57.8109;

      if (get_relative_mass()<=0.4)
         z = -cnsts.mathematics(one_third);
      else if (low_mass_star())
         z = A + x*(B + x*(C + x*(D + x*E)));
      else if (medium_mass_star())
         z = 2.25;                 // 15 according to Pols & Marinus 1994
      else                         // We do it differently.
         z = 2.25;                 // lekker puh.
#endif
      real z = 4; // this is neede to prevent Thermal in extreme mass
                     // ratio systems ... 
// (GN+SPZ Apr 29 1999) Pols & Marinus 1994 were maybe right: not!

      return z;
   }

// Thermal response function for hertzsprung_gap star.
real hertzsprung_gap::zeta_thermal() {

      real z;

      if (get_relative_mass()<=0.4)
         z = 0;         // no better estimate present.
      else if (low_mass_star())
         z = -2; 	// -10 according to Pols & Marinus 1994
      else              // Changed to current values
         z = -2;	// by (SPZ+GN: 1 Oct 1998)

      return z;
   }

void hertzsprung_gap::adjust_next_update_age() {
   
      real t_ms = main_sequence_time();
      if (relative_age<t_ms)
 	 relative_age = t_ms;
      next_update_age = t_ms + hertzsprung_gap_time(t_ms);
}

void hertzsprung_gap::detect_spectral_features() {

// 		Use standard spectral feature detection.
      single_star::detect_spectral_features();

      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;
   }

#if 0
real hertzsprung_gap::stellar_radius(const real mass, const real age_max) {

      real alpha, beta, gamma, delta;
      real t_ms   = main_sequence_time(mass);
      real age = Starlab::max(t_ms, age_max);            // Safety
      real t_hg   = hertzsprung_gap_time(t_ms);

      real log_mass = log10(mass);

      if (mass>1.334) {
         alpha = 0.1509 + 0.1709*log_mass;
         beta  = 0.06656 - 0.4805*log_mass;
         gamma = 0.007395 + 0.5083*log_mass;
         delta = (0.7388*pow(mass, 1.679)
               - 1.968*pow(mass, 2.887))
               / (1.0 - 1.821*pow(mass, 2.337));
       }
       else {
          alpha = 0.08353 + 0.0565*log_mass;
          beta  = 0.01291 + 0.2226*log_mass;
          gamma = 0.1151 + 0.06267*log_mass;
          delta = pow(mass, 1.25)
                * (0.1148 + 0.8604*mass*mass)
                / (0.04651 + mass*mass);
       }

       real l_g    = giant_luminosity(mass);

       real rt = delta*pow(10., (alpha + beta + gamma));
       real rg = (0.25*pow(l_g, 0.4) + 0.8*pow(l_g, 0.67))
               / pow(mass, 0.27);
       real ff = (age - t_ms)/t_hg;
       real r_hg = rt*pow(rg/rt, ff);

       return r_hg;
   }
#endif


// Stellar Gyration radii squared for detmination of
// angular momentum.
// Implemented by (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

// Helium core at the terminal-age main-sequence
// (SPZ+GN: 1 Oct 1998)
real hertzsprung_gap::TAMS_helium_core_mass() {

  // (SPZ+GN: 28 Jul 2000)
  // old Eggleton 2000 (ToBeBook) fit to core mass.
  real mc = (0.11*pow(relative_mass,1.2)
	  + 7.e-5*pow(relative_mass,4))
          / (1+ 2e-4*pow(relative_mass, 3));

  // (SPZ+GN: 28 Jul 2000)
  // produces lower mass core betwee 0 and 20 Msun.
  // is fine for low mass (<5Msun) stars but too small for 
  // high mass (>8Msun) stars.
  //  real mc = (0.05 + 0.08*pow(relative_mass,1.2)
  //          + 7.e-5*pow(relative_mass,4))
  //          / (1+ 2e-4*pow(relative_mass, 3));

  // (SPZ+GN: 31 Jul 2000)
  // TAMS core mass must be smaller than get_total_mass to 
  // assure that it becomes a helium star decently.
  return Starlab::min(Starlab::max(mc,core_mass), 
	     get_total_mass()-cnsts.safety(minimum_mass_step));

}

void hertzsprung_gap::evolve_core_mass(const real dt) {

  // (GN Oct 11 1999) correct for hydrogen fraction envelope: X = 0.7
  real X = cnsts.parameters(hydrogen_fraction);
  real delta_mcore = luminosity*dt
                   / (X * cnsts.parameters(energy_to_mass_in_internal_units));

// (GN+SPZ Apr 29 1999) no type change by core growth
  delta_mcore = Starlab::min(delta_mcore, envelope_mass);
  //(SPZ+GN: 31 Jul 2000)
  // removed the safety, nothing against type change by core growth.
  //  delta_mcore = min(delta_mcore, 
  //		    envelope_mass-cnsts.safety(minimum_mass_step));

  if (delta_mcore<0) {
      PRC(luminosity);PRC(dt);PRL(delta_mcore);
      cerr << "Core mass increas is negative in hertzsprung_gap:evolve_core_mass()"<<endl;
      delta_mcore = 0;
      dump(cerr, false);
      exit(-1);
      return;
  }

  envelope_mass -= delta_mcore; 
  core_mass += delta_mcore;
}

// wind_constant is fraction of envelope lost in nuclear lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 1 Oct 1998)
void hertzsprung_gap::update_wind_constant() {
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
  else { // (GN+SPZ Apr 29 1999) 1% loss on hg

    wind_constant = 0.01*relative_mass;
  }

  wind_constant = Starlab::max(wind_constant, 0.0);
}
