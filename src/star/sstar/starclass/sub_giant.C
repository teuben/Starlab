//
// sub_giant.C
//

#include "sub_giant.h"
#include "hertzsprung_gap.h"

// ANSI C++ first creates the base class before the dreived classes are
// created. 
 
     sub_giant::sub_giant(hertzsprung_gap & h) : single_star(h) {

       delete &h; 

// (GN+SPZ May  4 1999) last update age is time of previous type change
      last_update_age = next_update_age;

       adjust_next_update_age();

       instantaneous_element();
       update();

       post_constructor();
}


#if 0
void sub_giant::adjust_initial_star() {

  if(relative_age<=0) {
    real t_ms = main_sequence_time();
    relative_age = Starlab::max(t_ms + hertzsprung_gap_time(t_ms), relative_age);
  }
}
#endif

void sub_giant::instantaneous_element() {
//cerr << "SG: Instant element"<<endl;

  real l_g = giant_luminosity();
  real t_ms = main_sequence_time();
  real tend_hg = t_ms + hertzsprung_gap_time(t_ms);
  real t_gs = 0.15*main_sequence_time();

  real tau = t_gs / (tend_hg + t_gs - relative_age);

//  PRC(tau);PRC(t_gs);PRC(tend_hg);PRC(t_gs);PRL(relative_age);
//  PRC(luminosity);

  luminosity = l_g*pow(tau, 1.17);
  luminosity = Starlab::min(luminosity, maximum_luminosity());

  // (SPZ+GN:  1 Aug 2000)
  // coming from the hertzsprung-gap the effective readius should 
  // remain the same.
  //  effective_radius = 
            radius = (0.25*pow(luminosity, 0.4) + 0.8*pow(luminosity, 0.67))
                   / pow(relative_mass, 0.27);
}

// evolve a sub_giant star upto time argument according to
// the model discribed by Eggleton et al. 1989.
void sub_giant::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {
         real l_g = giant_luminosity();
         real t_ms = main_sequence_time();
         real tend_hg = t_ms + hertzsprung_gap_time(t_ms);
         real t_gs = 0.15*main_sequence_time();

	 real tau = t_gs / (tend_hg + t_gs - relative_age);
	 real tau_prev = t_gs
	              / (tend_hg + t_gs - (relative_age-dt));
	 
         luminosity = l_g*pow(tau, 1.17);
         luminosity = Starlab::min(luminosity, maximum_luminosity());

         radius = (0.25*pow(luminosity, 0.4) + 0.8*pow(luminosity, 0.67))
                / pow(relative_mass, 0.27);

	 // Not a miss C++ contest winner but it might works.
	 // (SPZ+GN:23 Sep 1998)
	 real m_tot = get_total_mass();
	 real new_mcore;
         if (relative_mass < // but not <= see white_dwarf::get_element_type()
	     cnsts.parameters(upper_ZAMS_mass_for_degenerate_core)) {
	   new_mcore = Starlab::max(core_mass, 0.146*pow(luminosity, 0.143));
	 }
	 else {
           // (GN Oct 11 1999) correct for hydrogen fraction envelope: X = 0.7
           real X = cnsts.parameters(hydrogen_fraction);

	   new_mcore = core_mass 
                     + l_g*6*t_gs*(pow(tau,1./6)-pow(tau_prev,1./6.))
	             / (X*cnsts.parameters(energy_to_mass_in_internal_units));
	 }
	 
	 // (GN+SPZ Apr 29 1999) no type change by core growth
	 //(SPZ+GN: 31 Jul 2000)
	 // removed the safety, nothing against type change by core growth.
	 //core_mass = min(m_tot-cnsts.safety(minimum_mass_step), new_mcore);
	 core_mass = Starlab::min(m_tot, new_mcore);
	 envelope_mass = m_tot - core_mass;

				    
      }
      else {
//		sub_giant lifetime exceeded. Transform star into
//		horizontal branch star.
         star_transformation_story(Horizontal_Branch);
         new horizontal_branch(*this);
         return;
      }

      update();
      stellar_wind(dt);
   }

star* sub_giant::reduce_mass(const real mdot) {

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

star* sub_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot);

      if (envelope_mass<=mdot) {
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

void sub_giant::adjust_accretor_age(const real mdot,
				    const bool rejuvenate=true) {

     real t_ms = main_sequence_time();
     real tend_hg_old = hertzsprung_gap_time(t_ms) + t_ms;
     real t_bgb_old = base_giant_branch_time(t_ms);

     real m_tot_new = get_total_mass() + mdot;
     real m_rel_new = Starlab::max(m_tot_new, relative_mass);

     t_ms = main_sequence_time(m_rel_new);
     real tend_hg_new = t_ms + hertzsprung_gap_time(m_rel_new, t_ms);
     real t_bgb_new = base_giant_branch_time(m_rel_new, t_ms);

     real dtime = relative_age - tend_hg_old;

     // For relative_mass > 20.97184... (see hertzsprung_gap.C)
     // sub_giants can not exist. (SPZ+GN:10 Oct 1998)
     
     if (t_bgb_new>0) {

       // (GN+SPZ May  4 1999) update last_update_age
       last_update_age = tend_hg_new;
       relative_age = tend_hg_new 
                    + dtime*(t_bgb_new/t_bgb_old);
       if (rejuvenate) {
	  //	 cerr << "Rejuvinating"<<endl;
	  //         PRC(relative_age);PRC(mdot);PRC(m_tot_new);
	  // 	   PRL(rejuvenation_fraction(mdot/m_tot_new));
	 relative_age *= rejuvenation_fraction(mdot/m_tot_new);
	  //	   PRL(relative_age);
       }

       relative_age = Starlab::max(relative_age, 
			  last_update_age  + cnsts.safety(minimum_timestep));
     }
     else {
       // Relative_age should be set to the predicted next update age.
       // Instead use tend_hg_new, which is end point of HG.
       //       relative_age = next_update_age;
       relative_age = tend_hg_new;
     }

     // next_update_age should not be reset here
     // next_update_age = tend_hg_new+t_bgb_new;

   }

void sub_giant::adjust_next_update_age() {

      real t_ms = main_sequence_time();
      real t_giant = t_ms + hertzsprung_gap_time(t_ms)
                          + base_giant_branch_time(t_ms);
      
      next_update_age = t_giant;
   }

void sub_giant::detect_spectral_features() {

      single_star::detect_spectral_features();


      if (accreted_mass>=cnsts.parameters(B_emission_star_mass_limit))
	spec_type[Emission]=Emission;

}

#if 0
real sub_giant::stellar_radius(const real mass, const real new_age) {

      real t_ms = main_sequence_time(mass);
      real t_giant = t_ms + hertzsprung_gap_time(mass, t_ms)
                          + base_giant_branch_time(mass, t_ms);

      real age = new_age;
      if(new_age<t_giant)
        age = new_age;

      real l_g = giant_luminosity(mass);
      real tend_hg = t_ms + hertzsprung_gap_time(mass, t_ms);
      real t_gs = 0.15*t_ms;

      real l_sg = l_g*pow(t_gs/(tend_hg + t_gs - age), 1.17);
           l_sg = Starlab::min(l_sg, maximum_luminosity(mass));

      real r_sg = (0.25*pow(l_sg, 0.4) + 0.8*pow(l_sg, 0.67))
                / pow(mass, 0.27);

      return r_sg;
   }
#endif

real sub_giant::gyration_radius_sq() {

  return cnsts.parameters(convective_star_gyration_radius_sq); 
}


real sub_giant::zeta_adiabatic() {

// (GN+SPZ Apr 28 1999) fit from Lev Yungelson private communication
// for giants with not too convective envelope = radiative envelope

  real r_dconv = 2.4*pow(relative_mass,1.56);
  if (relative_mass > 10 )
    r_dconv = 5.24*pow(relative_mass,1.32);
  else if (relative_mass > 5)
    r_dconv = 1.33*pow(relative_mass,1.93);
    
  if (radius < r_dconv) {

    return 4;
  }
  else {
//   		Hjellming and Webbink 1987 ApJ, 318, 804
    real x = core_mass/get_total_mass();
    real A = -0.220823;
    real B = -2.84699;
    real C = 32.0344;
    real D = -75.6863;
    real E = 57.8109;

    return A + x*(B + x*(C + x*(D + x*E)));

  }
}

// Values of zeta are changed (SPZ+GN:28 Sep 1998)
real sub_giant::zeta_thermal() {
  //cerr<<"real single_star::zeta_thermal()"<<endl;

  real z;
  if (low_mass_star())
    z = 0;
  else 
    z = 0; // (GN+SPZ Apr 28 1999) radius determined by core only (was -1) 

  return z;
}

// wind_constant is fraction of envelope lost in nuclear lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 1 Oct 1998)
void sub_giant::update_wind_constant() {
  
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
  else { // no wind for low mass ms stars
// (GN Apr 16 1999) 0.2 Msun loss for degenerate core stars
//    real t_ms = main_sequence_time();

    // (SPZ+GN: 26 Jul 2000) see Nelemans, YPZV 2000 (A&A Submitted)
    // wind_constant = 0.2; is a slight improvement on
    // Though a minimum of 0.0 is somewhat on the low side.
    wind_constant = Starlab::max(0., (2.5 - relative_mass)/7.5);

// (GN+SPZ May  4 1999) not needed: single_star::stellar_wind changed
//                  /(1 - 
//		   pow((t_ms + hertzsprung_gap_time(t_ms))
//                      /next_update_age,
//                     cnsts.parameters(massive_star_mass_loss_law)));
  }
  
  wind_constant = Starlab::max(wind_constant, 0.0);

}


