//
// super_giant.C
//

#include "super_giant.h"
#include "horizontal_branch.h"

// ANSI C++ first creates the base class before the dreived classes are
// created.

super_giant::super_giant(horizontal_branch & h) : single_star(h) {

    delete &h;

// (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;

    adjust_next_update_age();

    // (SPZ+GN: 27 Jul 2000) To 'fudge' in the second dredge-up.
    COcore_mass = initial_CO_core_mass(core_mass);
//    PRL(COcore_mass);
    second_dredge_up_time = last_update_age 
                          + (next_update_age-last_update_age) 
                          * min(1., relative_mass
			  / cnsts.parameters(super_giant2neutron_star));
//    PRL(second_dredge_up_time);

    instantaneous_element(); update();

    post_constructor();
}

#if 0
void super_giant::adjust_initial_star() {

  if(relative_age<=0) {
    real t_ms = main_sequence_time();
    real t_giant = t_ms + hertzsprung_gap_time(t_ms)
      + base_giant_branch_time(t_ms);
    real t_he = helium_giant_time(t_ms);
    relative_age = max(t_giant + t_he, relative_age);
   }
}
#endif

void super_giant::instantaneous_element() {

  real l_g  = giant_luminosity();
  real t_ms = main_sequence_time();
  real t_gs = 0.15*t_ms;
  real t_b  = base_giant_time(t_ms);

  if (next_update_age + t_b - relative_age<=0)
    luminosity = l_g*pow(t_gs/t_b, 1.17);
  else
    luminosity = l_g*pow(t_gs/(next_update_age 
			     + t_b - relative_age), 1.17);

  luminosity = min(luminosity, maximum_luminosity());

  // (SPZ+GN:  1 Aug 2000)
  // coming from previous type the effective readius should 
  // remain the same.
  //  effective_radius =	
  radius = (0.25*pow(luminosity, 0.4)
	 + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);
}

// evolve a super_giant star upto time argument according to
// the model discribed by Eggleton et al. 1989.
void super_giant::evolve_element(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {
         real l_g = giant_luminosity();
         real t_ms = main_sequence_time();
         real t_gs = 0.15*t_ms;
         real t_b  = base_giant_time(t_ms);

	 real tau = t_gs
	          / (next_update_age + t_b - relative_age);
	 real tau_prev = t_gs
	              / (next_update_age + t_b - (relative_age-dt));
		
         luminosity = l_g*pow(tau, 1.17);
         luminosity = min(luminosity, maximum_luminosity());

         radius = (0.25*pow(luminosity, 0.4)
                + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);

	 real m_tot = get_total_mass();

	 real new_mcore;
	 if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
	    core_mass     >= cnsts.parameters(helium2neutron_star)) {
	   // Helium core mass does not increase but CO core increases.
	     new_mcore = core_mass;
//	     PRC(second_dredge_up_time);PRL(relative_age);

	     real dmdt  = (core_mass - COcore_mass)
	                / (second_dredge_up_time-(relative_age-dt));
	     COcore_mass += dmdt * dt;
//	     PRC(dmdt);PRL(COcore_mass);
	 }
	 else {
	   if(relative_age <= second_dredge_up_time) {
	   // Helium core mass does not increase but CO core increases.
	     new_mcore = core_mass;
	     // PRC(relative_age);PRC(second_dredge_up_time);
	     // PRC(new_mcore);PRL(core_mass);
	     real dmdt  = (core_mass - COcore_mass)
	                / (second_dredge_up_time-(relative_age-dt));
	     COcore_mass += dmdt * dt;
//	     PRC(dmdt);PRL(COcore_mass);

	   }
	   else {
	     if (luminosity < 15725.) {                     // Mc < 0.73
	       new_mcore = sqrt(luminosity/47488. + 0.1804)+0.015;
	     }
	     else {
	       // See single_star::final_core_mass()
	       new_mcore = luminosity/(46818*pow(relative_mass, 0.25)) + 0.46;
	     }
	   }
	 }

	 //(SPZ+GN: 31 Jul 2000)
	 // removed the safety, nothing against type change by core growth.
	 // core_mass = min(m_tot-cnsts.safety(minimum_mass_step), new_mcore);
	 core_mass = min(m_tot, new_mcore);
	 envelope_mass = m_tot - core_mass;

#if 0 // (SPZ+GN: 27 Jul 2000)
	 // core mass for high mass stars
	 if (relative_mass > cnsts.parameters(super_giant2neutron_star)) { 
 	   real X = cnsts.parameters(hydrogen_fraction);
	   new_mcore = core_mass 
                     + l_g*6*t_gs*(pow(tau,1./6.)-pow(tau_prev,1./6.))
	             / (X*cnsts.parameters(energy_to_mass_in_internal_units));
	 }
	 else {  // Groenewegen & de Jong 1993, A&A 267,410 
	   if (luminosity < 15725.) {                     // Mc < 0.73
	     new_mcore = sqrt(luminosity/47488. + 0.1804)+0.015;
	   }
	   else {
	     // See single_star::final_core_mass()
	     new_mcore = luminosity/(46818*pow(relative_mass, 0.25)) + 0.46;
	   }
	 }

	 new_mcore = max(new_mcore, core_mass);
	 
	 core_mass = min(m_tot-cnsts.safety(minimum_mass_step), new_mcore);
	 envelope_mass = m_tot - core_mass;
#endif

      }
      else {
         create_remnant();
         return;
      }
  
      update();
      stellar_wind(dt);
}

real super_giant::initial_CO_core_mass(const real initial_mass) {

//  return final_CO_core_mass(initial_mass)
//         * cnsts.parameters(helium_star_lifetime_fraction);

  // Implementation of Nelemans YPZV 2000 (A&A submitted)
  // Itroduces discontinuity at relative_mass = 2.2
  // bases on Habets 1986 & IT85
  // (SPZ+GN: 27 Jul 2000)
  real final_coremass_fraction;
  if(initial_mass <= 0.8) 
    final_coremass_fraction = 1;
  else if(initial_mass >= cnsts.parameters(helium2neutron_star)) 
    final_coremass_fraction = 0.65;
  else 
    final_coremass_fraction = 1 - 0.32 * (initial_mass - 0.8);

  return cnsts.parameters(helium_star_lifetime_fraction)
         * final_coremass_fraction*initial_mass;
}

#if 0
void super_giant::evolve_without_wind(const real end_time) {

      real dt = end_time - current_time;
      current_time = end_time;
      relative_age += dt;

      if (relative_age<=next_update_age) {
         real l_g = giant_luminosity();
         real t_ms = main_sequence_time();
         real t_gs = 0.15*t_ms;
         real t_b  = base_giant_time(t_ms);

         luminosity = l_g*pow(t_gs/(next_update_age
                    + t_b - relative_age), 1.17);
         real l_max = maximum_luminosity();
         if(luminosity>l_max) luminosity = l_max;

         radius = (0.25*pow(luminosity, 0.4)
                + 0.8*pow(luminosity, 0.67))/pow(relative_mass, 0.27);
      }
      else {
         create_remnant();
         return;
      }

      update();
   }
#endif

#if 0

void super_giant::stellar_wind(const real dt) {

      real kappa = wind_constant;
      real wind_mass = 1.e6*dt*pow(kappa, 1.43)/pow(10., 13.6);

//              No stellar wind induced helium star creation allowed.
      if (wind_mass>envelope_mass) 
         wind_mass = 0.9*envelope_mass;

      if (is_binary_component())
         get_binary()->adjust_binary_after_wind_loss(this, 
                     wind_mass, dt);
      else
         reduce_mass(wind_mass);
   }
#endif

#if 0
// Determine size of stellar core according to
// Iben, I, Jr., \& Tutukov, A.V., 1993, ApJ 418, 343. Eq. 5.
// and
// Beveren, van D., 1980, PhD thesis, Free University Brussels.
// and
// Iben, I., Jr., and Tutukov, A., V., 1985, Ap. J. Suppl 58, 661.
real super_giant::helium_core_mass() {

//		Iben and Tutukov 1993
      real m_core = (log10(luminosity) - 1.55)/3.56;
//	 	Beveren 1980
      real m1_core = 0.073*(1 + cnsts.parameters(core_overshoot))
	           * pow(relative_mass, 1.42);
//	 	Iben and Tutukov 1985
      real m2_core = 0.058*(1 + cnsts.parameters(core_overshoot))
	           * pow(relative_mass, 1.57);

      m_core = max(max(m_core, m1_core), m2_core);

// 		Dewey, R.J., Cordes, J.M., 1987, ApJ 321, 780.
//		Used a minimum helium core mass of \sim 3 M_\odot.
//  		But this is not realistic.
// 		if(m_core>=2.8 && m_core<3.1) m_core = 3.1;

//              Limit core_mass by realistic limits.
      if (m_core>=cnsts.parameters(helium2neutron_star) &&
	  m_core<cnsts.parameters(minimum_helium_star))
	m_core = m_core<cnsts.parameters(minimum_helium_star);
      
      m_core = min(m_core, get_total_mass());

      return max(m_core, core_mass);
   }
#endif

void super_giant::create_remnant() {

     if (is_binary_component()) 
       get_binary()->dump("binev.data", false);
     else
       dump("binev.data", false);

     real COcore_mass = 0.65 * core_mass;
     stellar_type type;
     if (COcore_mass >= cnsts.parameters(COcore2black_hole)) 
       type = Black_Hole;
     else if(core_mass >= cnsts.parameters(Chandrasekar_mass))
       type = Neutron_Star;
     else
       type = Carbon_Dwarf;

     switch (type) {
       case Black_Hole : star_transformation_story(Black_Hole);
                         new black_hole(*this); 
			 return;
       case Neutron_Star : star_transformation_story(Neutron_Star);
		           new neutron_star(*this);
			   return;
       case Disintegrated : star_transformation_story(Disintegrated);
		            new disintegrated(*this);
			    return;
       case Carbon_Dwarf : star_transformation_story(Carbon_Dwarf);
	                   new white_dwarf(*this);
			    return;
       default :   cerr << "super_giant::create_remnant()" <<endl;
                   cerr << "star_type not recognized." << endl;
                   exit(-1);
     }


#if 0 // (SPZ+GN: 27 Jul 2000)
	if (relative_mass >= cnsts.parameters(super_giant2black_hole) ||
	    core_mass     >= cnsts.parameters(helium2black_hole)) {
	
	star_transformation_story(Black_Hole);
	new black_hole(*this);
	return;
      }
      else if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
	      core_mass     >= cnsts.parameters(Chandrasekhar_mass)) {
	// core_mass is CO core and not Helium core since (SPZ+GN: 27 Jul 2000)
	// core_mass     >= cnsts.parameters(helium2neutron_star)) {
	
	star_transformation_story(Neutron_Star);
	new neutron_star(*this);
	return;
      }
      else if(cnsts.parameters(super_giant_disintegration)) {
	if (relative_mass >= 6 &&
	    relative_mass <= cnsts.parameters(super_giant2neutron_star)) {

	  star_transformation_story(Disintegrated);
	  new disintegrated(*this);
	  return;
	}
      }
      else {
	if (core_mass<=cnsts.parameters(helium_dwarf_mass_limit)) {
	  cerr << "ERROR in super_giant::create_remnant()" << endl;
	  cerr << "Supergiants shoud not become helum dwarfs" << endl;
	  star_transformation_story(Helium_Dwarf);
	}
	else 
	  star_transformation_story(Carbon_Dwarf);

	new white_dwarf(*this);
	return;
      }
#endif // (SPZ+GN: 27 Jul 2000)

}

star* super_giant::reduce_mass(const real mdot) {

  if (envelope_mass<=mdot) {
	envelope_mass = 0;

// (GN+SPZ Apr 29 1999) stripped super_giants become 
// white_dwarf or helium stars
	 if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
	    core_mass     >= cnsts.parameters(helium2neutron_star)) {

	   // (SPZ+GN: 27 Jul 2000)
	   // Initialize core_mass as CO core mass and envelope mass
	   // as helium core mass to make super giant ready to become
	   // a helium giant.
	   real m_tot = core_mass;
	   core_mass = COcore_mass;
	   envelope_mass = m_tot - core_mass;

	   star_transformation_story(Helium_Giant);
	   return dynamic_cast(star*, new helium_giant(*this));
	 }
	 else {

	   // (SPZ+GN: 27 Jul 2000)
	   if(relative_age <= second_dredge_up_time) {

	   real m_tot = core_mass;
	   core_mass = COcore_mass;
	   envelope_mass = m_tot - core_mass;

	     star_transformation_story(Helium_Giant);
	     return dynamic_cast(star*, new helium_giant(*this));
	   }
	   else {
	     star_transformation_story(Carbon_Dwarf);	   
	     return dynamic_cast(star*, new white_dwarf(*this));
	   }
	 }
  }

  envelope_mass -= mdot;
  return this;
}

star* super_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      real mdot_temp = relative_mass*dt/get_binary()->get_donor_timescale();
      mdot = mass_ratio_mdot_limit(mdot_temp);

      if (envelope_mass<=mdot) {
         mdot = envelope_mass;
         envelope_mass = 0;
	 
// (GN+SPZ Apr 29 1999) stripped super_giants become 
// white_dwarf or helium stars
	 if(relative_mass >= cnsts.parameters(super_giant2neutron_star) ||
	    core_mass     >= cnsts.parameters(helium2neutron_star)) {

	   real m_tot = core_mass;
	   core_mass = COcore_mass;
	   envelope_mass = m_tot - core_mass;

	   star_transformation_story(Helium_Giant);
	   return dynamic_cast(star*, new helium_giant(*this));
	 }
	 else {
	   
	   // (SPZ+GN: 27 Jul 2000)
	   if(relative_age <= second_dredge_up_time) {

	   real m_tot = core_mass;
	   core_mass = COcore_mass;
	   envelope_mass = m_tot - core_mass;

	     star_transformation_story(Helium_Giant);
	     return dynamic_cast(star*, new helium_giant(*this));
	   }
	   else {
	     star_transformation_story(Carbon_Dwarf);	   
	     return dynamic_cast(star*, new white_dwarf(*this));
	   }
	 }

      }

// (GN+SPZ Apr 29 1999)
      adjust_donor_radius(mdot);

      envelope_mass -= mdot;
      return this;
}


void super_giant::adjust_accretor_age(const real mdot, const bool rejuvenate=true) {

      real m_tot_new = get_total_mass() + mdot;
      real m_rel_new = max(m_tot_new, relative_mass);

      real t_ms = main_sequence_time();
      real t_hg = hertzsprung_gap_time(t_ms);
      real t_bgb = base_giant_branch_time(t_ms);
      real t_he = helium_giant_time(t_ms);
      real t_super_old = t_ms + t_hg + t_bgb + t_he;
      real t_nuc = nucleair_evolution_time();
      real t_sg_old = t_nuc - t_super_old; 

           t_ms = main_sequence_time(m_rel_new);
           t_hg = hertzsprung_gap_time(m_rel_new, t_ms);
           t_bgb = base_giant_branch_time(m_rel_new, t_ms);
           t_he = helium_giant_time(m_rel_new, t_ms);
      real t_super_new = t_ms + t_hg + t_bgb + t_he;
           t_nuc = nucleair_evolution_time(m_rel_new);
      real t_sg_new = t_nuc - t_super_new; 

      real dtime = relative_age - t_super_old;

// (GN+SPZ May  4 1999) update last_update_age
      last_update_age = t_super_new;
      relative_age = t_super_new
                   + dtime*(t_sg_new/t_sg_old);
      if (rejuvenate)
         relative_age *= rejuvenation_fraction(mdot/m_tot_new);

       relative_age = max(relative_age, 
			  last_update_age + cnsts.safety(minimum_timestep));
      

      // next_update_age should not be reset here
      // next_update_age = t_nuc;
   }

void super_giant::adjust_next_update_age() {

      next_update_age = nucleair_evolution_time();
   }

#if 0
// Condensed polytropes
// see Hjellming and Webbink 1987 ApJ, 318, 804
real super_giant::zeta_adiabatic() {

      real x = core_mass/get_total_mass();
      real a = -0.220823;
      real b = -2.84699;
      real c = 32.0344;
      real d = -75.6863;
      real e = 57.8109;

      real z = a + b*x + c*x*x + d*x*x*x + e*x*x*x*x;

      return z;
   }
#endif

real super_giant::zeta_thermal() {

  real z = 0.; // (GN+SPZ Apr 29 1999) was -10 in 1998; 
               // was -0.64 somewhere in the past (~1992 or so).

      return z;
   }

#if 0
real super_giant::stellar_radius(const real mass, const real age) {

      real t_nuc = nucleair_evolution_time(mass);

      real l_g = giant_luminosity();
      real t_ms = main_sequence_time();
      real t_gs = 0.15*t_ms;
      real t_b  = base_giant_time(t_ms);

      real l_agb = l_g*pow(t_gs/(t_nuc + t_b - age), 1.17);
      l_agb = min(l_agb, maximum_luminosity(mass));

      real r_agb = (0.25*pow(l_agb, 0.4)
             + 0.8*pow(l_agb, 0.67))/pow(mass, 0.27);

      return r_agb;
   }
#endif


real super_giant::gyration_radius_sq() {

  return cnsts.parameters(convective_star_gyration_radius_sq); 
}

void super_giant::update_wind_constant() {
  
// (GN Apr  1 1999) fit for massive stars to Maeder (but not so much wind loss)
// (GN Apr 16 1999) low mass stars need separate treatment

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    wind_constant = meader_fit_dm;

  }
  else {
// (GN+SPZ May  4 1999) nor neede: see single_star::stelar_wind
//    real factor = 1- pow(relative_age/next_update_age,cnsts.parameters(
//                                      massive_star_mass_loss_law));
//    wind_constant = 0.2*(get_total_mass() - final_core_mass())/factor;

    // (SPZ+GN: 27 Jul 2000) 0.8 of inital envelope is lost on AGB
    // see Nelemans YPZV 2000
    wind_constant = 0.8*(relative_mass - final_core_mass());
  }
  
  wind_constant = max(wind_constant, 0.0);

}
