//
// helium_giant.C
//

#include "super_giant.h"
#include "helium_star.h"
#include "helium_giant.h"


// (GN+SPZ May  3 1999) only for stars more massive than 
// super_giant2neutron_star, below become white_dwarf
// to be done (SPZ+GN: 27 Jul 2000)
helium_giant::helium_giant(super_giant & g) : single_star(g) {

    delete &g;

// (GN+SPZ May  4 1999) last update age is time of previous type change
  last_update_age = next_update_age;

    real second_dredge_up_time = next_update_age 
                          * Starlab::min(1., relative_mass
			  / cnsts.parameters(super_giant2neutron_star));

    real remaining_time = second_dredge_up_time - relative_age;
    next_update_age = helium_time();
    relative_age = next_update_age - remaining_time;
    last_update_age = relative_age;

    // (SPZ+GN: 27 Jul 2000)
    // core mass has been set in super_giant before construcor is called.

    instantaneous_element();
    update();
    
    post_constructor();
 }

helium_giant::helium_giant(helium_star & h) : single_star(h) {

    delete &h;
    
// (GN+SPZ May  4 1999) last update age is time of previous type change
    last_update_age = next_update_age;

    adjust_next_update_age();

    update_wind_constant();
    
    instantaneous_element();
    update();

    post_constructor();

}

void helium_giant::instantaneous_element() {

  real m_tot = get_total_mass();

  real temp = log10(m_tot)
            * (0.4509 - 0.1085*log10(m_tot)) + 4.7143;
  luminosity  = 8.33*temp - 36.8;
  temp = 1.e-3*pow(10., temp);
  luminosity  = pow(10., luminosity);
  effective_radius = radius = 33.45*sqrt(luminosity)/(temp*temp);
  core_radius = 0.2*radius;
  
  if (envelope_mass > 0) {
    if (m_tot>=2.7)
      radius = 12./get_total_mass();
    else if (m_tot>=2.0) 
      radius *= 250;
    else if (m_tot>=1.0) 
      radius *= 150;
    else radius *= 6;
  }
  else {
    effective_radius = radius = core_radius;
  }
}

void helium_giant::evolve_element(const real end_time) {

        real dt = end_time - current_time;
        current_time = end_time;
        relative_age += dt;
  
        real m_tot = get_total_mass();

        if (relative_age<=next_update_age) {
            real tmp = log10(m_tot)
                * (0.4509 - 0.1085*log10(m_tot)) + 4.7143;
            luminosity  = pow(10., 8.33*tmp - 36.8);
            tmp = pow(1.e-3*pow(10., tmp), 2);
            radius = 33.45*sqrt(luminosity)/tmp;
            core_radius = 0.2*radius;

//		Helium giant branch.
//		Imitate: Habets, GMHJ, 1986, A&A 167, 61.
            if (envelope_mass > cnsts.safety(minimum_mass_step)) {
                  if (m_tot>=2.7)
                     radius = 12./get_total_mass();
                  else if (m_tot>=2.0) 
                     radius *= 250;
                  else if (m_tot>=1.0) 
                     radius *= 150;
                  else radius *= 6;
            }
            else {
               radius = core_radius;
            }
         }
         else {
            stellar_wind(dt);
            create_remnant();
            return;
         }

         update();
         stellar_wind(dt);

}

void helium_giant::update() {

  real m_tot = get_total_mass();
  core_mass = COcore_mass = CO_core_mass();
  envelope_mass = m_tot - core_mass;

  // removed (SPZ+GN:10 Nov 1998)
  //core_radius = helium_core_radius();

  detect_spectral_features();
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;
  effective_radius = radius;
}

// should only be used in constructors.
// (SPZ+GN:30 Sep 1998)
void helium_giant::adjust_next_update_age() {

  next_update_age /= cnsts.parameters(helium_star_lifetime_fraction);

}


void helium_giant::create_remnant() {

     if (is_binary_component()) 
       get_binary()->dump("binev.data", false);
     else
       dump("binev.data", false);

        stellar_type type = NAS;
     //if (get_total_mass() >= cnsts.parameters(helium2black_hole)) 
//           type = Black_Hole;

     if (relative_mass >= cnsts.parameters(maximum_main_sequence) &&
	 relative_mass < 300)
       type = Disintegrated;
     else if (core_mass >= cnsts.parameters(COcore2black_hole)) 
       type = Black_Hole;
     else if(core_mass >= cnsts.parameters(Chandrasekar_mass))
       type = Neutron_Star;
     else
       type = Carbon_Dwarf;

#if 0 //(SPZ+GN: 27 Jul 2000)
     // (GN+SPZ May  3 1999) core mass > 1.4 or total mass > 2.2 
        else if(core_mass >= cnsts. //(GN+SPZ May  3 1999) was 1.3 * cnsts
		             parameters(Chandrasekar_mass) ||
		get_total_mass() >= cnsts.parameters(helium2neutron_star))
	  type = Neutron_Star;
	
	// else if(get_total_mass()>=
	//         1.3*cnsts.parameters(kanonical_neutron_star_mass) &&
	//         relative_mass>=8) 
	// type = Disintegrated;
        else
           type = Carbon_Dwarf;
#endif // (SPZ+GN: 27 Jul 2000)

	switch (type) {
	  case Black_Hole : star_transformation_story(Black_Hole);
                            new black_hole(*this); 
                            break;
          case Neutron_Star : star_transformation_story(Neutron_Star);
			    new neutron_star(*this);
                            break;
          case Disintegrated : star_transformation_story(Disintegrated);
			       new disintegrated(*this);
			       break;
          case Carbon_Dwarf : star_transformation_story(Carbon_Dwarf);
	                      new white_dwarf(*this);
                              break;
          default :   cerr << "helium_star::create_remnant()" <<endl;
                      cerr << "star_type not recognized." << endl;
       }

     }


star* helium_giant::subtrac_mass_from_donor(const real dt, real& mdot) {

      mdot = relative_mass*dt/get_binary()->get_donor_timescale();

      mdot = mass_ratio_mdot_limit(mdot);
      
      if (mdot<envelope_mass)
	 envelope_mass -= mdot;
      else {
	 mdot = envelope_mass;
	 envelope_mass = 0;
	
	 if (core_mass <= cnsts.parameters(minimum_helium_star)) {
	    star_transformation_story(Carbon_Dwarf);
	    return dynamic_cast(star*, new white_dwarf(*this));
	 }
      }

      return this;
}

star* helium_giant::reduce_mass(const real mdot) {

      if (mdot < envelope_mass)
	envelope_mass -= mdot;
      else {
	real rest_mdot = mdot;
	rest_mdot -= envelope_mass;
	envelope_mass = 0;
	   
	if (core_mass>rest_mdot) {
	  if (core_mass-rest_mdot<=
	      cnsts.parameters(minimum_helium_star)) {
	    core_mass -= rest_mdot;
	    COcore_mass = core_mass;
	  }
	  else {
	    core_mass -= rest_mdot;
	    COcore_mass = core_mass;
	  }
	}
	else {
	  cerr << "ERROR!:"<<endl;
	  cerr << "void helium_giant::reduce_mass(mdot="
	       << rest_mdot<<")"<<endl;
	  cerr << "mdot exceeds helium core mass ("<<core_mass
	       << ")"<<endl;
	  cerr << "Decision: Disintegrate helium star!"<<endl;
	  
	  star_transformation_story(Disintegrated);
	  return dynamic_cast(star*, new disintegrated(*this));
	}
      }

      return this;
}


real helium_giant::add_mass_to_accretor(const real mdot) {

        if (mdot<0) {

           cerr << "helium_star::add_mass_to_accretor(mdot=" 
		<< mdot << ")"<<endl;
           cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

	   return 0;

        }

        adjust_accretor_age(mdot);
        envelope_mass += mdot;
	relative_mass = Starlab::max(relative_mass, get_total_mass());

	// next_update_age should nog be altered here (SPZ+GN: 3 Oct 1998)

// (GN+SPZ May  5 1999) wind_constant must be zero-age envelope mass, 
// not current
//	update_wind_constant();
	wind_constant += mdot;

	set_spec_type(Accreting);
	
        return mdot;

     }

real helium_giant::add_mass_to_accretor(real mdot, const real dt) {

        if (mdot<0) {

	  cerr << "helium_star::add_mass_to_accretor(mdot=" 
	       << mdot << ")" <<endl;
	  cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
        }

        mdot = accretion_limit(mdot, dt);
        adjust_accretor_age(mdot);
        envelope_mass += mdot;
// (GN+SPZ May  5 1999) wind_constant must be zero-age envelope mass, 
// not current
//	update_wind_constant();
	wind_constant += mdot;
	relative_mass = Starlab::max(relative_mass, get_total_mass());

	set_spec_type(Accreting);
	
        return mdot;

     }

real helium_giant::accretion_limit(const real mdot, const real dt) {
       
     real mdot_limit = mdot;

     real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;
     if (mdot>=eddington)
       mdot_limit =  eddington;

     if (envelope_mass<=0)
       mdot_limit = 0;

     return mdot_limit;

}

void helium_giant::adjust_accretor_age(const real mdot,
				       const bool rejuvenate) {

      real frac = (1-pow(mdot/(get_total_mass()+mdot),
			     cnsts.parameters(rejuvenation_exponent)));
      last_update_age *= frac;
      relative_age *= frac;

     }

// see helium_star.C
real helium_giant::zeta_adiabatic() {

     real z = 0;
//      Hjellming and Webbink 1987 ApJ, 318, 804
     real x = core_mass/get_total_mass();
     real A = -0.220823;
     real B = -2.84699;
     real C = 32.0344;
     real D = -75.6863;
     real E = 57.8109;

     if (get_total_mass()<=0.4)
        z = -cnsts.mathematics(one_third);
     else
        z = A + x*(B + x*(C + x*(D + x*E)));

     return z;

   }

real helium_giant::zeta_thermal() {

    real z = -2;

    if (get_core_mass()<=0.4)  // like a white dwarf
	z = -cnsts.mathematics(one_third);

    return z;

}

real helium_giant::CO_core_mass() {
      // C/O core of helium star grows linearly with time
      // (SPZ+GN:26 Sep 1998)

// (GN+SPZ May  4 1999) core groth totally in helium_star phase: core constant
//  real m_core = get_total_mass()*relative_age/next_update_age;
//  m_core = max(core_mass,m_core);
//    return min(m_core, get_total_mass());

  return Starlab::min(core_mass, get_total_mass());
}

void helium_giant::stellar_wind(const real dt) {

//  PRL(last_update_age);
//  PRL(next_update_age);
//  PRL(relative_age);
//  PRL(previous.relative_age);
//  PRL(dt);
// (GN+SPZ Apr 28 1999) wind for low mass stars per phase
    real end_time = next_update_age - last_update_age;
    real prev_rel_time = Starlab::max(0.,previous.relative_age - last_update_age);
    real relative_time = Starlab::min(relative_age - last_update_age, end_time);

//    PRL(end_time);
//    PRL(relative_time);
//    PRL(prev_rel_time);
//    PRL(wind_constant);
    
    real wind_mass = wind_constant 
                   * (pow(relative_time/end_time,
			cnsts.parameters(massive_star_mass_loss_law))
	           -  pow((relative_time-dt)/end_time,
			cnsts.parameters(massive_star_mass_loss_law)));

// (GN+SPZ May  6 1999) try low wind: 
//    wind_mass = 0.;

//    PRL(wind_mass);
//    PRL(envelope_mass);

  if (wind_mass>=envelope_mass) {
    wind_mass = envelope_mass;
    effective_radius = radius = core_radius;
  }

  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else
    reduce_mass(wind_mass);
  return;
}

real helium_giant::gyration_radius_sq() {

  return cnsts.parameters(radiative_star_gyration_radius_sq); 
}

// (GN+SPZ May  3 1999) helium_giants loose complete envelope in wind
void helium_giant::update_wind_constant() {

//  wind_constant = (1 - cnsts.parameters(helium_star_final_core_fraction))
//                * get_total_mass(); 

// (GN+SPZ May  7 1999) envelope is about 30% of total mass,
// we loose 10% of total mass ....
  wind_constant = 0.3*envelope_mass;

}

stellar_type helium_giant::get_element_type() {

     stellar_type type = Helium_Giant;
     if (envelope_mass <= 0)
	 type = Carbon_Star;

     return type;
     }

real helium_giant::temperature() {

  real T_eff = cnsts.parameters(Tsun)
             * sqrt(sqrt(luminosity)/effective_radius);
  
  return T_eff;

  //return sqrt(33.45*sqrt(luminosity)/effective_radius);
}


bool helium_giant::giant_star() {

  if (envelope_mass > cnsts.safety(minimum_mass_step))
    return TRUE;
  else
    return FALSE;

}

bool helium_giant::remnant() {

  if (envelope_mass > cnsts.safety(minimum_mass_step))
    return FALSE;
  else
    return TRUE;

}









