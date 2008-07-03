// 
// single_star.C
//
#include <fstream>
#include "single_star.h"
#include "proto_star.h"
//#include "main_sequence.h"

// Only to make SPZDCH star known to this file 
#include "SPZDCH_star.h"
#include "star_cluster.h"

#define REPORT_MASS_TRANSFER_TIMESCALE true

single_star * new_single_star(stellar_type type,	// All defaults are
			      int  id,			// now specified in
			      real t_cur,
			      real t_rel,
			      real m_rel,
			      real m_tot,
			      real m_core,
			      real co_core,
			      real p_rot,
			      real b_fld,
			      node* n)
{
  cerr << "Initialize from: "<< type_string(type)<<endl;
  PRC(id);PRC(t_cur);PRC(t_rel);PRC(m_rel);PRC(m_tot);
  PRL(m_core); 

  single_star* element = NULL;
  switch(type) {
  case Star_Cluster: element = new star_cluster(n);
    cerr << "new star_cluster in new_single_star()" << endl;
    break;
  case SPZDCH_Star: element = new SPZDCH_star(n);
    break;
  case Proto_Star: element = new proto_star(n);
    break;
  case Planet:
  case Brown_Dwarf: element = new brown_dwarf(n);
    break;
  case Main_Sequence: element = new main_sequence(n);
    break;
  case Hertzsprung_Gap: element = new hertzsprung_gap(n);
    break;
  case Sub_Giant: element = new sub_giant(n);
    break;
  case Horizontal_Branch: element = new horizontal_branch(n);
    break;
  case Super_Giant: element = new super_giant(n);
    break;
  case Carbon_Star:  
  case Helium_Star:  
  case Helium_Giant: element = new helium_star(n);
    break;
  case Hyper_Giant: element = new hyper_giant(n);
    break;
  case Helium_Dwarf:
  case Carbon_Dwarf:
  case Oxygen_Dwarf: element = new white_dwarf(n);
    break;
  case Thorn_Zytkow: element = new thorne_zytkow(n);
    break;
  case Xray_Pulsar:
  case Radio_Pulsar:
  case Neutron_Star: element = new neutron_star(n);
    element->set_rotation_period(p_rot);
    element->set_magnetic_field(b_fld);
    break;
  case Black_Hole: element = new black_hole(n);
    break;
  case Disintegrated: element = new disintegrated(n);
    break;
    //  default: element = new main_sequence(n);
  default: cerr << "No default stellar type in new_single_star()"<<endl;
    exit(1);
		    
  }
      
  element->initialize(id, t_cur, t_rel, m_rel, m_tot, m_core, co_core);

  element->post_constructor();

  return element;
}

single_star::single_star(node* n) : star(n) {

  star_type = NAS;
  for (int i=NAC; i<no_of_spec_type; i++)
    spec_type[i] = NAC;
  current_time=relative_age=0;
  last_update_age = next_update_age=0;
  relative_mass = accreted_mass=0;
  envelope_mass=core_mass=0;
  core_radius=effective_radius=radius=0;
  COcore_mass = 0;
  luminosity=0;
  velocity=0;
  wind_constant=0;
  magnetic_field = rotation_period = 0;
  birth_mass=0;
}


single_star::single_star(single_star & rv) : star(rv) {

  identity      = rv.identity;

  for (int i=NAC; i<no_of_spec_type; i++)
    spec_type[i] = rv.spec_type[i];

  //		copy star
  current_time  = rv.current_time;
  relative_age  = rv.relative_age;
  relative_mass = rv.relative_mass;
  envelope_mass = rv.envelope_mass;
  core_mass     = rv.core_mass;
  COcore_mass   = rv.COcore_mass;

  core_radius   = rv.core_radius;
  radius        = rv.radius;
  luminosity    = rv.luminosity;
  effective_radius = rv.effective_radius;
  last_update_age = rv.last_update_age;
  next_update_age = rv.next_update_age;
  velocity        = rv.velocity;
  wind_constant   = rv.wind_constant;
  accreted_mass   = rv.accreted_mass;
  magnetic_field  = rv.magnetic_field;
  rotation_period = rv.rotation_period;
  birth_mass      = rv.birth_mass;
			
  //              copy stellar history.
  previous.current_time    = rv.previous.current_time;
  previous.last_update_age = rv.previous.last_update_age;
  previous.next_update_age = rv.previous.next_update_age;
  previous.relative_age    = rv.previous.relative_age;
  previous.relative_mass   = rv.previous.relative_mass;
  previous.envelope_mass   = rv.previous.envelope_mass;
  previous.core_mass       = rv.previous.core_mass;
  previous.COcore_mass     = rv.previous.COcore_mass;

  previous.radius          = rv.previous.radius;
  previous.effective_radius= rv.previous.effective_radius;
  previous.magnetic_field  = rv.previous.magnetic_field;
  previous.rotation_period = rv.previous.rotation_period;
  previous.birth_mass      = rv.previous.birth_mass;

}

void single_star::post_constructor() {

  //(GN+SPZ Apr 28 1999) stars with M < 8 Msun have wind per phase:
  update_wind_constant();

  // MEmory refrsh to prevent helium star to become white dwarf with
  // zero-mass core.
  // Not clear yet if companion and binary must also be updated?
  // (SPZ+GN, 10 Nov 1998)
  if (is_binary_component())
    get_binary()->refresh_memory();
  else
    refresh_memory();
  //    refresh_memory();
    
  //    if (is_binary_component() && get_binary()->get_bin_type()==Detached)
  if (is_binary_component() && 
      get_binary()->get_bin_type() != Merged &&
      get_binary()->get_bin_type() != Disrupted) {
    get_binary()->set_bin_type(Detached);
    get_binary()->set_first_contact(false);
  }

  if (is_binary_component() && 
      get_binary()->roche_radius(this) > radius) {
    get_binary()->dump("SeBa.data", true);
  }

  if (remnant()) {
    if (is_binary_component())
      get_binary()->dump("binev.data", false);
    else
      dump("binev.data", false);
  }

}

void single_star::initialize(int id, real t_cur,
			     real t_rel, real m_rel, real m_tot,
			     real m_core, real co_core) {

  //  cerr << "Initialize from: ";
  //  PRC(id);PRC(t_cur);PRC(t_rel);PRC(m_rel);PRC(m_tot);PRL(m_core);
  
  real m_env = m_tot - m_core;
  if (m_rel<=0 || m_rel>cnsts.parameters(maximum_main_sequence)) {
    cerr << "Mass of initial star is prefered to be "
	 << "\nwithin reasonable limits <0, "
	 << cnsts.parameters(maximum_main_sequence)
	 << "] \nbut found (mass = " << m_rel << ")" << endl;
    if (m_rel<=0)
      exit (-1);
  }

  identity = id;
  luminosity = 0.01;
  current_time = t_cur;
  relative_age = t_rel;
  relative_mass = Starlab::max(m_rel, m_env+m_core);
  previous.envelope_mass = envelope_mass = m_env;
  core_mass = m_core;
  COcore_mass = co_core;

  instantaneous_element();

  // (SPZ+GN: 26 Jul 2000) 
  // update_wind_constant() must be called before
  update_wind_constant();
  adjust_next_update_age();

  // Was previously done after adjust_next_update_age, which may be wrong?
  instantaneous_element();

  update();
}

// Determine size of stellar core from fits to
// temperature and luminosity from figure of
// Iben, I.Jr., and Tutukov, A.V., 1985, ApJSS 58, 661.
real single_star::helium_core_radius() {

  // Safety. Should ne be required.
  if (core_mass==0) {
    cerr << "WARNING single_star::helium_core_radius():" << endl;
    cerr << "        core_radius == 0" << endl;
    cerr << "        action: return 0;" << endl;
    return 0;
  }

  real r_he_core;
  real tmp = log10(core_mass)
    * (0.4509 - 0.1085*log10(core_mass)) + 4.7143;
  real lum  = 8.33*tmp - 36.8;
  tmp = 1.e-3*pow(10., tmp);
  lum  = pow(10., lum);

  // helium core is smaller than helium star of same mass.
  // reasoning: after spiral-in class helium_star is supposed to
  // handle further mass transfer.
  //  real fraction = 0.5; // (SPZ+GN: 6 Oct 1998)
  real fraction =1.; // But gives unrealistic mass transfer after spiral-in
  r_he_core = fraction*33.45*sqrt(lum)/pow(tmp, 2);
 
  return Starlab::min(radius, r_he_core);
}

real single_star::main_sequence_time(const real mass) {

  real t_ms = (2550. + 667.*pow(mass,2.5) + pow(mass,4.5) )
    / (0.0327*pow(mass,1.5) + 0.346*pow(mass,4.5));

  return t_ms;
}

real single_star::main_sequence_time() {

  real t_ms = (2550. + 667.*pow(relative_mass,2.5) 
	       + pow(relative_mass,4.5) )
    / (0.0327*pow(relative_mass,1.5) 
       + 0.346*pow(relative_mass,4.5));

  return t_ms;
}

real single_star::hertzsprung_gap_time(const real mass, const real t_ms) {

  return 0.54*t_ms/(mass*(mass - 2.1) + 23.3);
}

real single_star::hertzsprung_gap_time(const real t_ms) {
 
  return 0.54*t_ms/(relative_mass*(relative_mass - 2.1) + 23.3);
}

real single_star::helium_giant_time(const real mass, const real t_ms) {

  real l_bms = base_main_sequence_luminosity(mass);
  real l_he = helium_giant_luminosity(mass);

  return t_ms*l_bms/(l_he*(pow(mass, 0.42) + 0.8));
}

real single_star::helium_giant_time(const real t_ms) {

  real l_bms = base_main_sequence_luminosity();
  real l_he  = helium_giant_luminosity();

  return t_ms*l_bms/(l_he*(pow(relative_mass, 0.42) + 0.8));
}

real single_star::base_giant_branch_time(const real mass, const real t_ms) {

  real t_g    = 0.15 * t_ms;
  real l_g    = giant_luminosity(mass);
  // (SPZ+GN: 26 Jul 2000) mass added to function call:
  real l_bagb = base_agb_luminosity(l_g, mass);

  return t_g - t_g*pow((l_g/l_bagb), 0.855);
}

real single_star::base_giant_branch_time(const real t_ms) {

  real t_g    = 0.15 * t_ms;
  real l_g    = giant_luminosity();
  real l_bagb = base_agb_luminosity(l_g);

  return t_g - t_g*pow((l_g/l_bagb), 0.855);
}

real single_star::base_giant_time(const real mass, const real t_ms) {

  real t_gs  = 0.15*t_ms;
  real l_g   = giant_luminosity(mass);
  real l_agb = agb_luminosity(mass);

  return t_gs*pow(l_g/l_agb, 0.855);
}

real single_star::base_giant_time(const real t_ms) {

  real t_gs  = 0.15*t_ms;
  real l_g   = giant_luminosity();
  real l_agb = agb_luminosity();

  return t_gs*pow(l_g/l_agb, 0.855);
}

real single_star::nucleair_evolution_time(const real mass) {

  real l_g = giant_luminosity(mass);
  //(SPZ+GN: 26 Jul 2000) mass should be passed through
  real l_bagb = base_agb_luminosity(l_g, mass); 
  real l_agb = agb_luminosity(mass);

  real t_ms = main_sequence_time(mass);
  real t_g = 0.15*t_ms;
  real t_hg = hertzsprung_gap_time(mass, t_ms);
  real t_bgb = t_g - t_g*pow((l_g/l_bagb), 0.855);
  real t_he = helium_giant_time(mass, t_ms); //(SPZ+GN: 26 Jul 2000) add 'mass'
  real t_gs = t_g;
  real t_b = t_gs*pow(l_g/l_agb, 0.855);
  real t_agb = t_gs - t_bgb - t_b;

  return t_ms + t_hg + t_bgb + t_he + t_agb;
}


real single_star::nucleair_evolution_time() {

  real l_g = giant_luminosity();
  real l_bagb = base_agb_luminosity(l_g);
  real l_agb = agb_luminosity();

  real t_ms = main_sequence_time();
  real t_g = 0.15*t_ms;
  real t_hg = hertzsprung_gap_time(t_ms);
  real t_bgb = t_g - t_g*pow((l_g/l_bagb), 0.855);
  real t_he = helium_giant_time(t_ms); 
  real t_gs = t_g;
  real t_b = t_gs*pow(l_g/l_agb, 0.855);
  real t_agb = t_gs - t_bgb - t_b;

  return t_ms + t_hg + t_bgb + t_he + t_agb;
}

real single_star::base_main_sequence_luminosity() {

  real l_bms;
  if (relative_mass<=1.093) {
    l_bms =  (1.107*pow(relative_mass, 3.)
	      +   240.7*pow(relative_mass, 9.))
      /  (1. + 281.9*pow(relative_mass, 4.));
  }
  else {
    l_bms =  (13990.*pow(relative_mass, 5.))
      /             (pow(relative_mass, 4.)   
		     +  2151.*pow(relative_mass, 2)
		     +  3908.*relative_mass +9536.);
  }
  
  return l_bms;
}

real single_star::base_main_sequence_luminosity(const real mass) {

  real l_bms;
  if (mass<=1.093) {
    l_bms =  (1.107*pow(mass, 3.)+ 240.7*pow(mass, 9.))
      /  (1. + 281.9*pow(mass, 4.));
  }
  else {
    l_bms =  (13990.*pow(mass, 5.))
      / (pow(mass, 4.) + 2151.*pow(mass, 2)+3908.*mass +9536.);
  }
   
  return l_bms;
}

real single_star::giant_luminosity(const real mass) {

  real l_g = (2.15 + 0.22*pow(mass, 3.))*pow(mass, 2)
    / (5.0e-6*pow(mass, 4.)+1.4e-2*mass*mass + 1.);

  return l_g;
}

real single_star::giant_luminosity() {
  
  real l_g = (2.15 + 0.22*pow(relative_mass, 3.))
    *  pow(relative_mass, 2)
    / (5.0e-6*pow(relative_mass, 4.)
       +  1.4e-2*relative_mass*relative_mass + 1.);

  return l_g;
}

real single_star::helium_giant_luminosity(const real mass) {

  // Prescription according to EFT89.
  //real l_bms = base_main_sequence_luminosity(mass);
  //return 0.763*pow(mass, 0.46)*l_bms
  //  + 50.0/pow(mass, 0.1);

  // Changed (SPZ+GN:09/98) as in Tout et all. 97.
  real l_g = base_giant_branch_luminosity(mass);
  real l_HeI = base_agb_luminosity(l_g, mass);
  real l_He = 49*pow(mass, 0.364) + 0.86*pow(mass, 4.);

  return Starlab::min(l_He, l_HeI);
}

real single_star::helium_giant_luminosity() {

  //	Helium core bunring giant.
  //real l_bms = base_main_sequence_luminosity();
  //return 0.763*pow(relative_mass, 0.46)*l_bms
  //  + 50.0/pow(relative_mass, 0.1);

  // Changed (SPZ+GN:09/98) as in Tout et all. 97.
  real l_g = base_giant_branch_luminosity(relative_mass);
  real l_HeI = base_agb_luminosity(l_g, relative_mass);
  real l_He = 49*pow(relative_mass, 0.364)
    +  0.86*pow(relative_mass, 4.);

  return Starlab::min(l_He, l_HeI);

}

real single_star::base_giant_branch_luminosity(const real mass) {

  real kappa, lambda;
  real log_mass = log10(mass);

  if (mass<1.334) {
    kappa  = 0.2594 + 0.1348*log_mass;
    lambda = 0.144 - 0.8333*log_mass;
  }
  else {
    kappa  = 0.092088 + 0.059338*log_mass;
    lambda = -0.05713 + log_mass*(0.3756 -0.1744);
  }

  real l_bms = base_main_sequence_luminosity(mass);

  return l_bms*pow(10., kappa + lambda);
}


real single_star::base_giant_branch_luminosity() {

  real kappa, lambda;
  real log_mass = log10(relative_mass);

  if (relative_mass<1.334) {
    kappa  = 0.2594 + 0.1348*log_mass;
    lambda = 0.144 - 0.8333*log_mass;
  }
  else {
    kappa  = 0.092088 + 0.059338*log_mass;
    lambda = -0.05713 + log_mass*(0.3756 -0.1744);
  }

  real l_bms = base_main_sequence_luminosity();

  return l_bms*pow(10., kappa + lambda);
}

real single_star::base_agb_luminosity(const real l_g) {

  // EFT89 method.
  //return l_g + 2.e3;
  
  // Changed (SPZ+GN:09/98) as in Tout et all. 97.
  real l_agb=2454.7;
  if (relative_mass >
      cnsts.parameters(upper_ZAMS_mass_for_degenerate_core)) 
    l_agb = 6.5*pow(Starlab::min(relative_mass, 20.), 3.25);
  
  return Starlab::max(l_agb, l_g);
}

real single_star::base_agb_luminosity(const real l_g,
				      const real mass) {

  real l_agb=2454.7;
  if (mass >
      cnsts.parameters(upper_ZAMS_mass_for_degenerate_core)) 
    l_agb = 6.5*pow(Starlab::min(mass, 20.), 3.25);
  
  return Starlab::max(l_agb, l_g);
}

real single_star::agb_luminosity(const real mass) {

  return mass*(4.0e3 + 5.0e2*mass);
}

real single_star::agb_luminosity() {
 
  return relative_mass*(4.0e3 + 5.0e2*relative_mass);
}

real single_star::maximum_luminosity(const real mass) {

  return 4000*mass + 500*pow(mass,2);
}

real single_star::maximum_luminosity() {
  return 4000*relative_mass + 500*pow(relative_mass,2);
}
     
// Helium core lifetime.
// Tabulation method is better since it depends on the mass
// of the core only.
real single_star::helium_time() {

  real m = get_total_mass();

  // Optional but rather old
  // Iben and Tutukov, 1985, ApJSS, 58, 661.
  // t_he =  0.56*t_ms/pow(relative_mass, 0.52);
  
  // Helium Star lifetime from 
  // Pols O.R., 1993,
  // PhD Thesis, University of Amsterdam, P13, Eq. 2.15

  real t_he = 0;
  if (m<=0.7) 
    t_he = 10.76*pow(m, -3.75); 
  else if(m<=1.6) 
    t_he = 17.1*pow(m, -2.45); 
  else if(m<=4.8) 
    t_he = 11.48*pow(m, -1.6); 
  else    
    t_he = 2.37*pow(m, -0.6);

  return t_he;
}

real single_star::temperature() {
  //  cerr<<"single_star::temperature()"<<endl;
  //  PRC(luminosity);PRL(effective_radius);

  // Effective radius is the radius of the star as it really is.
  // radius is the equilibrium radius of the star.
  // return pow(1130.*luminosity/(radius*radius), 0.25); // in [1000K]


  real T_eff = cnsts.parameters(Tsun)
    * pow(luminosity/pow(effective_radius, 2), 0.25); // in [K]
  
  return T_eff;
}

real single_star::magnitude() {

  return -2.5*log10(luminosity) + 4.75 - bolometric_correction();
}

// real function bolometric_correction()
// calculates bolometric correction for giant stars
// input: stellar effective surface temperature in kKelvin.
// output: bolometric correction for giant star in cgs.
real single_star::bolometric_correction() {

  // temperature() is defined in Kelvin.
  // here we should use old 10^3K implementation 
  // (SPZ+GN: 1 Oct 1998)
  real temp_in_kK = 0.001 * temperature();
  real bc;

  if (temp_in_kK < 4.195)
    bc = 2.5*log10((1.724e-7*pow(temp_in_kK,11.) + 1.925e-2)
		   / (1. + 1.884e-9*pow(temp_in_kK,14.)));
  else if (temp_in_kK >=  4.195 &&
	   temp_in_kK <= 10.89)
    bc = 2.5*log10((7.56e-2*pow(temp_in_kK,1.5))
		   / (1. + 6.358e-5*pow(temp_in_kK, 4.5)));
  else
    bc = 2.5*log10((2728/pow(temp_in_kK,3.5) + 1.878e-2*temp_in_kK)
		   /(1. + 5.362e-5*pow(temp_in_kK,3.5)));

  return bc;
}

// wind_velocity(): Calculates stellar wind velocoty.
// Steller wind velocity is 2.5 times stellar escape velocity
real single_star::wind_velocity() {

  real v_esc2 = cnsts.physics(G)
    * cnsts.parameters(solar_mass)*get_total_mass()
    / (effective_radius*cnsts.parameters(solar_radius));
  real v_wind = 2.5*sqrt(v_esc2)/cnsts.physics(km_per_s);

  return v_wind;
}

real single_star::kelvin_helmholds_timescale() {

  // Equilibrium radius is 'single_star::radius'
  // Effective radius may be puffed up by accretion or subtraction of mass.

  return 31.56*pow(relative_mass,2.)/(radius*luminosity); // [Myr]
}

real single_star::nucleair_evolution_timescale() {
  // overloaded for main_sequence:: as
  // t_nuc = 10^10 [years] Msun/Lsun.
  // Assumed that 0.1 Msun is thermalized.

  // but for giants etc. we use the following, based on
  // dimensional grounds and empirical comparison with
  // Pols, 1994, A&A 290, 119
  // (SPZ+GN:30 Sep 1998)
  
  //  real fused_mass = 0.1*relative_mass;
  //
  //  real t_nuc = cnsts.parameters(energy_to_mass_in_internal_units)
  //             * fused_mass/luminosity;
  //
  //  real t_kh = kelvin_helmholds_timescale();
  //
  //  return sqrt(t_nuc * t_kh);
  

  // (GN+SPZ Apr 28 1999) giant lifetime is ~ 10% of ms life time

  return 0.1*main_sequence_time();

}

real single_star::dynamic_timescale() {

  return 5.08e-11*sqrt(pow(radius, 3.)/relative_mass);    // [Myr]
}


bool single_star::low_mass_star() {

  bool is_low_mass_star = false;

  if (!remnant()) {
    if(relative_mass <= cnsts.parameters(low_mass_star_mass_limit)) 
      is_low_mass_star = true;
  }
  else if (get_total_mass() <=
	   cnsts.parameters(low_mass_star_mass_limit)) {
    is_low_mass_star = true;
  }

  return is_low_mass_star;
}

bool single_star::medium_mass_star() {
  
  return (!low_mass_star() && !high_mass_star())?true:false;
}

bool single_star::high_mass_star() {
  
  if(remnant())
    return (get_total_mass()>cnsts.parameters(medium_mass_star_mass_limit))
      ?true :false;
  else
    return (get_relative_mass()>cnsts.parameters(medium_mass_star_mass_limit))
      ?true :false;
}

void single_star::read_element() {

  real total_mass, log_temp, log_lum;
  real bol_corr, temp, mv;
  int type = (int)star_type;
  cin >> type >> relative_age >> relative_mass
      >> total_mass >> core_mass >> radius >> log_temp
      >> log_lum >> bol_corr >> mv >> velocity
      >> magnetic_field >>  rotation_period;

  envelope_mass = total_mass - core_mass;
  temp = pow(log_temp, 10);
  luminosity = pow(log_lum, 10);
}

void single_star::put_element() {

  ofstream outfile("binev.data", ios::app|ios::out);
  if (!outfile) cerr << "error: couldn't create file binev.data"<<endl;

  real bol_corr = bolometric_correction();
  real mv = magnitude();
  outfile << star_type << " " << relative_age << " "
	  << relative_mass << " " << get_total_mass() << " "
	  << core_mass << " " << radius << " "
	  << log10(temperature()) << " "
	  << log10(luminosity) << " " << bol_corr << " "
	  << mv << " " << " " << velocity << " "
	  << magnetic_field <<  rotation_period << endl;

  outfile.close();
}

void single_star::dump(ostream & s, bool brief) {

  if (brief) {
    
    s << get_node()->format_label() << " "
      << get_current_time() << " "
      << "0 1   ";
    s << get_node()->format_label() << " "	  
      << get_element_type() << " "
      << get_total_mass() << " "
      << radius << "    ";
    s << "0 0 0 0" << endl;
  }
  else
    {
      real bol_corr = bolometric_correction();
      int s1, s2, s3, s4, s5, s6;
      s1 = get_spec_type(Emission);
      s2 = get_spec_type(Blue_Straggler);
      s3 = get_spec_type(Barium);
      s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
      s5 = get_spec_type(Runaway);
      s6 = get_spec_type(Merger);

      s << "1\n ";
      s << get_element_type() << " " << s1 << " " << s2 << " "
        << s3 << " " << s4 << " " << s5 << " " << s6 << endl;
      s <<" " << identity
	<< " " << current_time
	<< " " << relative_age 
	<< " " << relative_mass
	<< " " << get_total_mass()
	<< " " << core_mass
	<< " " << radius
	<< " " << effective_radius
	<< " " << log10(temperature())
	<< " " << log10(luminosity)
	<< " " << bol_corr
	<< " " << magnitude()
	<< " " << velocity
	<< " " << magnetic_field
	<< " " << rotation_period
	<< endl;
    }
}

void single_star::dump(const char * filename, bool brief) {

  ofstream s(filename, ios::app|ios::out);
  if (!s) cerr << "error: couldn't create file "<<filename<<endl;

  if (brief) {
    
    s << get_node()->format_label() << " "
      << get_current_time() << " "
      << "0 1   ";
    s << get_node()->format_label() << " "	  
      << get_element_type() << " "
      << get_total_mass() << " "
      << radius << " ";
    s << "0 0 0 0" << endl;

  }
  else
    {

      real bol_corr = bolometric_correction();
      int s1, s2, s3, s4, s5, s6;
      s1 = get_spec_type(Emission);
      s2 = get_spec_type(Blue_Straggler);
      s3 = get_spec_type(Barium);
      s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
      s5 = get_spec_type(Runaway);
      s6 = get_spec_type(Merger);

      s << "1\n ";
      s << get_element_type() << " " << s1 << " " << s2 << " "
        << s3 << " " << s4 << " " << s5 << " " << s6 << endl;
      s <<" " << identity
	<< " " << current_time
	<< " " << relative_age
	<< " " << relative_mass
	<< " " << get_total_mass()
	<< " " << core_mass
	<< " " << radius
	<< " " << effective_radius
	<< " " << log10(temperature())
	<< " " << log10(luminosity)
	<< " " << bol_corr
	<< " " << magnitude()
	<< " " << velocity
	<< " " << magnetic_field
	<< " " << rotation_period
	<< endl;

    }
  s.close();
}

void single_star::put_state() {

  star_state str;
  str.init_star_state(dynamic_cast(single_star*, this));
  str.make_star_state(dynamic_cast(single_star*, this));

  str.put_star_state();
}

void single_star::put_hrd(ostream & s) {

  int s1, s2, s3, s4, s5, s6;
  s1 = get_spec_type(Emission);
  s2 = get_spec_type(Blue_Straggler);
  s3 = get_spec_type(Barium);
  s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
  s5 = get_spec_type(Runaway);
  s6 = get_spec_type(Merger);

  s << "1\n ";
  s << s1 << " " << s2 << " " << s3 << " " << s4 << " "
    << s5 << " " << s6 << " ";
  s << get_total_mass()
    << " " << log10(temperature()) 
    << " " << log10(luminosity) 
    << endl;
}

void single_star::refresh_memory() {

  previous.star_type = get_element_type();
  previous.current_time = current_time;
  previous.last_update_age = last_update_age;
  previous.next_update_age = next_update_age;
  previous.relative_age = relative_age;
  previous.relative_mass = relative_mass;
  previous.envelope_mass = envelope_mass;
  previous.core_mass = core_mass;
  previous.COcore_mass = COcore_mass;

  previous.radius = radius;
  previous.effective_radius = effective_radius;
  
  //             Pulsars
  previous.magnetic_field  = magnetic_field;
  previous.rotation_period = rotation_period;
  previous.birth_mass      = birth_mass;
}

void single_star::recall_memory() {

  star_type        = previous.star_type;
  current_time     = previous.current_time;
  last_update_age  = previous.last_update_age;
  next_update_age  = previous.next_update_age;
  relative_age     = previous.relative_age;
  relative_mass    = previous.relative_mass;
  envelope_mass    = previous.envelope_mass;
  core_mass        = previous.core_mass;
  COcore_mass      = previous.COcore_mass;

  radius           = previous.radius;
  effective_radius = previous.effective_radius;

  // Pulsars
  magnetic_field   = previous.magnetic_field;
  rotation_period  = previous.rotation_period;
  birth_mass       = previous.birth_mass;
}

// Mass transfer timescale is checked and updated at the moment the
// mass-ratio is reversed.

// in version 3.3, the mass transfer timescale is updated each
// timestep in ::recursive(), therefore this function is not required.
// It hangs the code
// because mass loss by stellar wind occurs after
// ::mass_ratio_mdot_limit().
// therefore the mass ratio > 1 after ::recursive step.
// (SPZ+GN:11 Oct 1998)
real single_star::mass_ratio_mdot_limit(real mdot) {

    
  // No limit for the moment.
  return mdot;
    
  real accretor_mass = 0;

  if (is_binary_component()) 
    get_companion()->get_total_mass();

  if (accretor_mass<get_total_mass()) {
    real mdot_max = get_total_mass() - accretor_mass;
    if (mdot>mdot_max) 
      mdot = mdot_max;
  }

  int p = cerr.precision(HIGH_PRECISION);
  PRC(accretor_mass);PRL(get_total_mass());
  cerr.precision(p);
  
  return mdot;
}

// Computes expansion of acceptor due to mass accretion.
// At moment accretor fills own Roche-lobe mass transfer becomes
// inconservative.
real single_star::accretion_limit(const real mdot, const real dt) {

  // Conservative mass transfer.
  // return mdot;

  // Non-conservative mass transfer.
  real r_rl = get_binary()->roche_radius(this);
  real mdot_kh = dt*relative_mass/kelvin_helmholds_timescale();
  real accretion = log10(r_rl/effective_radius)
    / pow(10, expansionB(relative_mass));

  accretion = Starlab::max(accretion, 0.);
  real mdot_max = mdot_kh*pow(accretion, 1./expansionA(relative_mass));
  mdot_max = Starlab::max(mdot_max, 0.);	

  return Starlab::min(mdot, mdot_max);

  // (SPZ+GN: 26 Jul 2000) Test Kelvin Helmholtz accretion
  // (GN Jul 28 1999) test non conservative Algol evolution
  // return min(mdot, mdot_kh);

}

void single_star::adjust_donor_radius(const real delta_m) {


  real zeta_th = zeta_thermal();

  // -1 because delta_m is defined positive and star loses mass.
  real delta_r = -1 * effective_radius * zeta_th * delta_m/relative_mass;

  effective_radius += delta_r;
  
  // Effective radius return to equilibrium radius when mass transfer
  // is on nuclear timescale
  // (SPZ+GN:30 Sep 1998)
  if (is_binary_component() &&
      get_binary()->get_current_mass_transfer_type()==Nuclear) {
    effective_radius = radius;
  }


}


// Adding mass to a star causes it to expand.
void single_star::adjust_accretor_radius(const real mdot, const real dt) {

  // function returns directly: effective radius is used to determine RLOF.
  // Allowing accretor to grow in size results in contact systems.
  // do not allow bloating
  // (SPZ+GN:28 Sep 1998)
  //(GN+SPZ Apr 28 1999) star do bloat however... bloating on again
  //  return;
  
  //cerr<<"void star::adjust_accretor_radius()"<<endl;
  //cerr<<"pre radius: "<<radius<<" "<<effective_radius<<endl;

  if (mdot>0) {
    real mdot_kh = relative_mass*dt/kelvin_helmholds_timescale();

    real r_fr = expansionA(relative_mass)*log10(mdot/mdot_kh)
      + expansionB(relative_mass);
    if (r_fr>0.5) // (GN+SPZ Apr 28 1999) was 1 
      r_fr = 0.5;
    // (GN+SPZ Apr 28 1999) radius is equilibrium radius
    //    effective_radaius = radius = radius*pow(10, pow(10, r_fr));

    real r_l = get_binary()->roche_radius(this);
    effective_radius = Starlab::max(Starlab::min(r_l, effective_radius), 
				    radius*pow(10, pow(10, r_fr)));
             
  }
}

real single_star::expansionA(const real m) {

  // Lineair interpolation.
  real rc, inc;
  if (m<=3) {
    rc = 0.120;  	//base on: 	m1 = 2; A1 = 0.599;
    inc = 0.359; 	//	  	m2 = 3; A2 = 0.719;
  }
  else if (m<=5) {
    rc = 0.144;	        //based on:	m1 = 3; A1 = 0.719;
    inc = 0.289;	//		m2 = 5; A2 = 1.006;
  }
  else if (m<=10) {
    rc = 0.123;	        //based on:	m1 = 5; A1 = 1.006;
    inc = 0.393;	//		m2 = 10; A2 = 1.619;
  }
  else {
    rc = 0.085;	        //based on:	m1 = 10; A1 = 1.619;
    inc = 0.772;	//		m2 = 17; A2 = 2.212;
  }
     
  return inc + rc*m;
}

real single_star::expansionB(const real m) {

  //              Lineair interpolation.
  real rc, inc;
  if (m<=3) {
    rc = -0.273;     //base on:      m1 = 2; B1 = -1.374;
    inc = -0.828;    //              m2 = 3; B2 = -1.647;
  }
  else if (m<=5) {
    rc = -0.192;     //based on:     m1 = 3; B1 = -1.647;
    inc = -1.071;    //              m2 = 5; B2 = -2.030;
  }
  else if (m<=10) {
    rc = -0.106;    //based on:     m1 = 5; B1 = -2.030;
    inc = -1.50;    //              m2 = 10; B2 = -2.560;
  }
  else {
    rc = -7.08e-2;   //based on:     m1 = 10; B1 = -2.560;
    inc = -1.852;    //              m2 = 17; B2 = -3.056;
  }

  real value = inc + rc*m;

  return Starlab::min(1., value);
}

// Merges cores and envelopes of two stars.
// Star is updated.
star* single_star::merge_elements(star* str) {

  merge_two_stars_story(str->get_element_type());

  star* merged_star = this;

  real m_conserved = get_total_mass() + str->get_total_mass();

  if (str->get_element_type()!=Main_Sequence) {

    // adding two cores of giants together should not result in
    // rejuvenation.
    // previous method appeared to make mergers too old.
    // (SPZ+GN:11 Oct 1998)
      
    add_mass_to_core(str->get_core_mass());

    if ((str->get_element_type()==Neutron_Star ||
	 str->get_element_type()==Black_Hole)   &&
	core_mass < cnsts.parameters(helium2neutron_star)) {
      real dm = cnsts.parameters(helium2neutron_star) - core_mass;
      core_mass = cnsts.parameters(helium2neutron_star);
      if (envelope_mass<dm)
	envelope_mass = 0;
      else
	envelope_mass -= dm;
    }

    //		What to do here is put in SPH!
    if (str->get_envelope_mass()>0) 
      add_mass_to_accretor(0.5*str->get_envelope_mass());

  }
  else {

    add_mass_to_accretor(str->get_total_mass());
  }

  if (get_total_mass()-m_conserved > 1.e-11) {
    cerr.precision(HIGH_PRECISION);
    cerr << "ERROR: Mass is not conserved in single_star::merge elements()"
	 << endl;

    PRC(get_total_mass());PRC(m_conserved);
    PRL(get_total_mass()-m_conserved);
    
    exit(-1);
  }

  spec_type[Merger]=Merger;
  adjust_next_update_age();

  
  // Redundant: in core mass addition star is no rejuvenated but aged.
  // addition of envelope mass causes rejuvenation.
  // (SPZ+GN:27 Sep 1998)
  //if (str->get_element_type()!=Main_Sequence) {
  //real aged = 0.01*(core_mass/agb_core_mass(relative_mass));
  //relative_age = (1+aged)*next_update_age;
  //}


  if (get_relative_mass() >= cnsts.parameters(massive_star_mass_limit) &&
      hydrogen_envelope_star()) 
    merged_star =  reduce_mass(envelope_mass);
  else {

    instantaneous_element();
    // calling evolve element may cause segmentation fault if star changes.
    //evolve_element(current_time);
  }

  cerr << "Merged star: "<<endl;
  merged_star->dump(cerr, false);

  //++++++Identify primary after merger
  //  if (is_binary_component())
  //    get_companion()->set_element_type(NAS);

  return merged_star;
}

// Determine mass transfer stability according to
// `zeta' discription.
real single_star::mass_transfer_timescale(mass_transfer_type &type) {

  type = Unknown;
  
  real z_l=0, mass_ratio = 1;
  if (is_binary_component()) {
    mass_ratio = get_companion()->get_total_mass()/get_total_mass();
    z_l = get_binary()->zeta(this, get_companion());
  }

  real z_ad = zeta_adiabatic();
  real z_th = zeta_thermal();

  real mtt;
  if (z_ad>z_l && z_th>z_l) {
    //         mtt = 11.*kelvin_helmholds_timescale();
    mtt = nucleair_evolution_timescale();
    type = Nuclear;
  }
  else if (z_ad>z_l && z_l>z_th) {
    mtt = kelvin_helmholds_timescale();

    type = Thermal;
  }
  else if (z_l>z_ad) {
    mtt = sqrt(kelvin_helmholds_timescale()*dynamic_timescale());
    //if (medium_mass_star())
    //  mtt = 0.09*kelvin_helmholds_timescale();
    //else
    //  mtt = sqrt(kelvin_helmholds_timescale()*dynamic_timescale());
    type = Dynamic;
  }
  else {
    cerr << "No clear indication for mass transfer timescale: "
	 << "Kelvin-Helmholds time-scale assumed."<<endl;

    mtt = kelvin_helmholds_timescale();
    type = Unknown;
  }

  if (low_mass_star()) {

    real mdot = get_binary()
      ->mdot_according_to_roche_radius_change(this,
					      get_companion());
    if (mdot>0) {
      real mtt_rl = get_relative_mass()/mdot;
      if(mtt>mtt_rl) {
	mtt = mtt_rl;
	type = AML_driven;
      }
      
    }
  }

  if (REPORT_MASS_TRANSFER_TIMESCALE) {
    cerr << "single_star::mass_transfer_timescale()" << endl;
    cerr << "    star id = " << identity
	 << "  Zeta (lobe, ad, th) = ("
	 << z_l <<", "<<z_ad<<", "<<z_th<<") : " << endl;
    cerr << type_string(type);
    cerr<<":    dm/dt=" <<get_relative_mass()/(mtt*1.e+6)
	<< " [Msun/yr]" << endl;
  }

  return mtt;
}

//		Stellar stability functions.
real single_star::zeta_adiabatic() {

  // (GN+SPZ Apr 28 1999) this is used by sub_giant: all stars with HW87
  // all stars should have their own actually...(when we have time)

  // (GN+SPZ Apr 28 1999) fit from Lev Yungelson private communication
  // for giants with not too convective envelope = radiative envelope

  real r_dconv = 2.4*pow(relative_mass,1.56);
  if (relative_mass > 10 )
    r_dconv = 5.24*pow(relative_mass,1.32);
  else if (relative_mass > 5)
    r_dconv = 1.33*pow(relative_mass,1.93);
    
  if (radius < r_dconv) {
    return 12.25;
    //    cerr << "radius < r_dconv" << endl;
  }
  else {
    //   		Hjellming and Webbink 1987 ApJ, 318, 804
    real x = core_mass/get_total_mass();
    real A = -0.220823;
    real B = -2.84699;
    real C = 32.0344;
    real D = -75.6863;
    real E = 57.8109;

    // (GN+SPZ Apr 28 1999) not for (sub) giants    
    //  if (low_mass_star())
    //    z = -cnsts.mathematics(one_third);
    //  else 
    return A + x*(B + x*(C + x*(D + x*E)));

  }

}

// Values of zeta are changed (SPZ+GN:28 Sep 1998)
real single_star::zeta_thermal() {
  //cerr<<"real single_star::zeta_thermal()"<<endl;

  real z;
  if (low_mass_star())
    z = 0;
  else 
    z = 0; // (GN+SPZ Apr 28 1999) radius determined by core only (was -1) 

  return z;
}

// add two cores. Is performed in ::merge_elements();
void single_star::add_mass_to_core(const real mdot) {

  if (mdot<=0) {
    cerr << "single_star::add_mass_to_core(mdot=" << mdot << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
    
  }
  else {
    adjust_accretor_age(mdot, false);
    core_mass += mdot;
    accreted_mass += mdot;
    
    if (relative_mass<get_total_mass()) 
      update_relative_mass(get_total_mass());
    
  }

}


//		general mass transfer utilities.
// Increase donor mass and possibly relative_mass of donor.
// Check mass-transfer timescales before use.
real single_star::add_mass_to_accretor(const real mdot) {

  //  cerr << "add_mass_to_accretor"<<endl;
  //  PRL(mdot);
  //  dump(cerr, false);
  
  if (mdot<=0) {
    cerr << "single_star::add_mass_to_accretor(mdot=" << mdot << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;

    set_spec_type(Accreting, false);
    
    return 0;
  }
  else {
    adjust_accretor_age(mdot);
    envelope_mass += mdot;
    accreted_mass += mdot;
    
    if (relative_mass<get_total_mass()) 
      update_relative_mass(get_total_mass());
    
    set_spec_type(Accreting);
    
  }

  
  return mdot;
}

real single_star::add_mass_to_accretor(real mdot, const real dt) {

  //  cerr << "add_mass_to_accretor"<<endl;
  //  PRL(mdot);
  //  dump(cerr, false);

  if (mdot<0) {
    cerr << "single_star::add_mass_to_accretor(mdot=" << mdot 
	 << ", dt=" << dt << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
    return 0;
  }

  mdot = accretion_limit(mdot, dt);

  adjust_accretor_age(mdot);
  envelope_mass += mdot;
  accreted_mass += mdot;

  if (relative_mass<get_total_mass()) 
    update_relative_mass(get_total_mass());

  adjust_accretor_radius(mdot, dt);
  
  set_spec_type(Accreting);

  return mdot;
}

// Post-evolution stellar update.
// Makes sure age and radius are updated.
void single_star::update() {

  // New core mass determination occurs in ::evolve_element.
  // (SPZ+GN:09/1998)
  // real m_tot = get_total_mass();
  // core_mass = helium_core_mass();
  // envelope_mass = m_tot - core_mass;

  core_radius = helium_core_radius();

  // (GN+SPZ Apr 28 1999)
  // effective_radius can be larger than  radius
  effective_radius = Starlab::max(radius,effective_radius);


  // last update age is set after stellar expansion timescale is set.
  // (GN+SPZ May  4 1999) last_update_age now used as time of last type change
  //  last_update_age = relative_age;

  detect_spectral_features();

}


void single_star::detect_spectral_features() {

  //		Clean old spectral specifications.
  set_spec_type(Emission, false);
  set_spec_type(Blue_Straggler, false);
  set_spec_type(Barium, false);

  // set new
  if (is_binary_component())
    if (!remnant() && 
	!get_companion()->hydrogen_envelope_star() &&
	accreted_mass >= cnsts.parameters(Barium_star_mass_limit)
	* relative_mass)
      set_spec_type(Barium);
}

// In the case of a binary, the companion star might accrete from a
// stellar wind or post-AGB mass-shell.
// Bondi, H., and Hoyle, F., 1944, MNRAS 104, 273 (wind accretion.
// Livio, M., Warner, B., 1984, The Observatory 104, 152.
real single_star::accrete_from_stellar_wind(const real mdot, const real dt) {

  //  PRC(mdot);PRL(dt);

  real alpha_wind = 0.5;
  real v_wind = get_companion()->wind_velocity();

  real acc_radius = pow(cnsts.physics(G)*cnsts.parameters(solar_mass)
			* get_total_mass()
			/ pow(v_wind*cnsts.physics(km_per_s), 2),2)
    / cnsts.parameters(solar_radius);
  real wind_acc = alpha_wind/(sqrt(1-pow(get_binary()->get_eccentricity(), 2))
			      * pow(cnsts.parameters(solar_radius)
				    * get_binary()->get_semi(),2));
  real v_factor = 1/pow(1+pow(velocity/v_wind, 2), 3./2.);

  real mass_fraction = acc_radius*wind_acc*v_factor;
  // (GN+SPZ May  4 1999) Do not multiply with dt*cnsts.physics(Myear)!

  //  PRC(v_wind);PRC(acc_radius);PRC(wind_acc);PRL(v_factor);
  //  PRL(mass_fraction);PRL(mdot);
  
  mass_fraction = Starlab::min(0.9, mass_fraction);

  return add_mass_to_accretor(mass_fraction*mdot, dt);
}

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real single_star::tf2_energy_diss(const real eta) {

  const real coeff2[] = {-0.397, 1.678, 1.277, -12.42, 9.446, -5.550};

  real y = log10(eta);
  real logT = ((((coeff2[5]*y + coeff2[4])*y + coeff2[3])*y 
		+ coeff2[2])*y + coeff2[1])*y + coeff2[0];

  return pow(10., logT);
}

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real single_star::tf3_energy_diss(const real eta) {

  const real coeff3[] = {-0.909, 1.574, 12.37, -57.40, 80.10, -46.43};

  real y = log10(eta);
  real logT = ((((coeff3[5]*y + coeff3[4])*y + coeff3[3])*y  
		+ coeff3[2])*y + coeff3[1])*y + coeff3[0];

  return pow(10., logT);
}

//	pretty-print
void single_star::print_status() {
  
  star_state str;
  str.init_star_state((single_star*)this);
  str.make_star_state((single_star*)this);
  
  cout << " [M = " << get_total_mass() 
       << ", R = " << effective_radius 
       << ", v = " << velocity << "] "
       << type_string(get_element_type()) << " ";
  str.put_star_state();

}

//	print data for EPJ Roche program.
void single_star::print_roche() {

  real r = effective_radius;
  if (effective_radius<radius) r = 1.e6;

  cerr << get_total_mass() << " " << r << " ";

}

void  single_star::set_spec_type(star_type_spec s,
				 bool on) { 	// default =true

  if (on) spec_type[s] = s;
  else    spec_type[s] = NAC;
}


// Angular momentum of homogeneous sphere.
real single_star::angular_momentum() {
       
  real m = get_total_mass()*cnsts.parameters(solar_mass);
  real r = effective_radius*cnsts.parameters(solar_radius);

  // (GN+SPZ May  5 1999) effective_radius may increase when rl shrinks
  if (is_binary_component()) {
    r = Starlab::min(r, 
		     get_binary()->roche_radius(this)*cnsts.parameters(solar_radius));
  }

  real o = 0;                            // default rotation.
  if(rotation_period>0)
    o = 2*PI/rotation_period;
  else if (is_binary_component()) 
    o = 2*PI/(get_binary()->get_period()*cnsts.physics(days));

  return gyration_radius_sq()*m*o*pow(r, 2);
	
}

real single_star::rejuvenation_fraction(const real mdot_fr) {

  real rejuvenation = (1-pow(mdot_fr,
			     cnsts.parameters(rejuvenation_exponent)));

  // no rejuvenation if companion has no hydrogen envelope.
  if(is_binary_component() &&
     !get_companion()->hydrogen_envelope_star()) 
    rejuvenation = 1;

  return rejuvenation;
}

void single_star::stellar_wind(const real dt) {
  //  cerr << "void single_star::stellar_wind(const real dt=" << dt << ")" << endl;

  //    PRL(last_update_age);
  //    PRL(next_update_age);
  //    PRL(relative_age);
  // (GN+SPZ Apr 28 1999) wind for low mass stars per phase
  real end_time = next_update_age - last_update_age;
  //    real prev_rel_time = max(0.,previous.relative_age - last_update_age);
  //    real relative_time = min(relative_age - last_update_age, end_time);
  real relative_time = relative_age - last_update_age;

  //    PRL(end_time);
  //    PRL(prev_rel_time);
  //    PRL(relative_time);
  //    PRL(wind_constant);

  // for high mass stars over whole evolution
  // (GN May 12 1999)
  // except for stars more massive than 85 that become WR after ms immediately
  // (GN Apr 19 2004) temporary fix for Hyper_Giants with M > 85 through merger
  if (get_element_type() == (int)Hyper_Giant ||
      (relative_mass >= cnsts.parameters(super_giant2neutron_star) &&
      relative_mass < 85.)) {
    end_time = nucleair_evolution_time();
    //      prev_rel_time = previous.relative_age;
    relative_time = relative_age;
  }

  //    PRC(end_time);
  //    PRC(dt);
  //    PRC(prev_rel_time);
  //    PRC(relative_time);
  //    PRL(wind_constant);

  real wind_mass = wind_constant 
    * (pow(relative_time/end_time,
	   cnsts.parameters(massive_star_mass_loss_law))
       -  pow((relative_time-dt)/end_time,
	      cnsts.parameters(massive_star_mass_loss_law)));

  //  PRC(wind_mass);
  // Previous second term according to GN.
  //	           -  pow((prev_rel_time)/end_time,
  //			cnsts.parameters(massive_star_mass_loss_law)));
  // (GN+SPZ Apr 28 1999) wind induced helium star formation can happen
  // because stellar_wind is last function in evolve_element
  if (wind_mass>=envelope_mass) {
    wind_mass = envelope_mass;
    radius = core_radius;
  }
  
  //  PRL(wind_mass);
  //  dump(cerr, false);
  if (is_binary_component())
    get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
  else {
    // (GN Oct 11 1999) for single stars: previous used for stellar wind! (?)
    // refresh_memory();

    reduce_mass(wind_mass);
  }

  return;
}

void single_star::update_relative_mass(const real new_relative_mass) {

  relative_mass = new_relative_mass;
  adjust_next_update_age();
  update_wind_constant();

}

void single_star::lose_envelope_decent() {

  //  cerr << "single_star::lose_envelope_decent" << endl;
  //  dump(cerr, false);
  if (envelope_mass>0) {
    if (is_binary_component()) {
      get_binary()->adjust_binary_after_wind_loss(
						  this, envelope_mass, POST_AGB_TIME);
    } 
    else 
      reduce_mass(envelope_mass);
  }

  // (GN Apr 19 2004) added check for Merged & Disrupted, otherwise merged black hole radius 
  // is set to 0, causing inf temperature, nan bol_cor and magnitude
  if (is_binary_component()  && 
      get_binary()->get_bin_type() != Merged &&
      get_binary()->get_bin_type() != Disrupted) {
    get_companion()->set_effective_radius(get_companion()->get_radius());
    get_companion()->set_spec_type(Accreting, false);
  }
}

// wind_constant is fraction of envelope lost in lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 3 Oct 1998)
void single_star::update_wind_constant() {
  
#if 0
  if (relative_mass >= cnsts.parameters(massive_star_mass_limit))
    wind_constant = (get_relative_mass()-core_mass)
      * cnsts.parameters(massive_star_envelope_fraction_lost);
  else 
    wind_constant = get_relative_mass()
      * cnsts.parameters(non_massive_star_envelope_fraction_lost);
    
#endif
  // (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense

  //  cerr << "update_wind_constant"<<endl;

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_constant = meader_fit_dm;
    else {// constant
      real final_mass = 43; // final mass after ms
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // no wind for low mass ms stars
    wind_constant = 0;
  }

  wind_constant = Starlab::max(wind_constant, 0.0);

  //  cerr << "single_star" << endl;
  //  PRL(wind_constant);

}


real single_star::potential_energy() {
  
  real GM2_R = cnsts.physics(G)*pow(cnsts.parameters(solar_mass), 2)
    / cnsts.parameters(solar_radius);
  real p = GM2_R*get_total_mass()*get_total_mass()
    / get_effective_radius();
     
  return -p;
}

real single_star::kinetic_energy() {
  
  real Mkm_s2 = cnsts.parameters(solar_mass)
    * pow(cnsts.physics(km_per_s), 2);
  real k = 0.5*Mkm_s2*get_total_mass()*pow(velocity, 2);
     
  return k;
}

real single_star::total_energy() {
  return kinetic_energy() + potential_energy();
}

real single_star::get_evolve_timestep() {

  // (GN+SPZ Apr 28 1999) was a bit too small
  //  return max(next_update_age - relative_age
  //	     -0.5*cnsts.safety(minimum_timestep),
  //	     cnsts.safety(minimum_timestep));

  // (GN+SPZ May  5 1999) type change time must be small because of rapid
  // growth of giants at end phase 0.0001 seems to be OK (?)
  //  return max(next_update_age - relative_age - (0.5*0.001), 0.001);

  real time_step = Starlab::max(next_update_age - relative_age, 0.0001);
  // (SPZ Sept 5 2004) old description gave problems with dynamics.
  //  time_step = Starlab::min(time_step, 0.05);
  // (SPZ Jan 1 2006) changed again from min() to max() to satisfy McScatter
  //  time_step = Starlab::max(time_step, 1.0);
  if (get_use_hdyn()) {
    time_step = Starlab::min(time_step, 0.05);
  }
  return time_step;

}

// Groenewegen & de Jong 1993, A&A 267,410 
real single_star::final_core_mass() {

  if (maximum_luminosity() < 15725) {                     // Mc < 0.73
    return sqrt(maximum_luminosity()/47488. + 0.18) + 0.015;
  }
  else {
    // (SPZ+GN: 26 Jul 2000) was: pow(relative_mass,0.19), but makes
    //                       white dwarfs too massive.
    //                       see Nelemans, Y.PZ.V. 2000 A&A (submitted)
    return  0.46 + maximum_luminosity()/(46818*pow(relative_mass,0.25));
  }

}

// WARNING: this function spoils the memory of a single star
real single_star::get_dlogR_dT() {
  
  real dt = 1.e-4*get_next_update_age();
  real r = get_effective_radius();
  refresh_memory();
  evolve_element(dt);
  real rdr = get_effective_radius();
  recall_memory();

  real dlogr_dt = (r-rdr)/dt;

  return dlogr_dt; // in Rsun per Mur
}
