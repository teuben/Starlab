#include "constants.h"
#include "stdfunc.h"

// constructor for constant class.
// Currently there are not functionalities.

//stellar_evolution_constants::stellar_evolution_constants() {
//  //Dummy function in order to allow compilation.
//}


real stellar_evolution_constants::mathematics(mathematical_constant pm) {

    switch(pm) {
	case one_third:                        return 0.33333333333333333333;
             break;
	case two_third:                        return 0.66666666666666666666;
             break;
	case pi:                               return 3.14159265358979323846;
	  break;                               // of the orde of 1
	case two_pi:                           return 6.28318530717958647692;
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(mathematical_constant "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::physics(physics_constants pp) {

    // physics parameters 
    switch(pp) {
	case gravitational_constant:
	case G:                             return 6.67e-8;    // [cgs]
             break;
        case speed_of_light:					    
        case C:                             return 2.9979e+10; // [cm/s]
             break;
        case million_years:
        case Myear:                         return 3.15e+13;   // [s]
             break;
        case seconds_per_day:
        case days:                          return 8.6400e+4;   // [s]
             break;
        case kilometer_per_second:
        case km_per_s:                      return 1.0e+5;      // [cm/s]
             break;
        case kilometer_in_centimeters:
        case kilometer:                     return 1.0e+5;   // [cm]
             break;
	case nucleair_efficiency:                      return 0.007;
	      break;                          // nuc. energy production eff
	                                      // Delta E = 0.007 Mc^2
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(physics_constants "
	          << pp << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::super_nova_kick(
				  super_nova_kick_distribution pk,
				  const real v_disp){
  
  // Super Nova Kick functions.
  // kick velicity imparted to the remnant during the supernova event.
  // changing the kick velocity requires recompilation.

  //	Paczynsky's velocity distribution with a dispersion
  //	of 600 km/s is prefered by Hartman.

  //	Gauss kick velocity distribution.
  //	Lyne, A.G., \& Lorimer, D.R., 1994, Nature 369, 127
  //	propose v_k = 450 +- 90 km/s
  //	this is comparable to a one dimensional gaussion
  //	of 260 km/s.

  //    Hansen & Phinney (1997, MNRAS 291, 569) kick velocity distribution.
  //    They prever a Maxwellian with a velocity dispersion of
  //    270 km/s.
  
  // selected kick distribution imparted to a newly formed neutron star
  // in a supernova. 
    switch(pk) {
        case no_velocity_kick:              return 0;
             break;                                
	case Maxwellian_velocity_kick:      return
					    random_maxwellian_velocity(v_disp);
             break;
        case internally_decided_velocity_kick:
	case Paczynski_velocity_kick:       return
					    random_paczynski_velocity(v_disp);
             break;
	case delta_function_velocity_kick:  return v_disp;
             break;
    default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "super_nova_kick(super_nova_kick_distribution "
	          << pk << ")"
		  << endl;
             exit(1);
    }
}  

real stellar_evolution_constants::parameters(astronomical_scale_parameter pa) {

    // Astronomical distance scale parameters in centimeters [cm]
    switch(pa) {
	case PC:
	case parsec:                             return 3.0857e+18;
             break;                                
	case AU:                      
	case astronomical_unit:                  return 1.496e+13; 
             break;                                
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(astronomical_scale_parameter "
	          << pa << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(solar_parameter ps) {

    // Solar parameters in cgs.
    switch(ps) {
	case solar_mass:
        case Msun:                      return 1.989e+33; //[gram]
             break;                                
	case solar_radius:
	case Rsun:                      return 6.96e+10;//[centimeter]
             break;                                
	case solar_luminosity:
	case Lsun:                      return 3.862e+33; // [erg/s]
             break;                                
	case solar_temperature:
	case Tsun:                      return 5770;     // [Kelvin]
             break;
        case energy_to_mass_in_internal_units:         return 
				   cnsts.physics(nucleair_efficiency)
				 * pow(cnsts.physics(C), 2)
                                 * cnsts.parameters(solar_mass)
                                 / (cnsts.parameters(solar_luminosity) *
				    cnsts.physics(Myear));
	     break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(solar_parameter "
	          << ps << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(pulsar_initial_conditions pp) {

    // Netron star basic parameters.
    switch(pp) {
	case pulsar_magnetic_field:              return 12;     // [log Gs]
             break;                                
	case pulsar_pulse_period:                return 0.1;    // [s]
             break;                                
        case kanonical_neutron_star_radius:      return 1.5e-5; // [cm]
             break;                               
        case kanonical_neutron_star_mass:        return 1.34;   // [msun]
             break;                                
	case maximum_neutron_star_mass:          return 2.0;    // [Msun]
             break;                              // was 1.5
	case minimum_neutron_star_mass:          return 0.0925; // [Msun]
             break;                              // Newtonian polytrope 
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(pulsar_initial_conditions "
	          << pp << ")"
		  << endl;
             exit(1);
    }
}

real stellar_evolution_constants::parameters(stellar_mass_limits pm) {

    // Switch statement for certain limits. All values are in
    // solar units [Msun].
    switch(pm) {
	case low_mass_star_mass_limit:           return 1.5;
             break;                                
	case medium_mass_star_mass_limit:        return 15;
             break;                               
	case massive_star_mass_limit:            return 25;
             break;                               
	case upper_ZAMS_mass_for_degenerate_core: return 2.3;
             break;                                
	case minimum_main_sequence:              return 0.075;
             break;                              // for pop III: 0.095 [Msun] 
	case maximum_planet_mass:                return 0.00314;
             break;                              // Lynden-Bell&O'Dwyer2001
	case helium_dwarf_mass_limit:            return 0.45;
             break;                                
	case carbon_dwarf_mass_limit:            return 1.2;
             break;
	case Chandrasekar_mass:                  return 1.44;
             break;                                
	case helium2neutron_star:                return 2.2;
             break;                                
        case COcore2black_hole:                  return 5; // was 10
                                                      // (SPZ+GN: 27 Jul 2000)
             break;                                
	case super_giant2neutron_star:           return 8;
             break;                                
        case super_giant2black_hole:             return 25; // was 40
             break;                                
	case maximum_main_sequence:              return 100;
             break;                               
	case minimum_helium_star:                return 0.33;
             break;
	case helium_star_lifetime_fraction:      return 0.9;
             break;
	case helium_star_final_core_fraction:    return 0.80;
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(stellar_mass_limits "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}


bool stellar_evolution_constants::parameters(boolean_parameter pb) {

    switch(pb) {                                
	case hyper_critical:                        return false;
             break;                                // 10^8 x Eddington
                                                   // see: Chevalier 1993
        case super_giant_disintegration:               return false;
             break;
        case proto_star_to_binary:                     return false;
             break;						    
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(boolean_parameter "
	          << pb << ")"
		  << endl;
             exit(1);
    }
}
    
real stellar_evolution_constants::parameters(accretion_parameter pa) {

    switch(pa) {                        // Mass of accretion to core.
	case black_hole_accretion_limit:        return 0.1;
             break;                             // 
	case neutron_star_accretion_limit:      return 0.05;
             break;                             //
	case white_dwarf_accretion_limit:       return 0.01;
             break;                             //
        case thermo_nuclear_flash:              return 1;  // used to be 0.05;
             break;                             // mass fractioncauses flash
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(accretion_parameter "
	          << pa << ")"
		  << endl;
             exit(1);

      }
}


real stellar_evolution_constants::parameters(model_parameter pm) {

    switch(pm) {
	case star_formation_efficiency:                return 1.0;
             break;                             // [%] 
	case star_formation_timescale:                 return 1.0; // [Myr]
             break;                                 // if(<0) abs times KH timescale 
	case magnetic_mass_limit:                      return 0.7;
             break;                             // [msun] magnetic if(M<this)
	case magnetic_braking_exponent:                return 2.5;
	      break;                        
	case corotation_eccentricity:                  return 0.001;
	      break;                          // corotation if(ecc<this)
	case tidal_circularization_radius:             return 5.0;
	      break;                          
	case core_overshoot:                           return 0.125;
	      break;                          // overshoot fraction.
        case hydrogen_fraction:                        return 0.7;
              break;                          // X
 	case common_envelope_efficiency:               return 4;
	      break;                          // (alpha_ce)
        case envelope_binding_energy:                  return 0.5;
	      break;                          // (lambda)
	case specific_angular_momentum_loss:           return 3.;
	      break;                          // (beta)
        case dynamic_mass_transfer_gamma:               return 1.75;
              break;			      // (gamma)
        case non_massive_star_envelope_fraction_lost:  return 0.03;
	      break;                          
        case massive_star_envelope_fraction_lost:      return 0.9;
	      break;               
        case relaxation_driven_mass_loss_constant:     return 1;           
             break;						       
        case massive_star_mass_loss_law:               return 6.8;
	      break;                          
        case time_independent_mass_loss_law:           return 1;
	      break;
        case Darwin_Riemann_instability_factor:        return
	                                       cnsts.mathematics(one_third);
	      break;
        case homogeneous_sphere_gyration_radius_sq:    return 0.4;
	      break;
        case radiative_star_gyration_radius_sq:        return 0.03;
	      break;
        case convective_star_gyration_radius_sq:       return 0.2;
	      break;
        case rejuvenation_exponent:                    return 1;
	      break;
        case spiral_in_time:                           return 0.0005; // Myr
	      break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(model_parameter "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}


real stellar_evolution_constants::parameters(observational_parameter pm) {

    switch(pm) {
        case B_emission_star_mass_limit:               return 0.1;
             break; 
        case Barium_star_mass_limit:                   return 0.01;
	      break;                        
        case Blue_straggler_mass_limit:                return 0;
	      break;                        
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(observational_parameter "
	          << pm << ")"
		  << endl;
             exit(1);
    }
}

// Safetly parameters means that changing these may cause
// unexpected errors and changes the bahavious of the program.
real stellar_evolution_constants::safety(safety_parameter ps) {

  switch(ps) {
    case timestep_factor:                      return 0.01;       
          break;                             // 0.01 ok, but big for AM CVns
    case maximum_binary_update_time_fraction: return 0.9;
          break;                            // see also star_to_dyn
    case minimum_timestep:                   return 1.e-11; // Do not change!
          break;                                   // == 10*HubbleT*precision
    case minimum_mass_step:                  return 1.e-5;
          break;                                   // Tricky used to be 1.e-7
    case maximum_timestep:                   return 1000; // was 1, was 1000  
          break;                                   // 2.5 works fine but slow
    case maximum_recursive_calls:            return 100;
          break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "parameters(safety_parameter "
	          << ps << ")"
		  << endl;
             exit(1);
  } 
}

real stellar_evolution_constants::star_to_dyn(dynamics_update_parameter dup) {

    switch(dup) {
        case stellar_mass_update_limit:                return 0.001;
             break; 
        case semi_major_axis_update_limit:             return 0.001;
	      break;                        
        case binary_update_time_fraction:              return 0.001; //was 0.05
             break;
        default:
	     cerr << "\nNo recognized option in "
		     "stellar_evolution_constants::"
		     "star_to_dyn(dynamica_update_parameter "
	          << dup << ")"
		  << endl;
             exit(1);
    }
}






