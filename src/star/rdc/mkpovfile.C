
//// rdc_movie: read input snapshot and reduce it to simple readable
////            ASCII for movie program.
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "single_star.h"
#include "sstar_to_dyn.h"
#include "dyn.h"

#ifndef TOOLBOX

#define NNEBULAE 10

typedef struct nebulae {
  int id;
  real time;
  vector pos;
};
nebulae nebula[NNEBULAE];

#define NSN_REMNANT 4
typedef struct sn_remnants {
  int id;
  bool bh;
  real time;
  vector pos;
};
sn_remnants sn_remnant[NSN_REMNANT];

#define NCOLLISIONS 12
typedef struct collisions {
  int starid;
  real time;
};
collisions collision[NCOLLISIONS];

enum filter_type {Visual=0, Radio, X_rays}; 

local filter_type get_filter_type(char fltr) {

  switch(fltr) {
    case 'V': return Visual;
              break;
    case 'R': return Radio;
              break;
    case 'X': return X_rays;
              break;
    default: cerr << "No such filter"<<endl;
             exit(-1);
  };

}

void new_camera_position(vector &v_new,
			 vector first_pos,
			 vector last_pos,
			 int nsteps,
			 int nsnap) {

  real fstep = (1.0*nsnap)/(1.0*nsteps);
  v_new[0] = first_pos[0] + fstep * (last_pos[0] - first_pos[0]);
  v_new[1] = first_pos[1] + fstep * (last_pos[1] - first_pos[1]);
  v_new[2] = first_pos[2] + fstep * (last_pos[2] - first_pos[2]);
  
}

#if 0
void new_camera_position(vector &v_new,
			 vector v_old,
			 real theta_rotation,
			 real phi_rotation,
			 int nsteps,
			 int nsnap) {

  
  real theta = nsnap*theta_rotation/nsteps;
  real phi   = nsnap*phi_rotation/nsteps;

  PRC(nsnap);PRC(theta);PRL(phi);
  real r = abs(v_old);
  v_new[0] = r * sin(theta);
  v_new[1] = r * sin(theta) * cos(phi);
  v_new[2] = r * cos(theta);

  cerr << "vector = " << v_new << endl;

}
#endif

local void print_camera_on_star(dyn* b, vector com_pos,
				int camera_on_star_id,
                                real r_star,
				real aperture, int blur_samples) {

    vector cam_view = b->get_pos() + b->get_vel();
    vector cam_pos = b->get_pos() + r_star*b->get_vel();

    // Look at the clusters CoM for now.
    cam_view[0] = 0;
    cam_view[1] = 0;
    cam_view[2] = 0;
    
    vector normal;
    normal[0] = sqrt(pow(abs(cam_pos), 2)
              -      pow(cam_pos[1], 2)
              -      pow(cam_pos[2], 2));
    normal[1] = 1;
    normal[2] = (-cam_pos[0] * normal[0] - cam_pos[1])/cam_pos[2];

    cout << "// Normal to camera " << endl;
    cout << "   #declare normal_to_camera = < " << normal[0] << ", "
	                                        << normal[1] << ", "
	                                        << normal[2] << " >" << endl;
    
    cout << "// camera located on star #" << camera_on_star_id << endl;
    cout << "camera {" << endl
	 << "   location < " << cam_pos[0] << ", "
                             << cam_pos[1] << ", "
                             << cam_pos[2] << " >" << endl
	 << "   look_at  < " << cam_view[0] << ", "
	                     << cam_view[1] << ", "
	                     << cam_view[2] << " >" << endl
	 << "   blur_samples " << blur_samples << endl; 
    if (aperture>0)
	cout << "   focal_point < " << cam_view[0] << ", "
	                            << cam_view[1] << ", "
	                            << cam_view[2] << " >" << endl
	     << "   aperture " << aperture << endl;
    cout << "}" << endl << endl;

}


local bool print_camera_position_recursive(dyn* b, vector cam_pos,
					   int camera_on_star_id,
                                           real r_star,
					   real aperture, int blur_samples) {

  if (b->get_oldest_daughter()) {
    
    vector com_pos  = b->get_pos() - cam_pos;


    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() >= 2) 
	return print_camera_position_recursive(bb, com_pos,
                                               camera_on_star_id, r_star,
                                               aperture, blur_samples);

      else if (bb->get_index()==camera_on_star_id) {
	print_camera_on_star(bb, com_pos, 
                             camera_on_star_id, r_star,
                             aperture, blur_samples);

	return true;
      }

  }
  return false;
}

local void print_filename_counter(int counter, ostream& s) {

    if (counter>=10000) { 
      cout << "\nToo namy filenames in print_povray_header()" << endl;
      exit(1);
    }
    else if (counter<10)
      s << "000" << counter;
    else if (counter<100)
      s << "00" << counter;
    else if (counter<1000)
      s << "0" << counter;
    else
      s << counter;
}

local void print_pl_nebula(vector pos, real scale) {

    cout << "object { Pl_Nebula scale " << scale
	 << " translate < " << pos[0] << ", "
                            << pos[1] << ", "
		            << pos[2] << " > }"
	 << endl;
}

local void print_sn_nebula(vector pos, real scale) {

  cout << "object { SN_Remnant scale " << scale
       << " translate < " << pos[0] << ", "
                          << pos[1] << ", "
		          << pos[2] << " > }"
       << endl;
}

local void print_some_data(dyn *b) {

  if (b->get_oldest_daughter()) {

    real time = b->get_starbase()->conv_t_dyn_to_star(b->get_system_time());
    
    int nts=0, nbs=0, nss=0;
    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() > 2)
	nts++;
      else if (bb->n_leaves() >= 2)
	nbs++;
      else 
	nss++;
  }
}

local void print_hertzsprung_Russell_diagram(dyn* b, vector cam_pos) {

  if (b->get_oldest_daughter()) {

    print_some_data(b);

    int stype_s[no_of_star_type_summ];
    for (int i=0; i<no_of_star_type_summ; i++) 
      stype_s[i]=0;
    
    cout << "#declare Stellar_HRD = union {" << endl;
    cout << "   object { Y_Axis } " << endl;
    cout << "   object { X_Axis }\n " << endl;
    
    for_all_leaves(dyn, b, bi) {
      
      star_type_spec tpe_class = NAC;
      spectral_class star_class;
      stellar_type stype = NAS;
      stellar_type_summary sstype = ZAMS;
      real t_cur, t_rel, m_rel, m_env, m_core, co_core; 
      real T_eff, L_eff, p_rot, b_fld;
	//real T_eff=0, L_eff=0;
      if (bi->get_use_sstar()) {
	stype = bi->get_starbase()->get_element_type();
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	
	T_eff = bi->get_starbase()->temperature();
	L_eff = bi->get_starbase()->get_luminosity();
      }
      else if (bi->get_star_story()) {
	extract_story_chapter(stype, t_cur, t_rel, m_rel, m_env, 
			      m_core, co_core,
			      T_eff, L_eff, p_rot, b_fld,
			      *bi->get_star_story());
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	star_class = get_spectral_class(T_eff);
#if 0
	if (find_qmatch(bi->get_star_story(), "Type")) {
	cerr <<"Reading"<<endl;
	stype = extract_stellar_type_string(
                getsq(bi->get_star_story(), "Type"));
	sstype = summarize_stellar_type(stype);
	stype_s[dynamic_cast(int, sstype)]++;
	T_eff = getrq(bi->get_star_story(), "T_eff");
	star_class = get_spectral_class(T_eff);
	L_eff = getrq(bi->get_star_story(), "L_eff");
#endif	
      }
      else {
	cout << "    No stellar information found for: ";
	bi->pretty_print_node(cout);
	return;
      }
      
      real xt_pos = (4.78 - log10(T_eff))/4.78;
      real yl_pos = (log10(L_eff) + 1)/7.;

      if (xt_pos>0&&xt_pos<1 && yl_pos>0&&yl_pos<1)
	cout << "   object { Red_Sphere translate < "
	     << xt_pos << ", "
	     << yl_pos << ", "
	     << "0 > }" << endl;
    }

    cout << "\n} // End Stellar_HRD" << endl;
    cout << "object { Stellar_HRD translate < "
	 << cam_pos[0]+0.5 << ", "
	 << cam_pos[1]-1 << ", "
	 << cam_pos[2] + 2
	 << " > }\n" << endl;

  }
}


local void print_star(dyn *bi, vector pos,
		      real scale_L, filter_type filter) {

  // To solar radii
  //  vector pos = bi->get_pos() - dc_pos;
  pos[0] = bi->get_starbase()->conv_r_dyn_to_star(pos[0]);
  pos[1] = bi->get_starbase()->conv_r_dyn_to_star(pos[1]);
  pos[2] = bi->get_starbase()->conv_r_dyn_to_star(pos[2]);

  real time = bi->get_starbase()->conv_t_dyn_to_star(bi->get_system_time());

  // And now to parsec
  real Rsun_per_parsec = cnsts.parameters(solar_radius)
                       / cnsts.parameters(parsec);
  pos[0] *= Rsun_per_parsec;
  pos[1] *= Rsun_per_parsec;
  pos[2] *= Rsun_per_parsec;
  

     star_type_spec tpe_class = NAC;
     spectral_class star_class;
     stellar_type stype = NAS;
     stellar_type_summary sstype = ZAMS;
     real t_cur, m_rel, m_env, m_core, co_core, T_eff, L_eff, p_rot, b_fld;
     real t_rel=0, R_eff=0;
     real M_tot, U, B, V, R, I;	
     if (bi->get_use_sstar()) {
       	stype = bi->get_starbase()->get_element_type();
	M_tot  = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
        t_cur = bi->get_starbase()->get_current_time();
        t_rel = bi->get_starbase()->get_relative_age();
        T_eff = bi->get_starbase()->temperature();
        L_eff = bi->get_starbase()->get_luminosity();
        star_class = get_spectral_class(T_eff);
	R_eff = bi->get_starbase()->get_effective_radius();
	ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     U, B, V, R, I);

     }
     else if (bi->get_star_story()) {

       extract_story_chapter(stype, t_cur, t_rel, m_rel, m_env, 
			     m_core, co_core,
			     T_eff, L_eff, p_rot, b_fld,
			     *bi->get_star_story());
       T_eff *= 1000;

       M_tot = m_env + m_core;
       sstype = summarize_stellar_type(stype);
       star_class = get_spectral_class(T_eff);
       
       ltm_to_ubvri(log10(L_eff), log10(T_eff), M_tot,
		     U, B, V, R, I);
       
       if (find_qmatch(bi->get_star_story(), "Class"))
	 tpe_class = extract_stellar_spec_summary_string(
             getsq(bi->get_star_story(), "Class"));
       if (L_eff>0)
          R_eff = pow(T_eff/cnsts.parameters(solar_temperature), 2)
	       / sqrt(L_eff);
     }
     else {
       cout << "    No stellar information found for: ";
       bi->pretty_print_node(cout);
       return;
     }

     int known_at;
     bool is_known = false;
     bool should_be_known = false;
     
     if (remnant(stype) && t_rel <= 1) 
       should_be_known = true;
     
     if (remnant(stype)) {
       for (int i=0; i<NNEBULAE; i++)
	 if (bi->get_index() == nebula[i].id) {
	   is_known = true;
	   known_at = i;
	 }

     }

     if (should_be_known) {
       if (!is_known)
	 for (int i=0; i<NNEBULAE; i++) 
	   if (nebula[i].id < 0) {
	     nebula[i].id  = bi->get_index();
	     nebula[i].pos = pos;
	     nebula[i].time = t_rel;
	     known_at = i;
	   }
       
       if (stype== Helium_Dwarf || Carbon_Dwarf || Oxygen_Dwarf) 
	 print_pl_nebula(nebula[known_at].pos, scale_L * (0.5 + t_rel));
       else {
#if 0
	 print_sn_nebula(nebula[known_at].pos, scale_L * (0.5 + t_rel));
	 if (t_rel<0.01) 
	   cout << "object { single_star \n"
		<< "         finish { ambient <"
		<< 0.4 << ", " << 0.6 << ", " << 0.7 << ">}\n" 
		<< "         scale " << scale_L * pow(100 * (0.01-t_rel), 3)
		<< " translate < "
		<< pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
		<< endl;
#endif
	 
       }
     }
     else if(is_known)
       nebula[known_at].id = -1;

     //  Add SN remnant
     for (int i=0; i<NSN_REMNANT; i++)
       if(bi->get_index() == sn_remnant[i].id &&
	  t_cur >= sn_remnant[i].time && t_cur < sn_remnant[i].time+0.025) {
	 if (sn_remnant[i].bh)
	   print_sn_nebula(pos, scale_L * (0.5 + t_rel));
	 else
	   print_pl_nebula(pos, scale_L * (0.5 + t_rel));
	 if(t_cur >= sn_remnant[i].time && t_cur < sn_remnant[i].time+0.01) {
	   cout << "object { single_star \n"
		<< "         finish { ambient <";
	   if (sn_remnant[i].bh)
	     cout << 0.4 << ", " << 0.6 << ", " << 0.7 << ">}\n";
	   else
	     cout << 0.56 << ", " << 0.14 << ", " << 0.14 << ">}\n";
	   cout << "         scale " << scale_L * pow(100 * (0.01-t_rel), 3)
		<< " translate < "
		<< pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
		<< endl;
	 }
       }

     //  Add collisions
     for (int i=0; i<NCOLLISIONS; i++)
       if(bi->get_index() == collision[i].starid &&
	  t_cur >= collision[i].time && t_cur < collision[i].time+0.01) {
	 cout << "object { OStar scale "
	      << scale_L * pow(100 * (0.01 - (t_cur-collision[i].time)), 3)
	      << " translate < "
	      << pos[0] << ", " << pos[1] << ", " << pos[2]-0.5 << " > }"
	      << endl;
       }
     
     // The eye works as an exponential detector.
     real logL_eff = scale_L * sqrt(log10(1 + L_eff));
#if 0       
     cout << "object { "
	  << type_short_string(star_class) << "Star scale "
	  << logL_eff << " translate < "
	  << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	  << endl;

     if (tpe_class == Accreting)
       cout << "object { Accretion_Disc scale "
	    << logL_eff << " translate < "
	    << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	    << endl;
#endif

     real G;
     real apparent_size = -1;
     if (filter == Visual) {

	// Transform UBVRI into RGB on a rather clumsy way.
	real UM100 = -7.38974,  UM1 = 5.17096,  UM01 = 17.9608;
	real BM100 = -6.23550,  BM1 = 5.06279,  BM01 = 16.7472;
	real VM100 = -5.91088,  VM1 = 4.46967,  VM01 = 15.0682;
	real RM100 = -15.9099,  RM1 = 4.13746,  RM01 = 13.7278;
	real IM100 = -25.9089,  IM1 = 3.81839,  IM01 = 12.0133;

	RM100 = -16; RM01 = 14;
	R = max(0., min(1., (R-RM100)/(RM1-RM100)));

	G = 1 - min(1., max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));

	BM100 = -7; BM01 = 17;
	B = max(0., min(1., (BM1-B)/(BM01-BM100)));

	apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
     }
     else {

       real VM100 = -5.91088,  VM1 = 4.46967,  VM01 = 15.0682;
       apparent_size = 0.1 *scale_L * M_tot;

       switch(stype) {
          case Carbon_Star:  
          case Helium_Star:  
          case Helium_Giant: R = 0; B=0; // green
	    G = 1 - min(1., max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Carbon_Dwarf: 
          case Helium_Dwarf: 
          case Oxygen_Dwarf: G=0; B=0; // Red
	    R = 1 - min(1., max(0., sqrt(abs((V-VM1)/(VM01-VM100)))));
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Thorn_Zytkow: R = 1; G=0; B=0;
	    apparent_size = 0.1 * scale_L * (VM01-V)/(VM01-VM100);
	                      break;
          case Xray_Pulsar:  
          case Radio_Pulsar: 
          case Neutron_Star: R = 1; G=0; B=0; // Blue
	                       B = -0.3 * log10(min(1., p_rot));
			       apparent_size = 0.1 * scale_L * B;
 	                       break;
          case Black_Hole:   R = 1; G=1; B=1;
	                      break;
          default:
	       R = 0; G=0; B=0;
	    apparent_size = -1;
	                      break;
       
       };
     }
     
     cout << "object { single_star \n"
          << "         finish { ambient <"
          << R << ", " << G << ", " << B << ">}\n" 
	  << "         scale " << apparent_size << " translate < "
	  << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	  << endl;

     if (tpe_class == Accreting)
       cout << "object { Accretion_Disc scale "
	    << logL_eff << " translate < "
	    << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
	    << endl;
}

local void print_node(dyn *bi, vector pos, real mass_scale,
		      real mmax)  {

  real time = bi->get_system_time();

  real apparent_size = mass_scale * sqrt(bi->get_mass());
  real mmin = 0;
  real m = bi->get_mass();
  real R = sqrt((mmax-m)/(mmax-mmin));
  real G = 0.5;
  real B = sqrt(m/(mmax-mmin));
  if(m<0.001) {
    R = 0.5;
    G = 0.5;
    B = 0.5;
  }
  else {
    R = 1;
    G = 0;
    B = 0;
  }
  cout << "object { simple_star \n"
       << "         finish { ambient <"
       << R << ", " << G << ", " << B << ">}\n" 
       << "         scale " << apparent_size << " translate < "
       << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
       << endl;

//  cout << "object { single_star \n"
//       << "         finish { ambient <"
//       << R << ", " << G << ", " << B << ">}\n" 
//       << "         scale " << apparent_size << " translate < "
//       << pos[0] << ", " << pos[1] << ", " << pos[2] << " > }"
//       << endl;
}

local int print_povray_binary_recursive(dyn *b,
					vector dc_pos, 
					real mass_limit, real number_limit,
					bool povray, real scale_L,
					filter_type filter) {

  int nb = 0;
  if (b->get_oldest_daughter()) {
    
    vector r_com  = b->get_pos() - dc_pos;

    real m_tot = b->get_starbase()->conv_m_dyn_to_star(b->get_mass());

    for_all_daughters(dyn, b, bb)
      if (bb->n_leaves() >= 2) 
	nb += print_povray_binary_recursive(bb, dc_pos,
					    mass_limit, number_limit,
					    povray, scale_L, filter);
      else {
	print_star(bb, bb->get_pos()-r_com, scale_L, filter);
      }

  }
  return nb;
}

local void print_povtime(real time,
			 vector pos,
			 real scale=1,
			 real depth=0.25) {

    int p = cout.precision(LOW_PRECISION);
    cout << "text { ttf \"timrom.ttf\" ";
    if(time>0)
      cout << "\"Time = " << time << " Myr \" ";
    else
      cout << "\"Time = " << -time << " N-body \" ";

    cout << depth << ", 0\n"
	 << "       pigment { Red }\n"
         << "       translate < " << pos[0] << ", "
	                  << pos[1] << ", "
	                  << pos[2] << " >\n"
         << "       scale " << scale
	 << "}\n" << endl;
    cout.precision(p);
}

local void print_start_text(vector cam_pos) {
}

local void print_povray_stars(dyn *b, real mass_limit,
			      real number_limit,
			      bool povray,
			      real scale_L,
			      filter_type filter) {

  for (int i=0; i<NNEBULAE; i++) 
    nebula[i].id  = -1;
    
  bool cod = false;

  vector dc_pos = 0;
  bool try_com = false;
  if(abs(dc_pos) == 0) {
    if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {
      
      if (getrq(b->get_dyn_story(), "density_center_time")
	  != b->get_system_time()) {
	warning("mkpovfile: neglecting out-of-date density center");
	try_com = true;
      } else
	cod = true;
      
      dc_pos = getvq(b->get_dyn_story(), "density_center_pos");
    }

    if (try_com && find_qmatch(b->get_dyn_story(), "com_pos")) {

      if (getrq(b->get_dyn_story(), "com_time")
	  != b->get_system_time()) {
	warning("lagrad: neglecting out-of-date center of mass");
      } else
	dc_pos = getvq(b->get_dyn_story(), "com_pos");
    }
  }


  // For now: put denxity center in geometric origin
  dc_pos = 0;
  
  int ns=0, nb=0;
  for_all_daughters(dyn, b, bi) 
    if (bi->is_leaf()                    &&
	((bi->get_index()>0 && bi->get_index() <= number_limit) ||
	 bi->get_starbase()->get_total_mass() >= mass_limit)) {
	
      print_star(bi, bi->get_pos() - dc_pos, scale_L, filter);
      ns++;
    }

  for_all_daughters(dyn, b, bi) {
    if (!bi->is_leaf())                   
      nb += print_povray_binary_recursive(bi, dc_pos,
					  mass_limit, number_limit,
					  povray, scale_L, filter);
  }
  
if (!povray && ns+nb==0) 
    cout << "            (none)\n";
}

local void print_povray_bodies(dyn *b, real mass_limit,
			      real number_limit, real mmax, 
			      real scale_M) {


/*
  real mmin = VERY_LARGE_NUMBER, mmax = -1;
  for_all_leaves(dyn, b, bi) {
    real m = bi->get_mass();
    mmin = min(mmin, m);
    mmax = max(mmax, m);
  } 
*/
    
  int ns=0;
  for_all_leaves(dyn, b, bi) 
    if(bi->get_mass() >= mass_limit) {
	
      print_node(bi, bi->get_pos(), scale_M, mmax);
      ns++;
    }
  
  if (ns==0) 
    cout << "            (none)\n";
}

void print_povray_header(dyn* b, vector cam_pos,
                         int camera_on_star_id, 
			 real aperture, real gamma, 
			 int blur_samples,
                         int counter, real scale_L,
			 int horizontal, int vertical,
			 bool print_hrd) {

  cout << "\n\n// End of POVRay scenary." << endl;
  cout << "\n\n// Start new POVRay scenary." << endl;

  cout << "//       filename = starlab_";
  print_filename_counter(counter, cout);
  cout << ".pov\n";
    
  cout << "#include \"colors.inc\"" << endl
       << "#include \"shapes.inc\"" << endl
       << "#include \"textures.inc\"" << endl
       << "#include \"glass.inc\"" << endl;
  cout << "#include \"astro.inc\"" << endl;

  cerr<< endl;

  cout << "#version 3.0" << endl << endl;

  cout << "global_settings { " << endl
       << "  assumed_gamma " << gamma << endl 
       << "  max_trace_level 15" << endl
       << "  ambient_light White" << endl
       << "}"
       << endl << endl;

  if (camera_on_star_id<=0) {

    vector normal;
    normal[0] = sqrt(pow(abs(cam_pos), 2)
              -      pow(cam_pos[1], 2)
              -      pow(cam_pos[2], 2));
    normal[1] = 1;
    normal[2] = (- cam_pos[0] * normal[0] - cam_pos[1])/cam_pos[2];

    cout << "// Normal to camera " << endl;
    cout << "   #declare normal_to_camera = < " << normal[0] << ", "
	                                        << normal[1] << ", "
	                                        << normal[2] << " >" << endl;
    
    cout << "camera {" << endl
	 << "   location < " << cam_pos[0] << ", "
                             << cam_pos[1] << ", "
                             << cam_pos[2] << " >" << endl
	 << "   look_at  <0.0, 0.0,  0.0>" << endl
	 << "   blur_samples " << blur_samples << endl; 
    if (aperture>0)
      cout << "   focal_point <0, 0,  0>" << endl
	   << "   aperture " << aperture << endl;
    cout << "}" << endl << endl;

    //
    cout << "light_source { < " << cam_pos[0] << ", "
                                << cam_pos[1] << ", "
                                << cam_pos[2] << " > White }" << endl;

//    cout << "#object { Initial_text }\n" << endl;
    //
    
    real time = b->get_starbase()->conv_t_dyn_to_star(b->get_system_time());
    vector time_pos;
    if(vertical<=400) {         //320x240
      time_pos[0] = -10;
      time_pos[1] = 6;
    }
    else {        //640x480
      time_pos[0] = -16;
      time_pos[1] = -11;
    }
//    time_pos[0] = -4 - vertical * 0.0188;
//    time_pos[1] = -4 - horizontal * 0.0188;
//    time_pos[2] = 2.5*(5+cam_pos[2]);
    time_pos[0] = -16;
    time_pos[1] = -12;
    time_pos[2] = 1;
//
    real scale = 0.2;
    real depth = 0.25;
    print_povtime(time, time_pos, scale, depth);    
      
    if (print_hrd)
      print_hertzsprung_Russell_diagram(b, cam_pos);
    
  }
  else {
    cam_pos = 0;
    if(!print_camera_position_recursive(b, cam_pos, camera_on_star_id, scale_L,
					aperture, blur_samples)) {
      cout << "No star id = " << camera_on_star_id << " found." << endl;
      cout << " in print_camera_position_recursive() " << endl;
      exit(-1);
    }
  }

  cout << "// Here follow the definitions of the stars..." << endl
       << endl;
}

void print_mpegplayer_param_file(ostream& s,
				 int first_frame = 1,
				 int last_frame = 100,
				 int horizontal = 120,
				 int vertical = 90,
				 int GOP_size = 10) {
				 

     s <<"# mpeg_encode parameter file\n"
	  << "PATTERN         IBBBPBBBBP\n"
	  << "OUTPUT          starlab.mpg\n"
	  << endl;

     s << "YUV_SIZE        "<<vertical<<"x"<<horizontal<<"\n"
	  << "BASE_FILE_FORMAT        PPM\n"
	  << "INPUT_CONVERT   *\n"
	  << "GOP_SIZE        "<<GOP_size<<"\n" 
	  << "SLICES_PER_FRAME  1\n"
	  << endl;

     s << "INPUT_DIR       .\n"
	  << "INPUT\n"
	  << "starlab_*.ppm   [";
          print_filename_counter(first_frame, s);
      s   << "-";
          print_filename_counter(last_frame-1, s);
      s   << "]\n"
	  << "END_INPUT\n" << endl;

    s << "FORCE_ENCODE_LAST_FRAME\n"
//	 << "PIXEL           HALF\n"
	 << "PIXEL           WHOLE\n"
	 << "RANGE           10\n"
	 << endl;

    s << "PSEARCH_ALG     LOGARITHMIC\n"
	 << "BSEARCH_ALG     CROSS2\n"
	 << endl;

    s << "IQSCALE         8\n"
	 << "PQSCALE         10\n"
	 << "BQSCALE         25\n"
	 << endl;

    s << "REFERENCE_FRAME ORIGINAL" << endl;
}

void rdc_and_wrt_movie(dyn *b, bool povray, real scale_L, real mmax, 
		       char fltr) {

  // Fill in the collision database.
    collision[0].starid = 2; 
    collision[0].time = 50.8389;
#if 0
    collision[1].starid = 2;
    collision[1].time = 3.09252;
    collision[2].starid = 3;
    collision[2].time = 3.093;
  //  collision[0].starid = 1;
  //  collision[0].time = 1.65;
    //  collision[1].starid = 1;
    //  collision[1].time = 1.68;
    //  collision[2].starid = 1;
    //  collision[2].time = 1.73;
  collision[3].starid = 1;
  collision[3].time = 1.74;
  collision[4].starid = 1;
  collision[4].time = 2.2233267;
  collision[5].starid = 1;
  collision[5].time = 2.750325;
  collision[6].starid = 1;
  collision[6].time = 2.89577;
  collision[7].starid = 1;
  collision[7].time = 3.357066;
  collision[8].starid = 1;
  collision[8].time = 3.5157569;
  collision[9].starid = 1;
  collision[9].time = 3.84768766;
  collision[10].starid = 1;
  collision[10].time = 4.16801414;
  collision[11].starid = 1;
  collision[11].time = 4.6702693;

  sn_remnant[0].bh = true;
  sn_remnant[0].id = 4441;
  sn_remnant[0].time = 4.42;
  sn_remnant[1].bh = true;
  sn_remnant[1].id = 1;
  sn_remnant[1].time = 5.35;
  sn_remnant[1].bh = false;
  sn_remnant[2].id = 2056;
  sn_remnant[2].time = 6.14;
  sn_remnant[1].bh = false;
  sn_remnant[3].id = 3193;
  sn_remnant[3].time = 6.26;
#endif

    filter_type filter = get_filter_type(fltr);
    
    real mass_limit = 0;
    real number_limit = VERY_LARGE_NUMBER;
    if(b->get_use_sstar())
      print_povray_stars(b, mass_limit, number_limit, povray, scale_L, filter);
    else
      print_povray_bodies(b, mass_limit, number_limit, mmax, scale_L);
    
}

#else

//-----------------------------------------------------------------------------
//  main  --  driver to use  compute_mass_radii() as a tool
//-----------------------------------------------------------------------------
enum SF_base_type {No_base=0, StTr_station, SW_ISD}; 

main(int argc, char ** argv)
{
    bool  c_flag = false;      // if TRUE: a comment given on command line
    bool  v_flag = false;      // if TRUE: a comment given on command line
    bool  B_flag = false;  
    bool  print_hrd = false;

    char filter = 'V';
    real aperture    = -1;
    vector first_campos, last_campos, cam_pos;
    first_campos[2] = last_campos[2] = -10;
    real theta_rotate = 0;
    real phi_rotate   = 0;
    int blur_samples = 10;
    int camera_on_star_id = -1;
    int nstart = 0;

    char * mpeg_paramfile = "paramfile.mpeg";
    int horizontal =  90;
    int vertical   = 120;
    int first_frame = 1;
    int last_frame = 0;
    real scale_L   = 1;

    real gamma = 2.2;

    SF_base_type base_type = No_base;
    int nsteps = 32767;
    int nskip = 0; 
    int nint_skip = 0;
    
    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "a:B:b:F:f:g:H:hI:L:N:n:W:X:x:Y:y:Z:z:P:T:S:s:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'a': aperture = atof(poptarg);
	              break;
	    case 'B': base_type = (SF_base_type)atoi(poptarg);
	              break;
	    case 'b': blur_samples = atoi(poptarg);
	              break;
	    case 'F': filter = *poptarg;
	              break;
	    case 'f': mpeg_paramfile = poptarg;
	              break;
	    case 'g': gamma = atof(poptarg);
	              break;
	    case 'H': horizontal = atoi(poptarg);
	              break;
	    case 'h': print_hrd = true;
	              break;
	    case 'I': camera_on_star_id = atoi(poptarg);
	              break;
	    case 'W': vertical = atoi(poptarg);
	              break;
	    case 'L': scale_L *= atof(poptarg);
	              break;
	    case 'P': phi_rotate = atof(poptarg);
	              phi_rotate *= cnsts.mathematics(pi)/180;
	              break;
	    case 'T': theta_rotate = atof(poptarg);
	              theta_rotate *= cnsts.mathematics(pi)/180;
	              break;
	    case 'N': nsteps = atoi(poptarg);
	              break;
	    case 'n': nstart = atoi(poptarg);
	              break;
	    case 'X': first_campos[0] = atof(poptarg);
	              break;
	    case 'Y': first_campos[1] = atof(poptarg);
	              break;
	    case 'Z': first_campos[2] = atof(poptarg);
	              break;
	    case 'x': last_campos[0] = atof(poptarg);
	              break;
	    case 'y': last_campos[1] = atof(poptarg);
	              break;
	    case 'z': last_campos[2] = atof(poptarg);
	              break;
	    case 'S': nskip = atoi(poptarg);
	              break;
	    case 's': nint_skip = atoi(poptarg);
	              break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cout, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    dyn *b;

    int GOP_size = min(1000, nsteps);
    remove(mpeg_paramfile);
    ofstream os(mpeg_paramfile, ios::app|ios::out);
    if (!os) cerr << "\nerror: couldn't create file "
		  << mpeg_paramfile << endl;

    print_mpegplayer_param_file(os, first_frame, GOP_size,
				horizontal, vertical, GOP_size);
    os.close();

  int id;
  real time, mmax = -1;
  vector pos;
    
    int nread = nstart;
    int j=0;
    for (int i = 0; i < nskip; i++){
	cerr << " skipping snapshot " << i << endl;
	if (! forget_node(cin)) exit(1);
	if(j==nint_skip) {
	  nread++;
	  j=0;
	}
	else
	  j++;
    }
    
    int nsnap = 0;	
    while (b = get_dyn(cin)) {

      if(nsnap==0)
	for_all_leaves(dyn, b, bi) {
	  mmax = max(mmax, bi->get_mass());
	}

      //      addstar(b);                            
      //b->set_use_sstar(true); 
      // Do not make new stars but make sure the starbase is searched.
      //b->set_use_sstar(false);

	nsnap ++;
	nread ++;
	cerr << " reading snapshot " << nread << ", and "
	     << nsnap << " done." << endl;
	
	if (camera_on_star_id<=0) 
	  new_camera_position(cam_pos, first_campos, last_campos,
			      nsteps, nsnap);
	
	print_povray_header(b, cam_pos, camera_on_star_id,
			    aperture, gamma, blur_samples, nread, scale_L,
			    horizontal, vertical, print_hrd);

       if (base_type == StTr_station) {
	  // Plot the starbase at the origin.
	  cout << "#include \"starbase.inc\"" << endl;
	  cout << "object { Base_Station }\n" << endl;
	}
       else if (base_type == SW_ISD) {
	 cout << "object { StarWars_Scene }\n" << endl;
       }
	
	rdc_and_wrt_movie(b, true, scale_L, mmax, filter);
	last_frame++;

	if (v_flag)
	  put_dyn(cout, *b);
	
	rmtree(b);

	for (int i = 0; i < nint_skip; i++){
	  nsnap ++;
	  cerr << " skipping snapshot " << nsnap << endl;
	  if (! forget_node(cin)) exit(1);
	}

	if (nsnap==nsteps)
	  break;
	
    }

}

#endif

// endof: mkpovfile.C

