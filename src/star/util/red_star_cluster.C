/*
 *  red_stellar_system.C: reduce some useful information from the
 *  			 stellar system
 *.............................................................................
 *    version 1:  Jan 1997   Simon Portegies Zwart   email: spz@astro.uva.nl
 *.............................................................................
 *  non-local function: 
 *  all functions are local
 *.............................................................................
 *  Input options:
 *     -B:	Binning option of the Lagrangian radii.
 *     -C:	Cut off criterion, mass, luminosity, number of stars.
 *     -l:      Lower luminosity limit (see below).
 *     -n:	Number of Lagrangian radii bins.
 *		Might for some choices be forced.
 *     -o:	Output the story with data written in root.
 *     -S:      Sort option, on what parameter must information be sorted.
 *.............................................................................
 *  Concerning the Lower luminosity limit option (-l, see above).
 *    Currently this option is only applied to the sorting functions.
 *    This means that the binning (using Lagrangian radii) is applied
 *    on all cluster members.
 *    whether or not this is realistic depends on the application.
 *    Consequently, the luminosity cut-off does not affect the
 *    Lagrangian radii.
 *.............................................................................
 */

#include "stardyn_util.h"

#ifdef TOOLBOX

#define SCRATCH_PAD_LINE_LENGTH 255

char* stellar_type_summ_string = "ZAMS\tEarly_G\tLate_G\tHelium\tInert\tUnknown";

#define MAX_PREDEF_RADII 8
enum bin_option {lineair, predefined
};

enum sort_option {distance_to_com, stellar_mass, stellar_luminosity
};

enum cut_option {lagrangian_number_radii,
		 lagrangian_mass_radii, lagrangian_lumi_radii,
};

enum output_option {velocity_anisotropy, dispersion_velocity, number_of_stars,
		    number_of_star_types, total_mass_out, mass_over_light,
		    mean_mass, number_density, mass_density
};

static real to_kms=1;
static real to_myear=1;
static real to_msun=1;
static real to_rsun=1;
static real to_pc=1;

typedef  struct  {
  real sort_param;

  stellar_type_summary type;
  real age;
  real m_rel;
  real m_tot;
  real m_core;
  real T_eff;
  real L_eff;

  real  radius;  // sorted on this parameter.
  real  dmass;
  vector  vel;
  vector  pos;
} star_table, *star_table_ptr;

local bool check_input_option(bin_option binning, cut_option
			      radius_cut, sort_option sort_param, output_option option) {


  bool bin_choice  = FALSE;
  switch(binning) {
  case lineair: 
  case predefined: 
    bin_choice = true;
    break;
  case '?':   bin_choice = false;
  }

  bool cut_choice  = FALSE;
  switch(radius_cut) {
  case lagrangian_number_radii:
    cut_choice = TRUE;
    cerr << "Equal number of stars per bin.\n";
    break;
  case lagrangian_mass_radii:
    cut_choice = TRUE;
    cerr << "Equal mass per bin.\n";
    break;
  case lagrangian_lumi_radii:
    cut_choice = TRUE;
    cerr << "Equal luminosity per bin.\n";
    break;
  case '?':   cut_choice = FALSE;
  }

  bool sort_choice = false;
  switch(sort_param) {
  case distance_to_com: sort_choice = true;
    cerr << "Sorted on distance to cluster center.\n";
    break;
  case stellar_mass: sort_choice = true;
    cerr << "Sorted on stellar mass.\n";
    break;
  case stellar_luminosity: sort_choice = true;
    cerr << "Sorted on stellar luminosity.\n";
    break;
  case '?':   sort_choice = false;
  }

  bool output_choice = FALSE;
  switch(option) {
  case velocity_anisotropy:
  case dispersion_velocity:
  case number_of_stars:
  case number_of_star_types:
  case mass_over_light:
  case total_mass_out:
  case mean_mass:
    output_choice = TRUE;
    break;
  case '?':   output_choice = FALSE;
  }

  return (bin_choice && cut_choice && sort_choice);
}

local bool visible_object(star_table star, real l_min) {

  // Normally all stars with L>l_min are visible.
  // bool visible = (star.L_eff>l_min?1:0);

  bool visible = false;
  if (star.L_eff>=l_min
      && (int)star_type != (int)Inert_Remnant)	// casts added by Steve!!
    visible = true;

  // However, remnants are only visible if l_min is very small.
  // Neutron stars are intrinsically bright but extreme blue!
  // only if l_min==0, all stars are visible.
  if (l_min <=0 && star.type == Inert_Remnant)
    visible = true;

  return visible;
}

local real predifinced_lagrangian_fraction(int i) {

  real bin_radius;
  switch(i) {
  case 0: bin_radius = 0.02;
    break;
  case 1: bin_radius = 0.05;
    break;
  case 2: bin_radius = 0.10;
    break;
  case 3: bin_radius = 0.20;
    break;
  case 4: bin_radius = 0.30;
    break;
  case 5: bin_radius = 0.50;
    break;
  case 6: bin_radius = 0.70;
    break;
  case 7: bin_radius = 0.90;
    break;
  default: bin_radius = 1;
  }   

  return bin_radius;
}



local void accumulate_potential_energy(dyn* bj, dyn*bi,
	    			       real& epot, real& rmin, dynptr& bmin)

// Determine the potential energy of bi with respect to bj, along
// with bi's nearest neighbor and minimum neighbor distance.

{
    if (bj->get_oldest_daughter() != (dyn*)NULL)
	for (dyn * bb = bj->get_oldest_daughter(); bb != (dyn*)NULL;
	    bb = bb->get_younger_sister()) {
	    accumulate_potential_energy(bb, bi, epot, rmin, bmin);
	}
    else
	if (bi != bj) {
	    vector d_pos = bi->get_pos() - bj->get_pos();
	    real mi = bi->get_mass();
	    real mj = bj->get_mass();
	    real r = sqrt(d_pos * d_pos);
	    epot += -mi*mj/r;
	    if (r < rmin) {rmin = r; bmin = bj;}
	}
}

// compute_energies: calculate the appropriate color code for particle bi
//		     relative to "particle" bj (usually the root node).

local bool bound_object(dyn* bj, dyn* bi)
{
    bool bound = true;

    real   epot = 0, rmin = VERY_LARGE_NUMBER;
    dynptr bmin = bi;

    accumulate_potential_energy(bj, bi, epot, rmin, bmin);

    real   mi = bi->get_mass();
    real   mj = bmin->get_mass();
    vector d_vel = bi->get_vel() - bmin->get_vel();
    vector d_pos = bi->get_pos() - bmin->get_pos();
    real   r = sqrt(d_pos * d_pos);
    real   e0 = (0.5 * d_vel * d_vel - (mi + mj)/r);

    if (e0 < 0) {
	real   epot1 = 0, rmin1 = VERY_LARGE_NUMBER;
	dynptr bmin1 = bi;

	accumulate_potential_energy(bj, bmin, epot1, rmin1, bmin1);

	if (bi == bmin1) {
	    real e  = - mi*mj / r;
	    vector R_vel = (mi*bi->get_vel()+mj*bmin->get_vel())/(mi+mj);
	    real ekin = 0.5*(mi+mj)*R_vel*R_vel;

	    if (epot + epot1 - 2*e + ekin < 0) bound = true;
	    else bound = false;

	} else bound = true;

    } else {
	vector vel = bi->get_vel();
	real ekin = 0.5*bi->get_mass()*vel*vel;
	
	if (ekin + epot > 0.0) bound = false;
    }

    return bound;
}

void compute_cod_velocity(dyn *b, vector& vel)
{
  real total_weight = 0;
  int prev_n_bound, n_bound = 0;

  vector old_vel;
  real dvel = 0;
  do {
      
    old_vel = vel;
    vel = 0;

    prev_n_bound = n_bound;
    n_bound = 0;

    for_all_daughters(dyn, b, d) {

      real density = getrq(d->get_dyn_story(), "density");
      real dens2 = density * density;                        // weight factor

      if (bound_object(b, d)) {
	n_bound++;
	total_weight += dens2;
	vel += dens2 * d->get_vel();
      }
    }	
    vel /= total_weight;

    dvel = abs(old_vel - vel);
    
    b->set_vel(vel);
    cerr << " number of bound objects: " << n_bound <<endl;
    cerr << " density center velocity: " << vel << endl;
  }
  while(dvel>0.01);
  
//  while(prev_n_bound != n_bound);
	
}

local 
void compute_core_parameters(dyn* b, int k, vector& center, vector &vcore,
			     real& rcore, int& ncore, real& mcore)
{
    compute_density(b, k);
    compute_mean_cod(b, center, vcore);

    center -= b->get_pos();		// mean_cod quantities
    vcore -= b->get_vel();		// include the root node

    b->set_vel(vcore);
    compute_cod_velocity(b, vcore); 
    cerr << " density center position: " << center << endl;

    //    for_all_daughters(dyn, b, d)
    //      if (!bound_object(b, d))
    //	cerr << d->get_index() <<" "
    //	     << d->get_mass()*to_msun << " "
    //	     << abs(d->get_pos()-center) * to_pc << " "
    //	     << abs(d->get_vel()-vcore) *to_kms << endl;	

    real total_weight = 0;
    real sum = 0;
    for_all_daughters(dyn, b, bi) {
	real density = getrq(bi->get_dyn_story(), "density");
	real dens2 = density * density;                        // weight factor

        total_weight += dens2;
	sum += dens2 * square(bi->get_pos() - center);
    }

    real rcore2 = sum/total_weight;
    rcore = sqrt(rcore2);
    ncore = 0;
    mcore = 0;

    for_all_daughters(dyn, b, bj)
	if (square(bj->get_pos() - center) <= rcore2) {
	    ncore++;
	    mcore += bj->get_mass();
	}
}

local void print_core_parameters(dyn* b, vector& density_center,
				         vector &vcore,
				 real& rcore) {
				 
    real mcore;
    int ncore;
    
    compute_core_parameters(b, 6, density_center, vcore, rcore, ncore, mcore);
    //PRI(4); PRL(density_center);
    //PRI(4); PRC(rcore); PRC(ncore); PRL(mcore);

    //place the data in the root dyn-story.
    putvq(b->get_dyn_story(), "density_center_pos", density_center);
    putrq(b->get_dyn_story(), "core_radius", rcore);
    putiq(b->get_dyn_story(), "n_core", ncore);
    putrq(b->get_dyn_story(), "m_core", mcore);
}

local real lagrangian_radius(int i, int nzones,
                             bin_option binning) {

  real bin;
  if (binning == predefined && nzones==MAX_PREDEF_RADII) 
    bin = predifinced_lagrangian_fraction(i);
  else
    bin = ((i + 1) / (real)nzones);

  return bin;
}

local real* the_whole_cluster(dyn* b, star_table_ptr table, int n,
			      real l_min) {

  if (n==0)
    n = b->n_daughters();
  real *bin_param = new real[1];

  cerr << "\t" << Starlab::max(0.0, b->get_system_time()) << "\t" << 1 << endl;

  int i=n-1;
  while(!visible_object(table[i--], l_min));

  bin_param[0] = table[i].sort_param;

  //    cerr << "rn_1 = " << bin_param[0] << endl;

  return bin_param;
}

local real* put_lagrangian_number_radii(dyn* b, star_table_ptr table,
					int nzones, int n, real l_min,
					bin_option binning) {

  // return if unresonable input is given.
  if (nzones < 2) 
    return the_whole_cluster(b, table, n, l_min);
      

  if (n==0)
    n = b->n_daughters();

  real* bin_percent = new real[nzones];
  real* bin_param   = new real[nzones];

  int i;
  real total_number=0;
  for (i=0; i<n; i++)
    if (visible_object(table[i], l_min)) 
      total_number += 1;

    // determine the lagrangian radii.
  for (i=0; i<nzones; i++)
    bin_percent[i] = lagrangian_radius(i, nzones, binning)*total_number;

    // write the lagrangian parameter to the *dyn story.
  char r_n[3];
  sprintf(r_n, "rn_");


    // read back from *dyn the lagrangian parameter.
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];
  sprintf(scratch_pad, "n_number_zones = %d", nzones);
  cerr << "\t" << Starlab::max(0.0, b->get_system_time()) << "\t" << nzones << endl;

  real cumulative = 0;
  i = 0;
  for (int k=0; k<nzones; k++) {
    while (cumulative < bin_percent[k])  	
      if (visible_object(table[i++], l_min))
	cumulative += 1;			

    bin_param[k] = table[i-1].sort_param;
    sprintf(scratch_pad, "%s%d = %f", r_n, k+1, table[i-1].sort_param);
    add_story_line(b->get_log_story(), scratch_pad);
    //       cerr << r_n << k+1  << " = " << bin_param[k] << endl;
  }

  delete []bin_percent;
  
  //  if (scale)
  //     for (int k=0; k<nzones; k++) 
  //	 bin_param[k] *= to_pc;

  return bin_param;
}

local real* put_lagrangian_mass_radii(dyn* b, star_table_ptr table,
				      int nzones, int n, real l_min,
				      bin_option binning) {

  // return if unresonable input is given.
  if (nzones < 2)
    return the_whole_cluster(b, table, n, l_min);

  if (n==0)
    n = b->n_daughters();

  real* bin_percent = new real[nzones];
  real* bin_param   = new real[nzones];

  int i;
  real m_tot=0;
  for (i=0; i<n; i++)
    if (visible_object(table[i], l_min)) 
      m_tot += table[i].sort_param;
   
    // determine the lagrangian parameter.
  for (i=0; i<nzones; i++)
    bin_percent[i] = lagrangian_radius(i, nzones, binning) * m_tot;

    // write the lagrangiant radii to the *dyn story.
  char r_n[3];
  sprintf(r_n, "rm_");

    // read back from *dyn the lagrangian radii.
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];
  sprintf(scratch_pad, "n_mass_zones = %d", nzones);
  cerr << "\t" << Starlab::max(0.0, b->get_system_time()) << "\t" << nzones << endl;

  real cumulative = 0;
  i = 0;
  for (int k=0; k<nzones; k++) {
    while (cumulative < bin_percent[k]) {
      if (visible_object(table[i], l_min)) 
	cumulative += table[i].sort_param;
      i++;
    };

    bin_param[k] = table[i-1].sort_param;
    sprintf(scratch_pad, "%s%d = %f", r_n, k+1, table[i-1].sort_param);
    add_story_line(b->get_log_story(), scratch_pad);
    //       cerr << r_n << k+1  << " = " << bin_param[k] << endl;
  }

  delete []bin_percent;

  //  if (scale)
  //     for (int k=0; k<nzones; k++) 
  //	 bin_param[k] *= to_msun;

  return bin_param;
}

local real* put_lagrangian_lumi_radii(dyn* b, star_table_ptr table,
				      int nzones, int n, real l_min,
				      bin_option binning) {

  // return if unresonable input is given.
  if (nzones < 2)
    return the_whole_cluster(b, table, n, l_min);

  if (n==0)
    n = b->n_daughters();

  real* bin_percent = new real[nzones];
  real* bin_param   = new real[nzones];

  int i;
  real total_lumi=0;
  for (i=0; i<n; i++)
    if (visible_object(table[i], l_min))
      total_lumi += table[i].sort_param;
   
    // determine the lagrangian parameters.
  for (i=0; i<nzones; i++)
    bin_percent[i] = lagrangian_radius(i, nzones, binning) * total_lumi;

    // write the lagrangiant parameter to the *dyn story.
  char r_n[3];
  sprintf(r_n, "rl_");

    // read back from *dyn the lagrangian radii.
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];
  sprintf(scratch_pad, "n_lumi_zones = %d", nzones);
  cerr << "\t" << Starlab::max(0.0, b->get_system_time()) << "\t" << nzones << endl;

  real cumulative = 0.0;
  i = 0;
  for (int k=0; k<nzones; k++) {
    while (cumulative < bin_percent[k]) {
      if (visible_object(table[i], l_min)) 
	cumulative += table[i].L_eff;
      i++;
    };

    bin_param[k] = table[i-1].sort_param;
    sprintf(scratch_pad, "%s%d = %f", r_n, k+1, table[i-1].sort_param);
    add_story_line(b->get_log_story(), scratch_pad);
    //       cerr << r_n << k+1  << " = " << table[i-1].sort_param << endl;
  }

  delete []bin_percent;

  //  if (scale)
  //     for (int k=0; k<nzones; k++) 
  //	 bin_param[k] *= to_pc;

  return bin_param;
}

local void extr_velocity_anisotopy(star_table_ptr table,
				   int n, real zone, int nstar, real l_min,
				   bool verbatim) {
  real v2, vr, vr2;
  real total_weight =0;
  real vr2_sum = 0;
  real vt2_sum = 0;

  real weight = 1;			// Heggie/Binney & Tremaine
  int i = nstar;
  int n_stars = 0;
  while (table[i++].sort_param < zone) 
    if (visible_object(table[i], l_min)) {
      n_stars++;

      //weight = bi->get_mass();

      v2 = table[i].vel * table[i].vel;
      vr = table[i].pos * table[i].vel;
      vr2 = vr * vr / square(table[i].pos);

      total_weight += weight;
      vr2_sum += weight * vr2;
      vt2_sum += weight * (v2 - vr2);
    }

  if (i>nstar) 
      cerr << 1 - 0.5 * vt2_sum / vr2_sum << " ";
  else
    cerr << 0 << " ";
}

local void extr_dispersion_velocity(star_table_ptr table, int n,
				    real* zones, int nzones, real l_min,
				    bool scale_flag,
				    bool verbatim = true) {

  int n_stars;
  real v2_cum, v_cum;
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];

  if (!scale_flag)
    cerr << "Zone \t\tdispersion velocity <v^2>^{1/2}:"<<endl;
  else
    cerr << "Zone \t\tdispersion velocity <v^2>^{1/2} [km/s]:"<<endl;

  int i = 0;
  for (int k=0; k<nzones; k++) {
    
    n_stars = 0;
    v2_cum = 0;
    while (table[i++].sort_param < zones[k]) 
      if (visible_object(table[i], l_min)) {
	v2_cum += pow(abs(table[i].vel), 2);
	n_stars++;
      }

    v_cum=0;
    if (n_stars>=1) {
      v_cum = sqrt(v2_cum/dynamic_cast(int, n_stars));
    }

    //    if (scale_flag) {
    //      v_cum *= to_kms;	
    //    }

    sprintf(scratch_pad, "v_disp = %lf", v_cum);
    if (scale_flag)
      cerr << k+1 << " " << zones[k]*to_pc
	   << " \t" << v_cum*to_kms << endl;
    else
      cerr << k+1 << " " << zones[k]
	   << " " << v_cum << endl;

  }           
}
local void extr_dispersion_velocity(star_table_ptr table, int n,
				    real zone, int nstar, real l_min,
				    bool scale_flag,
				    bool verbatim = true) {

  int n_stars;
  real v2_cum, v_cum;
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];

  int i = nstar;
    
    n_stars = 0;
    v2_cum = 0;
    while (table[i++].sort_param < zone) 
      if (visible_object(table[i], l_min)) {
	v2_cum += pow(abs(table[i].vel), 2);
	n_stars++;
      }

    v_cum=0;
    if (n_stars>=1) {
      v_cum = sqrt(v2_cum/dynamic_cast(int, n_stars));
    }

    if (scale_flag)
      cerr << v_cum*to_kms << " ";
    else
      cerr << v_cum << " ";

  }           

local void count_stars(star_table_ptr table, int n,
        	       real zone, int nstar, real l_min,
		       bool scale_flag = true,
		       bool verbatim = true) {

  int i = nstar;
  int n_stars = 0;
  while (table[i++].sort_param < zone) 
    if (visible_object(table[i], l_min)) 
      n_stars++;

  if (scale_flag)
	cerr << n_stars << " ";
  else
	cerr << n_stars << " ";

}           

local void count_stars(star_table_ptr table, int n,
        	       real* zones, int nzones, real l_min,
		       bool scale_flag = true,
		       bool verbatim = true) {

  int n_stars;
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];

  cerr << "Zone \t\tnumber of stars:"<<endl;


  int i = 0;
  for (int k=0; k<nzones; k++) {

    n_stars = 0;
    while (table[i++].sort_param < zones[k]) 
      if (visible_object(table[i], l_min)) 
	n_stars++;

   
    sprintf(scratch_pad, "n_stars = %d", n_stars);
    if (scale_flag)
      cerr << k+1 << " " << zones[k]*to_pc
	   << " \t" << n_stars << endl;
    else
      cerr << k+1 << "  " << zones[k]
	   << " " << n_stars << endl;

  }           
}

local void count_stellar_types(star_table_ptr table, int n,
			       real* zones, int nzones,
			       real l_min,
			       bool scale_flag = true, 
			       bool verbatim = true) {

  int* types = new int[no_of_star_type_summ];

  int i = 0;
  cerr << "Zone \t\t" << stellar_type_summ_string << endl;
  for (int k=0; k<nzones; k++) {
    while (table[i++].sort_param < zones[k]) 
      if (visible_object(table[i], l_min)) 
	types[(int)table[i].type]++;

    if (scale_flag)	
      cerr << k+1 << " " << zones[k]*to_pc
	   << " \t";
    else
      cerr << k+1 << " " << zones[k]
	   << " \t";

    for (int j=0; j<no_of_star_type_summ-1; j++) {
      cerr << types[j] << "\t";
      types[j] = 0;
    }
    cerr << endl;
  }

  delete []types;
}

// the luminosity limit affects the cumulative luminosity, not the mass!
local int mass_over_light_ratio(star_table_ptr table, int n,
        		         real zone, int nstar, real l_min,
				 bool verbatim = true) {

  real ml_ratio;
  real mass, light;

  int i = nstar;
  real m_cum, l_cum;

  m_cum = l_cum = 0;
    while (table[i++].sort_param < zone) {
      mass = table[i].m_tot;
      light = table[i].L_eff;
      if (visible_object(table[i], l_min)) 
	l_cum += light;
      m_cum += mass;
    };

    if (l_cum>0)
      ml_ratio = m_cum/l_cum;
    else
      ml_ratio = -1;
   
      cerr << ml_ratio << endl;

    return i - nstar;
  }           

local void total_stellar_mass(star_table_ptr table, int n, real* zones,
                              int nzones, real l_min,
			      bool scale_flag = true,
			      bool verbatim = true) {

  real m_tot;
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];

  cerr << "Zone \ttotal mass:"<<endl;

  int i = 0;
  for (int k=0; k<nzones; k++) {

    m_tot = 0;
    while (table[i++].sort_param < zones[k]) 
      if (visible_object(table[i], l_min)) 
	m_tot += table[i].m_tot;
   
    sprintf(scratch_pad, "m_tot = %f", m_tot);
    if (scale_flag)
      cerr << k+1 << " " << zones[k]*to_pc
	   << " \t" << m_tot << endl;
    else
      cerr << k+1 << " " << zones[k]
	   << " \t" << m_tot/to_msun << endl;
  }           
}
local void mean_stellar_mass(star_table_ptr table, int n, real zone,
                             int nstar, real l_min,
			     bool scale_flag, bool verbatim = true) {

  real m_tot, n_star, mmass;

  int i = nstar;
    m_tot = n_star = 0;
    while (table[i++].sort_param < zone) 
      if (visible_object(table[i], l_min)) {
	m_tot += table[i].m_tot;
	n_star += 1;
      }

    mmass = m_tot / n_star;
    if (!scale_flag)
      mmass /= to_msun;
     
    if (scale_flag)
      cerr << mmass << " ";
    else
      cerr << mmass/to_msun << " ";
  }           

local void mean_stellar_mass(star_table_ptr table, int n, real* zones,
                             int nzones, real l_min,
			     bool scale_flag, bool verbatim = true) {

  real m_tot, n_star, mmass;
  char scratch_pad[SCRATCH_PAD_LINE_LENGTH];

    if (scale_flag)
      cerr << "Zone \tmean mass: "<<endl;
    else
      cerr << "Zone \tmean mass [msun]:"<<endl;


  int i = 0;
  for (int k=0; k<nzones; k++) {

    m_tot = n_star = 0;
    while (table[i++].sort_param < zones[k]) 
      if (visible_object(table[i], l_min)) {
	m_tot += table[i].m_tot;
	n_star += 1;
      }

    mmass = m_tot / n_star;
    if (!scale_flag)
      mmass /= to_msun;
     
    sprintf(scratch_pad, "m_tot = %f", mmass);
    if (scale_flag)
      cerr << k+1 << " " << zones[k]*to_pc
	   << " \t" << mmass << endl;
    else
      cerr << k+1 << " " << zones[k]
	   << " \t" << mmass/to_msun << endl;
  }           
}

local void stellar_number_density(star_table_ptr table, int n,
                                  real zone1, real zone2, 
				  int nstar, real l_min,
				  bool scale_flag, bool verbatim = true) {

  real m_tot, n_star;
  real area, density;
  real to_area = 4*cnsts.mathematics(pi);

  real to_kubic_parsec = 1./pow(to_rsun/4.4e+7, 3);

  int i = nstar;
    m_tot = n_star = 0;
    while (table[i++].sort_param < zone2) 
      if (visible_object(table[i], l_min)) {
	m_tot  += table[i].m_tot;
	n_star += 1;
      }

    if (zone1>0)
      area = to_area * pow(zone1, 2) * (zone2 - zone1);
    else
      area = cnsts.mathematics(one_third) * to_area * pow(zone2, 3);

    density = n_star/area;

    if (scale_flag)
      cerr << density*to_kubic_parsec << " ";
    else
      cerr << density << " ";
  }           

local void stellar_mass_density(star_table_ptr table, int n,
				real zone1, real zone2,
				int nstar, real l_min,
				bool scale_flag, bool verbatim = true) {

  real m_tot, n_star;
  real area, density;
  real to_area = 4*cnsts.mathematics(pi);

  real to_kubic_parsec = 1./pow(to_rsun/4.4e+7, 3);

  int i = nstar;

    m_tot = n_star = 0;
    while (table[i++].sort_param < zone2) 
      if (visible_object(table[i], l_min)) {
	m_tot  += table[i].m_tot;
	n_star += 1;
      }

    if (zone1>0)
      area = to_area * pow(zone1, 2) * (zone2 - zone1);
    else
      area = cnsts.mathematics(one_third) * to_area * pow(zone2, 3);

    density = m_tot/area;

    if(scale_flag)
      cerr << density*to_kubic_parsec << " ";
    else
      cerr << density << " ";
  }           


//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)
{
  if (((star_table_ptr) pi)->radius > ((star_table_ptr) pj)->radius)
    return(1);
  else if (((star_table_ptr)pi)->radius < ((star_table_ptr)pj)->radius)
    return(-1);
  else
    return(0);
}

//-----------------------------------------------------------------------------
//  compare_parameters  --  compare parameters of two particles
//-----------------------------------------------------------------------------

local int compare_parameters(const void * pi, const void * pj)
{
  if (((star_table_ptr) pi)->sort_param > ((star_table_ptr) pj)->sort_param)
    return(1);
  else if (((star_table_ptr)pi)->sort_param < ((star_table_ptr)pj)->sort_param)
    return(-1);
  else
    return(0);
}

local real sorted_parameter(star_table table, sort_option sort_param) {

  real parameter;
  if (sort_param==distance_to_com) {
    parameter = table.radius;
  }
  if (sort_param==stellar_mass) {
    parameter = table.dmass;
  }
  else 
    parameter = table.L_eff;

  return parameter;
}

//-----------------------------------------------------------------------------
//  
//-----------------------------------------------------------------------------
local star_table_ptr sort_stellar_parameter(dyn * b, int nzones, int &n,
					    sort_option sort_param) {

  stellar_type type=NAS;
  story* st=NULL;
  real t_cur, t_rel, m_rel, m_env, m_core, co_core, t_eff, l_eff, p_rot, b_fld;

  if (n==0)
    n = b->n_daughters();

  star_table_ptr table = new star_table[n];
  if (table == NULL) {
    cerr << "sort_stellar_parameter: "
	 << "not enough memory left for table\n";
    exit(0);
  }

  // Use the geometric center if the density center is unknown.
  vector dc_pos = 0;
  vector dc_vel = 0;
  real rcore = 0;
  if (find_qmatch(b->get_dyn_story(), "density_center_pos"))
    dc_pos = getvq(b->get_dyn_story(), "density_center_pos");
  else
    print_core_parameters(b, dc_pos, dc_vel, rcore);
  
  int i;

  // fill table with usefull information.
  // note that this part might needs a flattening of the tree.
  dyn* bi;
  n=0;
  for (i = 0, bi = b->get_oldest_daughter(); bi != NULL;
       bi = bi->get_younger_sister()) {

    if (bound_object(b, bi)) {

    table[i].radius = abs(bi->get_pos() - dc_pos);
    table[i].dmass = bi->get_mass();

    st = bi->get_starbase()->get_star_story();
    extract_story_chapter(type, t_cur, t_rel, 
		      m_rel, m_env, m_core,
		      co_core,
		      t_eff, l_eff,
		      p_rot, b_fld,
		      *st);
    
    if (summarize_stellar_type(type)==Inert_Remnant)
       l_eff = 0;
    table[i].type = summarize_stellar_type(type);
    table[i].age = t_rel;
    table[i].m_rel = m_rel;
    table[i].m_tot = m_env+m_core;
    table[i].m_core = m_core;
	// XXX co_core too?
    table[i].T_eff = t_eff;
    table[i].L_eff = l_eff;
    table[i].vel = bi->get_vel() -dc_vel; 
    table[i].pos = bi->get_pos() - dc_pos; 

    /*
    table[i].type = (stellar_type_summary)1;
    table[i].age = table[i].dmass; 
    table[i].m_rel = table[i].dmass; 
    table[i].m_tot = table[i].dmass; 
    table[i].m_core = table[i].dmass; 
    table[i].T_eff = table[i].dmass; 
    table[i].L_eff = table[i].dmass; 
    */
    
    table[i].sort_param = sorted_parameter(table[i], sort_param);
    i++;
  }
  }
  n=i;
  qsort((void *)table, (size_t)n, sizeof(star_table), compare_parameters);

    // now the array table is sorted in radius.
    // we can do all kind of beautiful things with it.

  return table;
}

local void red_stellar_system(dyn* b, int nzones,
			      bin_option binning, cut_option
			      radius_cut, output_option option,
			      sort_option sort_param, real l_min,
			      bool scale_flag = false,
			      bool verbatim = true) {

  int n = b->n_daughters();

  star_table_ptr table = sort_stellar_parameter(b, nzones, n, sort_param);

  real* zones;
  if (nzones < 2) {
    zones = the_whole_cluster(b, table, n, l_min);
    nzones = 1;
  }
  else
    switch(radius_cut) {
    case lagrangian_number_radii:
      zones = put_lagrangian_number_radii(b, table,
					  nzones, n, l_min,
   				          binning);
      break;
    case lagrangian_mass_radii:
      zones = put_lagrangian_mass_radii(b, table, nzones,
					n, l_min,
	                                binning);
      break;
    case lagrangian_lumi_radii:
      zones = put_lagrangian_lumi_radii(b, table, nzones,
					n, l_min, binning);
      break;
    case '?':   cerr << "no option specified in red_stellar_system.";
      cerr << "use whole cluster instead.";
    default:    zones = the_whole_cluster(b, table, n, l_min);
      nzones = 1;
    };
  //       put_lagrangian_parameters(b, table,   nzones, n, l_min, binning);


  cerr << "i    R   N   <m>    n     rho    <v>    vt/vr    M/L" << endl;
  cerr << "    [pc]    [Msun]   [*/pc^3]      [km/s]        [Msun/Lsun]"
       << endl;
  
  int ncount, nstar = 0;
  for (int k=0; k<nzones; k++) {
    if (scale_flag)
      cerr << k+1 << "  " << zones[k]*to_pc << " ";
    else
      cerr << k+1 << "  " << zones[k] << " ";

    count_stars(table, n, zones[k], nstar,
		l_min, scale_flag, verbatim);
    mean_stellar_mass(table, n, zones[k], nstar,
		      l_min, scale_flag, verbatim);
    if (k==0) {
      stellar_number_density(table, n, 0, zones[k], nstar,
			     l_min, scale_flag, verbatim);
      stellar_mass_density(table, n, 0, zones[k], nstar,
			   l_min, scale_flag, verbatim);
    }
    else {
      stellar_number_density(table, n, zones[k-1], zones[k], nstar,
			     l_min, scale_flag, verbatim);
      stellar_mass_density(table, n, zones[k-1], zones[k], nstar,
			   l_min, scale_flag, verbatim);
    }
    extr_dispersion_velocity(table, n, zones[k], nstar,
			     l_min, scale_flag, verbatim);
    extr_velocity_anisotopy(table, n, zones[k], nstar,
			    l_min, verbatim);
    ncount = mass_over_light_ratio(table, n, zones[k], nstar,
			  l_min, verbatim);
    nstar += ncount;
  }

  //Print zero's to allow automatic ploting program to take over.
  cerr << "0 0 0 0 0 0 0 0 0" <<endl;
    
#if 0
  switch(option) {
  case velocity_anisotropy:
    extr_velocity_anisotopy(table, n, zones, nzones,
			    l_min, verbatim);
    break;
  case dispersion_velocity:
    extr_dispersion_velocity(table, n, zones, nzones,
			     l_min, scale_flag, verbatim);
    break;
  case number_of_stars:
    count_stars(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  case number_of_star_types:
    count_stellar_types(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  case mass_over_light:
    mass_over_light_ratio(table, n, zones, nzones,
			  l_min, scale_flag, verbatim);
    break;
  case total_mass_out:
    total_stellar_mass(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  case mean_mass:
    mean_stellar_mass(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  case number_density:
    stellar_number_density(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  case mass_density:
    stellar_mass_density(table, n, zones, nzones, l_min, scale_flag, verbatim);
    break;
  default:    cerr << "no option specified in red_stellar_system.";
    break;
  };
#endif

  delete []table;
  delete []zones;

}

local void project(dyn *b, int axis) {

  vector new_pos, new_vel;
  for_all_leaves(dyn, b, bi) {
      new_pos = bi->get_pos();
      new_vel = bi->get_vel();
      new_pos[axis] = 0;
      new_vel[axis] = 0;
      bi->set_pos(new_pos);
      bi->set_vel(new_vel);
  }
}

local void print_parameter_usage() {

  cerr << "Parameters usage for: red_stellar_system." << endl;
  cerr << "options:\n"
       << "     -B:  Binning option of the Lagrangian radii.\n"
       << "          0= lineair, 1= predefined.\n"
       << "     -C:  Cut off creterium, mass, luminosity, number of stars.\n"
       << "          0= by number, 1= by mass, 2= by luminosity.\n"
       << "     -l:  Lower luminosity limit.\n"
       << "     -N:  Number of Lagrangian radii bins.\n"
       << "	       Might for some choices be forced.\n"
       << "     -n:  Allow for n-squared operations.\n"
       << "     -o:  Output the story with data written in root.\n"
       << "     -S:  Sort option, on which parameter should be sorted.\n"
       << "          0= distance to com, 1= stellar mass.\n"
       << "     -s:  Scale to astrophysical units.\n"
       << "     -v:  give extended output\n"
       << endl;
}

/*-----------------------------------------------------------------------------
 *  main  --  driver to use  compute_mass_radii() as a tool
 *-----------------------------------------------------------------------------
 */
main(int argc, char ** argv)
{
  char  *comment;
  bool out = false;
  bool verbatim = false;
  bool nsq_opt = false;
  bool scale_flag = false;
  bin_option binning = lineair;
  cut_option radius_cut = lagrangian_mass_radii;
  output_option option = velocity_anisotropy;
  sort_option sort_param = distance_to_com;
  real l_min = 0.0;		// minimum visible luminosity in Lsun.
  int nzones = 50;
  bool  c_flag = FALSE;      /* if TRUE: a comment given on command line   */
  int axis = -1;
  
  extern char *poptarg;
  int c;
  char* param_string = "a:B:C:c:l:N:noS:sv";

  while ((c = pgetopt(argc, argv, param_string)) != -1)
    switch(c)
      {
      case 'a': axis = atoi(poptarg);
	        break;
      case 'B': binning = (bin_option)atoi(poptarg);
	        break;
      case 'C': radius_cut = (cut_option)atoi(poptarg);
	        break;
      case 'c': c_flag = TRUE;
	        comment = poptarg;
	        break;
      case 'l': l_min = atof(poptarg);
	        break;
      case 'N': nzones = atoi(poptarg);
	        break;
      case 'n': nsq_opt = true;
	        break;
      case 'o': out = true;
	        break;
      case 'S': sort_param = (sort_option)atoi(poptarg);
	        break;
      case 's': scale_flag = true;
	        break;
      case 'v': verbatim = true;
	        break;
      case '?': params_to_usage(cerr, argv[0], param_string);
	        print_parameter_usage();
	        exit(1);
      }            

  dyn *b;

  if (!check_input_option(binning, radius_cut, sort_param, option)) {
    cerr << "input option not recognized!" << endl;
    print_parameter_usage();
    exit(1);
  }
  if (binning==predefined)
    nzones = MAX_PREDEF_RADII;

  cerr.precision(LOW_PRECISION);
  vector density_center=0;
  vector density_velocity=0;
  real rcore;
  while (b = get_dyn()) {
    
    if (c_flag == TRUE)
      b->log_comment(comment);

    b->log_history(argc, argv);

    to_myear= 1./b->get_starbase()->conv_t_star_to_dyn(1);
    to_msun = 1./b->get_starbase()->conv_m_star_to_dyn(1);
    to_rsun = 1./b->get_starbase()->conv_r_star_to_dyn(1);
    to_pc   = to_rsun * cnsts.parameters(Rsun)/cnsts.parameters(parsec);
    to_kms  = (to_rsun* cnsts.parameters(Rsun)
            /           cnsts.physics(kilometer))
            / (to_myear*cnsts.physics(Myear));

    //cerr << "scalings" << to_myear << " "
    //                   << to_msun  << " "
    //			 << to_rsun  << " "
    //		         << to_kms   << endl;

    cerr << "Time = " << b->get_system_time() << endl;
    
    if (nsq_opt)
      print_core_parameters(b, density_center, density_velocity, rcore);

    if(axis>=0)
      project(b, 0);

    red_stellar_system(b, nzones, binning, radius_cut, option,
		       sort_param, l_min, scale_flag, verbatim);

    if (out) put_dyn(b);
    delete b;
  }
}

#endif
