
//  sstar_stats.C: Print out various diagnostic statistics on a system.

#include "sstar_to_dyn.h"
#include "dyn.h"
#include "single_star.h"

#ifndef TOOLBOX

#define OUTPUT_MASS_LIMIT   10		// Msun
#define OUTPUT_NUMBER_LIMIT 10
#define MINIMUM_LUMINOSITY 1.0e-4	// planets' luminosity

#define MAX_PREDEF_RADII 9
enum bin_option {lineair, predefined};

typedef  struct  {

  real sort_param;

  stellar_type_summary type;
  real Lu;
  real Lb;
  real Lv;
  real Lr;
  real Li;
} star_table, *star_ubvri_ptr;

struct pop {
    int stypes[no_of_stellar_type];

    pop() {
	for (int i=0; i<no_of_stellar_type; i++)
	    stypes[i] = 0;
    }
};

struct nm_bin {
    int   no_of_stars;
    real  mass_limit;
    real  total_mass;
    nm_bin() {
	no_of_stars = 0;
	mass_limit=total_mass = 0;
    }
};

local int which_zone(dyn* bi, vector& center, int n_lagr, real* r_lagr)
{
    vector pos = something_relative_to_root(bi, &dyn::get_pos);
    vector dr = pos - center;
    real dr2 = dr*dr;

    for (int j = 0; j < n_lagr; j++)
	if (r_lagr[j] > dr2) return j;

    return n_lagr;
}

local void print_stellar_content(dyn* b) {

    int n_lagr = 0;
    real *r_lagr = NULL;
    vector lagr_pos = 0;

    // As in lagrad, use the geometric center if 
    // density center is unknown.

    if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

	// Make a list of previously computed lagrangian radii.

	n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	r_lagr = new real[n_lagr];
	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	lagr_pos = getvq(b->get_dyn_story(), "lagr_pos");

    }
    else 
	r_lagr = new real[n_lagr];

    pop *content = new pop[n_lagr+1];

    // Note that r_lagr omits the 0% and 100% radii.
    
    // Initialization:

    for (int i = 0; i <= n_lagr; i++)
	if (i < n_lagr) r_lagr[i] *= r_lagr[i];

    int N_total = 0;
    stellar_type stype = NAS;

    for_all_daughters(dyn, b, bi) {

	// Count single stars and binaries separately.

	if (bi->get_oldest_daughter() != NULL)
	
	    stype = Double;
    
	else {

	    if (b->get_use_sstar()) 
		stype = bi->get_starbase()->get_element_type();
	    else if (find_qmatch(bi->get_star_story(), "Type")) 
		stype = extract_stellar_type_string(
					getsq(bi->get_star_story(), "Type")); 
	    else {
		cerr << "    No stellar information found for: ";
		bi->pretty_print_node(cerr);
		return;
	    }
	}

	// Find which zone we are in.

	int i = which_zone(bi, lagr_pos, n_lagr, r_lagr);

	// Update statistics.

	N_total++;
	content[i].stypes[(int)stype]++;
    }

    if (N_total > 0) {

	int j, k;
	char* ts;
	cerr << "           ("; 
	for (j = 0; j <= dynamic_cast(int, Super_Giant); j++) 
	    cerr << " " << type_short_string(dynamic_cast(stellar_type, j));

	cerr << "     ";
	for (j = dynamic_cast(int, Carbon_Star); 
	     j < dynamic_cast(int, no_of_stellar_type); j++) 
	    cerr << " " << type_short_string(dynamic_cast(stellar_type, j));

	cerr << ")" << endl;
      
	for (k = 0; k <= n_lagr; k += 5) {
	    if (k > 0) cerr << endl;
	    for (int i = k; i < k+5 && i <= n_lagr; i++) {
		cerr << "    zone " << i << " ";
		for (j=0; j<=Super_Giant; j++)
		    cerr << " " << content[i].stypes[j];
		cerr << "     ";
		for (j=Carbon_Star; j<no_of_stellar_type; j++)
		    cerr << " " << content[i].stypes[j];
		cerr << endl;
	    }
	}

    } else
	cerr << "            (none)\n";
  
    delete [] content;
    delete [] r_lagr;

//  }
//  else 
//      cerr << "     ---No Lagrangian radii specified--- " << endl;
  
}

local int print_special_stars(dyn* b) {

    if (!b->get_use_sstar())
	return -1;

//    real M_to = turn_off_mass(b->get_starbase()
//			       ->conv_t_dyn_to_star(b->get_system_time()));

    real current_time = b->get_starbase()
                         ->conv_t_dyn_to_star(b->get_system_time());
    int N_special = 0;
    real f_Bs, t_Bs;
    for_all_leaves(dyn, b, bi) {

	// Count all single stars 

	stellar_type stype = NAS;
	star* ss = NULL;

	if (has_sstar(bi)) {
	    star* ss = dynamic_cast(star*, bi->get_starbase());

	    if (ss!=NULL) {
		star_state sts(ss);
      
		if (sts.special()) {

		    N_special++;

		    cerr << "    " << bi->format_label() 
			 << " is ";
		    put_short_state(sts, cerr);
		    cerr << " (" << type_dominant_state(sts)
			 << ")";
	  
		    if (sts.class_spec[Blue_Straggler]) {
			//f_Bs = sts.mass/M_to;
			//cerr << " [f_Bs = " << f_Bs << "]";
		        t_Bs = current_time - ss->get_relative_age();
			cerr << " [t_Bs = " << t_Bs << " Myr]";
		    }
	    
		    if (ss->is_binary_component()) {
			cerr << " binary companion is "
			    << bi->get_binary_sister()->format_label() << " ";
			print_star(dynamic_cast(starbase*,
						ss->get_companion()));
		    }
		    cerr << endl;
		}
		else if (ss->get_element_type() == Disintegrated) {

		    N_special++;

		    cerr << "    " << bi->format_label() 
			 << " is disintegrated star" << endl;
		}
	    }
	}
    }
  
    if (N_special == 0)
	cerr << "     ---No special stars---" <<endl;
}

local void print_compact_stars(dyn* b) {

    //  for_all_daughters(dyn, b, bi) {
    if (!b->get_use_sstar())
	return;

    int N_compact = 0;

    for_all_leaves(dyn, b, bi) {

	// Count all single stars 

	stellar_type stype = NAS;
	star* ss = NULL;

	// if (bi->get_oldest_daughter() != NULL)
	if (has_sstar(bi)) {
      
	    ss = dynamic_cast(star*, bi->get_starbase());

	    bool accreting = false;
	    if (ss!=NULL) 
		if(ss->remnant()) {
		    N_compact++;

		    stype = ss->get_element_type();
		    accreting = ((ss->get_spec_type(Accreting))?true:false);

		    if (stype == Xray_Pulsar || stype == Radio_Pulsar) {
	  
			cerr << "    " << bi->format_label()
			     << " is " << type_string(stype)
			     << " (Ps = " << ss->get_rotation_period()
			     << ", log B = " << ss->get_magnetic_field()
			     << ") ";
			if (ss->is_binary_component()) {
			    cerr << " binary companion is "
				 << bi->get_binary_sister()->format_label()
				 << " ";
			    print_star(dynamic_cast(starbase*,
						    ss->get_companion()));
			}
			cerr << endl;
		    }
		    else if (stype == Helium_Giant) {
			
			cerr << "    " << bi->format_label()
			     << " is " << type_string(stype)
			     << " (M =" << ss->get_relative_mass()
			     << ", Reff = " << ss->get_effective_radius()
			     << ", L = " << ss->get_luminosity()
			     << ") ";
			if (ss->is_binary_component()) {
			    cerr << "binary companion is "
				 << bi->get_binary_sister()->format_label()
				 << " ";
			    print_star(dynamic_cast(starbase*,
						    ss->get_companion()));
			}
			cerr << endl;
		    }
		    else if (accreting) {

			cerr << "    " << bi->format_label()
			     << " is accreting " << type_string(stype);
			if (ss->is_binary_component()) {
			    cerr << " binary companion is "
				 << bi->get_binary_sister()->format_label()
				 << " ";
			    print_star(dynamic_cast(starbase*,
						    ss->get_companion()));
			}
			cerr << endl;
		    }
		}
		else if (stype == Thorn_Zytkow) {

		    cerr << "    " << bi->format_label()
			 << " is " << type_string(stype)
			 << " (M =" << ss->get_relative_mass()
			 << ", Reff = " << ss->get_effective_radius()
			 << ", L = " << ss->get_luminosity()
			 << ") ";
		    if (ss->is_binary_component()) {
			cerr << "binary companion is "
			     << bi->get_binary_sister()->format_label()
			     << " ";
			print_star(dynamic_cast(starbase*,
						ss->get_companion()));
		    }
		    cerr << endl;
		}
	}
    }
  
    if (N_compact == 0)
	cerr << "     ---No compact stars---" <<endl;
    else
	cerr << "        Number of compact stars = " << N_compact << endl;
}


local int compare_parameters(const void * pi, const void * pj) {

  if (((star_ubvri_ptr) pi)->sort_param > ((star_ubvri_ptr) pj)->sort_param)
    return(1);
  else if (((star_ubvri_ptr)pi)->sort_param < ((star_ubvri_ptr)pj)->sort_param)
    return(-1);
  else
    return(0);
}

local void sort_stellar_ubvri(dyn * b, int nzones, int n,
			      star_ubvri_ptr table, real *Ltot) {

  // Use the geometric center if the density center is unknown.
  vector dc_pos = 0;
  vector dc_vel = 0;
  real rcore = 0;
  if (find_qmatch(b->get_dyn_story(), "density_center_pos"))
    dc_pos = getvq(b->get_dyn_story(), "density_center_pos");
//  else
//    print_core_parameters(b, dc_pos, dc_vel, rcore);

  int i=0;

  // fill table with usefull information.
  // note that this part might needs a flattening of the tree.
  stellar_type stype;
  real Lu, Lb, Lv, Lr, Li;

  Ltot[0] = Ltot[1] = Ltot[2] = Ltot[3] = Ltot[4] = 0;
	
  vector rstar;
  for_all_leaves(dyn, b, bi) {

      rstar = (bi->get_top_level_node()->get_pos() - dc_pos);
      table[i].sort_param = abs(rstar[0]);

      stype = NAS;

      Lu=Lb=Lv=Lr=Li=0;
      get_Lubvri_star(bi, stype, Lu, Lb, Lv, Lr, Li);

    
      if (summarize_stellar_type(stype)==Inert_Remnant)
	  Lu=Lb=Lv=Lr=Li=0;

      table[i].type = summarize_stellar_type(stype);
      table[i].Lu = Lu;
      table[i].Lb = Lb;
      table[i].Lv = Lv;
      table[i].Lr = Lr;
      table[i].Li = Li;

//      Ltot[0] += table[i].Lu;
//      Ltot[1] += table[i].Lb;
//      Ltot[2] += table[i].Lv;
//      Ltot[3] += table[i].Lr;
//      Ltot[4] += table[i].Li;
      
      i++;
  }

  qsort((void *)table, (size_t)n, sizeof(star_table), compare_parameters);

  for(i=1; i<n; i++) {
      table[i].Lu += table[i-1].Lu;
      table[i].Lb += table[i-1].Lb;
      table[i].Lv += table[i-1].Lv;
      table[i].Lr += table[i-1].Lr;
      table[i].Li += table[i-1].Li;
  }
  for(i=0; i<n; i++) {
      table[i].Lu = log(table[i].Lu);
      table[i].Lb = log(table[i].Lb);
      table[i].Lv = log(table[i].Lv);
      table[i].Lr = log(table[i].Lr);
      table[i].Li = log(table[i].Li);
  }
  Ltot[0] = table[n-1].Lu;
  Ltot[1] = table[n-1].Lb;
  Ltot[2] = table[n-1].Lv;
  Ltot[3] = table[n-1].Lr;
  Ltot[4] = table[n-1].Li;

  PRC(Ltot[0]);PRC(Ltot[1]);PRC(Ltot[2]);PRC(Ltot[3]);PRC(Ltot[4]);

//  for(i=0; i<n; i++) {
//      Ltot[0] += table[i].Lu;
//      Ltot[1] += table[i].Lb;
//      Ltot[2] += table[i].Lv;
//      Ltot[3] += table[i].Lr;
//      Ltot[4] += table[i].Li;
//  }

    // now the array table is sorted in radius.
    // we can do all kind of beautiful things with it.
}


local real predifinced_lagrangian_fraction(int i) {

  real bin_radius;
  switch(i) {
  case 0: bin_radius = 0.005;
    break;
  case 1: bin_radius = 0.01;
    break;
  case 2: bin_radius = 0.02;
    break;
  case 3: bin_radius = 0.05;
    break;
  case 4: bin_radius = 0.01;
    break;
  case 5: bin_radius = 0.25;
    break;
  case 6: bin_radius = 0.50;
    break;
  case 7: bin_radius = 0.75;
    break;
  default: bin_radius = 0.90;
  }   

  return bin_radius;
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

local void get_lagrangian_ubvri_radii(dyn *b, int nzones, bin_option binning,
				      real *rLu, real *rLb, real *rLv, 
				      real *rLr, real *rLi) {

  // return if unresonable input is given.
  if (nzones < 2)
    return;

  int n = b->n_leaves();

  star_ubvri_ptr table = new star_table[n];
  if (table == NULL) {
    cerr << "sort_stellar_ubvri: "
	 << "not enough memory left for table\n";
    exit(0);
  }

  real Ltot[] = {0, 0, 0, 0, 0};
  sort_stellar_ubvri(b, nzones, n, table, Ltot);

  // real L_percent[5][nzones];	// g++ extension!

  if (nzones > 100) {
    cerr << "get_lagrangian_ubvri_radii: too many zones!" << endl;
    return;
  }

  real L_percent[5][100];	// ANSI

    // determine the lagrangian radii.
  for (int k=0; k<nzones; k++) {
      for (int l=0; l<5; l++) {
	  L_percent[l][k] = lagrangian_radius(k, nzones, binning)*Ltot[l];
      }
  }

  real Lcum[] = {0, 0, 0, 0, 0};
  int ipercent[] = {0, 0, 0, 0, 0};

  for (int k=0; k<nzones; k++) {
      for (int l=0; l<5; l++) {

	  while (Lcum[l] < L_percent[l][k]) {
	      switch(l) {
		  case 0: Lcum[l] = table[ipercent[l]].Lu;
			    break;
		  case 1: Lcum[l] = table[ipercent[l]].Lb;
			    break;
		  case 2: Lcum[l] = table[ipercent[l]].Lv;
			    break;
		  case 3: Lcum[l] = table[ipercent[l]].Lr;
			    break;
		  case 4: Lcum[l] = table[ipercent[l]].Li;
			    break;
			};
			    
	      ipercent[l]++;
	  }
      }

      rLu[k] = table[ipercent[0]-1].sort_param;
      rLb[k] = table[ipercent[1]-1].sort_param;
      rLv[k] = table[ipercent[2]-1].sort_param;
      rLr[k] = table[ipercent[3]-1].sort_param;
      rLi[k] = table[ipercent[4]-1].sort_param;
  }

  int p = cerr.precision(4);
  cerr << "    Total luminosities: " << exp(table[n-1].Lu) << " "
			             << exp(table[n-1].Lb) << " "
			             << exp(table[n-1].Lv) << " "
				     << exp(table[n-1].Lr) << " "
				     << exp(table[n-1].Li) 
       << " [Lsun]\n" << endl;
  cerr.precision(p);

  delete []table;

}    

local void print_lagrangian_ubvri_radii(dyn *b, int nzones, 
					bin_option binning) {

  cerr << "  (currently NON) Projected ubvri Lagrangian radii" << endl;
  cerr << "    [%]      u        b<        v        r        i"
       << endl;
  cerr << " Temporarily removed"<<endl;

#if 0
  real *rLu = new real[nzones];
  real *rLb = new real[nzones];
  real *rLv = new real[nzones];
  real *rLr = new real[nzones];
  real *rLi = new real[nzones];

  get_lagrangian_ubvri_radii(b, nzones, binning, rLu, rLb, rLv, rLr, rLi);

  cerr << "  (currently NON) Projected ubvri Lagrangian radii" << endl;
  cerr << "    [%]      u        b        v        r        i"
       << endl;
  
  int p = cerr.precision(4);
  int ncount, nstar = 0;
  for (int k=0; k<nzones; k++) 
      cerr << "     " << lagrangian_radius(k, nzones, binning) << "    " 
	   << rLu[k] << "   "<< rLb[k] << "   " << rLv[k] 
	             << "   " << rLr[k] << "   "<< rLi[k] << endl;
  cerr.precision(p);

  delete []rLu;
  delete []rLb;
  delete []rLv;
  delete []rLr;
  delete []rLi;
#endif
}

local void print_mass_over_light(dyn* b) {

    if (!b->get_use_sstar())
	return;

    // Use lagr_pos as center; do nothing if unavailable.

    if (find_qmatch(b->get_dyn_story(), "n_lagr")) {

	// Make a list of previously computed lagrangian radii.

	int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	real *r_lagr = new real[n_lagr];
	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	vector lagr_pos = 0;

	if (find_qmatch(b->get_dyn_story(), "lagr_pos"))
	    lagr_pos = getvq(b->get_dyn_story(), "lagr_pos");
	else
	    warning("print_mass_over_light: lagr_pos not found");

	int  *total_star = new int[n_lagr+1];
	real *total_mass = new real[n_lagr+1];
	real *total_lumi = new real[n_lagr+1];
	real *total_Ngt1Msun = new real[n_lagr+1];
	real *total_Ngt10Msun = new real[n_lagr+1];

	int i;
	for (i = 0; i <= n_lagr; i++)  {
	    total_star[i] = 0;
	    total_Ngt1Msun[i] = total_Ngt10Msun[i]
			      = total_mass[i] = total_lumi[i] = 0;
	}

	// Note that r_lagr omits the 0% and 100% radii.
    
	// Initialization:

	for (i = 0; i <= n_lagr; i++)
	    if (i < n_lagr) r_lagr[i] *= r_lagr[i];

	int N_total = 0, N_Mgt1Msun = 0, N_Mgt10Msun = 0;
	real L_star = 0, M_star = 0;
	bool Mgt1Msun, Mgt10Msun; 

	for_all_leaves(dyn, b, bi) {

	    // Count single stars and binaries separately.

	    Mgt1Msun = Mgt10Msun = false; 

	    if (bi->get_oldest_daughter() == NULL) 
		if (has_sstar(bi)) {

		    L_star = ((single_star*)bi->get_starbase())
					      ->get_luminosity();
		    M_star = bi->get_starbase()->get_total_mass();
		    if ( ((single_star*)bi->get_starbase())->remnant() )
			L_star = 0;
		    if(M_star>1) {
			Mgt1Msun = true;
			if(M_star>10)
			    Mgt10Msun = true;
		    }
		}
		else if (bi->get_star_story()!=NULL) {
	  
		    L_star = getrq(bi->get_star_story(), "L_eff"); 
		    M_star = getrq(bi->get_star_story(), "M_env")
				+ getrq(bi->get_star_story(), "M_core"); 
		    if(M_star>1) {
			Mgt1Msun = true;
			if(M_star>10)
			    Mgt10Msun = true;
		    }
		}

	    // Find which zone we are in.

	    i = which_zone(bi, lagr_pos, n_lagr, r_lagr);

	    // Update statistics.

	    N_total++;

	    total_mass[i] += M_star;
	    total_lumi[i] += L_star;
	    total_star[i]++;

	    if(Mgt1Msun) {
		total_Ngt1Msun[i]++;
		N_Mgt1Msun++;
	    }
	    if(Mgt10Msun) {
		total_Ngt10Msun[i]++;
		N_Mgt10Msun++;
	    }
	}

	cerr << "      Name   M(tot)     L(tot)  "
	     << "\tN(all)\t N(M>1)\t N(M>10Msun)\n";

	if (N_total > 0) {

	    for (i=0; i<=n_lagr; i++)
		if (total_lumi[i]>0 && total_star[i]>0)
		    cerr << "    zone " << i << "   "
			 << total_mass[i] <<"  \t"
			 << total_lumi[i] <<"  \t"
			 << total_star[i] <<"  \t"
			 << total_Ngt1Msun[i] <<"  \t"
			 << total_Ngt10Msun[i] << endl;
		else
		    cerr << "    zone " << i << " "
			 << " no mass   "
			 << " no luminosity    \tno stars"
			 << endl;

	} else
	    cerr << "            (none)\n";
  
	delete []r_lagr;
	delete []total_mass;
	delete []total_lumi;
	delete []total_star;
	delete []total_Ngt1Msun;
	delete []total_Ngt10Msun;
    }
}

local void print_mass_spectrum(dyn* b, int nzones, bool verbose) {

    if (!b->get_use_sstar() || nzones < 2)
	return;

    nm_bin *mass_bins = new nm_bin[nzones+1];
    if (mass_bins == NULL) {
	cerr << "not enough memory left for print_mass_spectrum\n";
	return;
    }

    real m_min = b->get_starbase()->conv_m_star_to_dyn(0.1);   // Msun
    real m_max = b->get_starbase()->conv_m_star_to_dyn(100);   // Msun
    
    real log_m_min = log10(m_min);
    mass_bins[1].mass_limit = log_m_min;
    real mass_int = (log10(m_max)-mass_bins[1].mass_limit)/(nzones-2);

    int i;
    for (i=0; i<nzones; i++) {
	mass_bins[i].mass_limit = pow(10, log_m_min + (i-1)*mass_int);
	mass_bins[i].no_of_stars = 0;
	mass_bins[i].total_mass = 0;
    }
    mass_bins[0].mass_limit = mass_bins[nzones].mass_limit = 0;
    
    int N_total = 0;
    real mi, m_total=0;
    for_all_leaves(dyn, b, bi) {
	mi = bi->get_mass();
	m_total += mi;
	N_total++;
	if (mi > m_min && mi <= m_max) {

	    for (i=2; i<nzones-1; i++)
	  
		if (mi > mass_bins[i-1].mass_limit   &&
		    mi <= mass_bins[i].mass_limit) {

		    mass_bins[i-1].no_of_stars++;
		    mass_bins[i-1].total_mass  += mi;
		}
	}
	else if (mi<=m_min) {

	    mass_bins[0].no_of_stars++;
	    mass_bins[0].total_mass  += mi;
	}
	else {

	    mass_bins[nzones].no_of_stars++;
	    mass_bins[nzones].total_mass  += mi;
	}
    }
    
    real scaling = 1;
    if (b->get_use_sstar()) {
	scaling = b->get_starbase()->conv_m_star_to_dyn(1);
	m_min = b->get_starbase()->conv_m_dyn_to_star(m_min);
	m_max = b->get_starbase()->conv_m_dyn_to_star(m_max);
    }
    int p = cerr.precision(LOW_PRECISION);

    if (N_total > 0) {

	cerr << endl;
	if (verbose)
	    cerr << "  Mass function (dlog_m = " << mass_int
		 << "  min = " << m_min
		 << "  max = " << m_max << " [Msun]): "
		 << endl;

	cerr << "            " << mass_bins[0].no_of_stars
	     << "  <";
	for (i=1; i<=nzones-2; i++)
	    cerr << "   " << mass_bins[i].no_of_stars;

	cerr << "   > " << mass_bins[nzones-1].no_of_stars
	     << endl;
    }
    else
	cerr << "            (none)\n";

    cerr.precision(p);    
    delete []mass_bins;
}


local void print_luminosity_function(dyn* b, int nzones, bool verbose) {

    if (!b->get_use_sstar() || nzones < 2)
	return;

    nm_bin *lstar_bins = new nm_bin[nzones+1];
    if (lstar_bins == NULL) {
	cerr << "not enough memory left for print_luminosity_function\n";
	return;
    }

    real l_min  = 0.01;
    real l_max  = 1.e+6;
    real log_l_min  = log10(l_min);
    real log_l_max  = log10(l_max);
    
    lstar_bins[1].mass_limit = log_l_min;
    real lstar_int = (log_l_max-log_l_min)/(nzones-2);

    int i;
    for (i=0; i<nzones; i++) {
	lstar_bins[i].mass_limit = pow(10, log_l_min + (i-1)*lstar_int);
	lstar_bins[i].no_of_stars = 0;
	lstar_bins[i].total_mass = 0;
    }
    lstar_bins[0].mass_limit = lstar_bins[nzones].mass_limit = 0;
    
    int N_total = 0;
    real li, l_total=0;
    for_all_leaves(dyn, b, bi) {

	li = ((star*)bi->get_starbase())->get_luminosity();
	l_total += li;
	N_total++;
	if (li > l_min && li <= l_max) {

	    for (i=2; i<nzones-1; i++)

		if (li > lstar_bins[i-1].mass_limit   &&
		    li <= lstar_bins[i].mass_limit) {
		    lstar_bins[i-1].no_of_stars++;
		    lstar_bins[i-1].total_mass  += li;
		}
	}
	else if (li <= l_min) {

	    lstar_bins[0].no_of_stars++;
	    lstar_bins[0].total_mass  += li;
	}
	else {

	    lstar_bins[nzones].no_of_stars++;
	    lstar_bins[nzones].total_mass  += li;
	}
    }
    
    if (N_total > 0) {

	cerr << endl;
	if (verbose)
	    cerr << "  Luminosity function (dlog_L = " << lstar_int
		 << "  min = 10^" << log_l_min
		 << "  max = 10^" << log_l_max << " [Lsun]): "
		 << endl;

	cerr << "            " << lstar_bins[0].no_of_stars
	     << "  <";
	for (i=1; i<=nzones-2; i++)
	    cerr << "   " << lstar_bins[i].no_of_stars;

	cerr << "   > " << lstar_bins[nzones].no_of_stars
	     << endl;
	  
    }
    else
	cerr << "            (none)\n";
    
    delete [] lstar_bins;
}

local void print_massive_star(dyn *bi, vector center_pos,
			      vector center_vel, bool verbose) {

    real mass  = bi->get_starbase()->get_total_mass();
    real r_com = abs(bi->get_pos() - center_pos);
    real v_com = abs(bi->get_vel() - center_vel);

    cerr << "    star#  " << bi->format_label()
	 << "      " << mass << "  ("; put_state(make_star_state(bi), cerr);
    cerr << ") \t  " << r_com << "   " << v_com << endl;
}
    
#if 0
local int print_massive_binary(dyn *bi,
			       vector center_pos, vector center_vel,
			       real mass_limit, real number_limit,
			       bool verbose) {


    real m_tot  = bi->get_starbase()->conv_m_dyn_to_star(bi->get_mass());
    real mass1  = bi->get_oldest_daughter()
		    ->get_starbase()->get_total_mass();
    int index1  = bi->get_oldest_daughter()->get_index();
    real mass2  = bi->get_oldest_daughter()->get_younger_sister()
		    ->get_starbase()->get_total_mass();
    int index2  = bi->get_oldest_daughter()->get_younger_sister()->get_index();

    if ((bi->get_index() > 0 && bi->get_index() <= number_limit) ||
	(mass1 >= mass_limit || (index1 >0 && index1 <= number_limit)) ||
	(mass2 >= mass_limit || (index2 >0 && index2 <= number_limit))) {
    
	real r_com  = abs(bi->get_pos() - center_pos);
	real v_com  = abs(bi->get_vel() - center_vel);

	cerr << "            " << bi->format_label()
	     << "      " << m_tot 
	     << "  (Multiple) \t  " << r_com << "   " << v_com << endl;

	cerr << "                     " << mass1 << "  (";
	put_state(make_star_state(bi->get_oldest_daughter()), cerr);
	cerr << ") " << endl;
	cerr << "                     " << mass2 << "  (";
	put_state(make_star_state(bi->get_younger_sister()), cerr);
	cerr << ") " << endl;
    }

    return 1;
}
#endif

local bool contains_massive_star(dyn *b,
				 real mass_limit,
				 real number_limit) {

    for_all_leaves(dyn, b, bi) 
	if ((bi->get_index()>0 && bi->get_index() <= number_limit) ||
	    bi->get_starbase()->get_total_mass() >= mass_limit) 
	    return true;

    return false;
}

local int print_massive_binary_recursive(dyn *b,
					 vector center_pos, vector center_vel,
					 real mass_limit, real number_limit,
					 bool verbose) {

    int nb = 0;
    if (b->get_oldest_daughter() &&
	contains_massive_star(b, mass_limit, number_limit)) {

	nb ++;

	real r_com  = abs(b->get_pos() - center_pos);
	real v_com  = abs(b->get_vel() - center_vel);

	real m_tot = b->get_starbase()->conv_m_dyn_to_star(b->get_mass());
	cerr << "            " << b->format_label()
	     << "      " << m_tot 
	     << "  (Multiple) \t  " << r_com << "   " << v_com << endl;

	for_all_daughters(dyn, b, bb)
	    if (bb->n_leaves() >= 2) 
		nb += print_massive_binary_recursive(bb,
						     center_pos, center_vel,
						     mass_limit, number_limit,
						     verbose); 
	    else
		print_massive_star(bb, center_pos, center_vel, verbose);

	// nb += print_massive_binary(b, center_pos, center_vel,
	//			      mass_limit, number_limit,  verbose);
    }
    return nb;
}

local void print_massive_star_header(bool cod) {

    cerr << endl;
    cerr << "  Massive stars:" << endl;
    if (cod)
	cerr << "            Name    M(sun)   (type)   \t r_cod    v_cod"
	     << endl;
    else
	cerr << "            Name    M(sun)   (type)   \t r_com    v_com"
	     << endl;
}

typedef  struct
{
    real mass;
    dyn *b;
} md_pair, *md_pair_ptr;

local int compare_mass(const void * pi, const void * pj)  // decreasing mass
{
    if (((md_pair_ptr) pi)->mass > ((md_pair_ptr) pj)->mass)
        return -1;
    else if (((md_pair_ptr)pi)->mass < ((md_pair_ptr)pj)->mass)
        return +1;
    else
        return 0;
}


local void print_most_massive_stars(dyn *b,
				    real mass_limit,
				    real number_limit,
				    vector center_pos,
				    bool verbose) {

    bool cod = false;
    vector center_vel = 0;

    if (abs(center_pos) == 0)
	cod = (get_std_center(b, center_pos, center_vel) == 1);

    // Want to print out stars more massive than mass_limit, but no more
    // than number_limit.  New code from Steve (8/01).

    int nl = 0;
    for_all_daughters(dyn, b, bi)
	if (bi->is_leaf()) nl++;
    
    md_pair_ptr md_list = new md_pair[nl];

    nl = 0;
    for_all_daughters(dyn, b, bi)
	if (bi->is_leaf()) {
	    md_list[nl].mass = bi->get_starbase()->get_total_mass();
	    md_list[nl].b = bi;
	    nl++;
	}

    // Sort the list, in order of decreasing mass.

    qsort((void *)md_list, (size_t)nl, sizeof(md_pair), compare_mass);

    int ns = 0;
//    if (nl < rint(number_limit)) nl = (int)rint(number_limit);

    for (ns = 0; ns < nl; ns++)
	if (md_list[ns].mass < mass_limit) break;

    if (ns < rint(number_limit)) ns = Starlab::min(nl, (int)rint(number_limit));

    print_massive_star_header(cod);

    // Stars:

    for (int i = 0; i < ns; i++)
	print_massive_star(md_list[i].b, center_pos, center_vel, verbose);

    delete [] md_list;

    // Binaries:

    int nb = 0;
    for_all_daughters(dyn, b, bi) {
	if (!bi->is_leaf())                   
	    nb += print_massive_binary_recursive(bi, center_pos, center_vel,
						 mass_limit, number_limit,
						 verbose); 
    }
  
    if (verbose && ns+nb==0) 
	cerr << "            (none)\n";
}

// Print stellar timescale and turn-off mass.

void print_sstar_time_scales(dyn *b) {
  
    real time = b->get_system_time();
    cerr << "     = " << b->get_starbase()->conv_t_dyn_to_star(time)
	 << " Myr (Turn_off_mass = "
	 << turn_off_mass(b->get_starbase()->conv_t_dyn_to_star(time))
	 << ")" << endl;

    b->get_starbase()->print_stellar_evolution_scaling(cerr);
}

void sstar_stats(dyn* b, bool mass_spectrum, vector center,
		 bool verbose) { 

  if (b->get_use_sstar()) {

      int axis = 0;   // x-axis
      if (verbose) cerr << "\n  Projected luminosity profile: ";
      real r1_r9 = compute_projected_luminosity_radii(b, axis, true, 10);

      if (verbose) cerr << "\n  Projected King profile:"; 
      print_fitted_king_model(r1_r9, proj_r10_r90_light);

      if (mass_spectrum) {
      
	  print_luminosity_function(b, 12, verbose);
	  print_mass_spectrum(b, 12, verbose);
	  print_most_massive_stars(b, OUTPUT_MASS_LIMIT, OUTPUT_NUMBER_LIMIT, 
				   center, verbose);
      }
    
      if (verbose)
	  cerr << endl
	       << "  M, L and massive stars by Lagrangian zone (all stars):"
	       << endl;

      print_mass_over_light(b);    

      if (verbose)
	  cerr << endl
	       << "  stellar ubvr and i Lagrangian zones (all stars):"
	       << endl;

      print_lagrangian_ubvri_radii(b, MAX_PREDEF_RADII, predefined);

      if (verbose)
	  cerr << "\n  Stellar content:\n";
      print_stellar_content(b);

      if (verbose)
	  cerr << "\n  Special star content:\n";
      print_special_stars(b);

      if (verbose)
	  cerr << "\n  Compact stars:\n";
      print_compact_stars(b);

  }
  else 
      cerr << "\n  No stellar evolution\n";
      }

#else

main(int argc, char **argv)
{
    bool binaries = true, verbose = true, out = false,
         n_sq = true, calc_e = true;
    int which_lagr = 0;
    real T_start = 0;

    extern char *poptarg;
    int c;
    char* param_string = "bneost";	// Note: "v" removed because only the
					// "verbose" option currently works.

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'b': binaries = !binaries;
		      break;
	    case 'e': calc_e = !calc_e;
		      break;
	    case 'n': n_sq = !n_sq;
		      break;
	    case 'o': out = !out;
		      break;
	    case 's': which_lagr = 2;
		      break;
	    case 't': which_lagr = 1;
		      break;
	    case 'v': verbose = !verbose;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Loop over input until no more data remain.

    dyn *b;
    int i = 0;
    vector zero = 0;
    bool mass_spectrum = true;
    while (b = get_dyn(cin)) {

        // NOTE:  get_xxx() reads NaN in as legal garbage...

      if(find_qmatch(b->get_oldest_daughter()
		 ->get_starbase()->get_star_story(), "Type")) {
	
	addstar(b,                            // Note that T_start and
		b->get_system_time(),          // Main_Sequence are
		Main_Sequence,                 // defaults. They are
		true);                         // ignored if a star
	b->set_use_sstar(true);                // is there.

      if (i++ > 0) cerr << endl;
      sstar_stats(b, mass_spectrum, zero,  verbose);

      }
      
      if (out) put_node(cout, *b);
      rmtree(b);
    }
}

#endif
