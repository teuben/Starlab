//  coagulation.C: SPZ:May 1998
//                 the dynamical evolution of the mass spectrum is computed 
//                 due to collisional evolution in a semi-analytic
//                 fashion.
//                 In first instance only a single step is performed
//                 with time T 

#include "dyn.h"

#ifdef TOOLBOX

// cerr << "    [R]: " << Rsun_pc / r_conv_star_to_dyn << " pc\n";

#define Rsun_pc 2.255e-8                 // R_sun/parsec = 6.960e+10/3.086e+18


#define ALL_i ni = n->get_oldest_daughter(); ni != NULL; \
					     ni = ni->get_younger_sister()
#define j_ABOVE_i nj = ni->get_younger_sister(); nj != NULL; \
					     nj = nj->get_younger_sister()


local real conv_v_dyn_to_star(real v, real rf, real tf) {

//              Internal velocity is km/s
//              Internal size is solar raii
//              Internal time is Myr.
//              km/s to Rsun/Myr
//      real to_Rsun_Myr = cnsts.physics(km_per_s) * cnsts.physics(Myear)
//	               / cnsts.parameters(solar_radius);

      real myear = 3.15e+13;
      real km_s = 1.0e+5;
      real R_sun = 6.960e+10;
      real to_Rsun_Myr = km_s * myear / R_sun;
      real to_dyn      = rf/tf;
      
      return v/(to_Rsun_Myr * to_dyn);
   }


typedef  struct
{
    real  radius;
    real  mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)
{
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return(1);
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return(-1);
    else
        return(0);
}

//-----------------------------------------------------------------------------
//  compute_general_mass_radii  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

static real nonlin_masses[9] = {0.005, 0.01, 0.02, 0.05, 0.1,
				0.25, 0.5, 0.75, 0.9};

local void compute_general_mass_radii(dyn * b, vector dc_pos, int nzones,
				      real r_lagr[],
				      real m_cutoff, real M_cutoff,
				      bool nonlin = false, 
				      boolfn bf=NULL)
{
    if (nzones < 2) return;
    if (nonlin && nzones != 10) return;		// Special case

    // Note: nzones specifies the number of radial zones to consider.
    //	     However, we only store radii corresponding to the nzones-1 
    //	     interior zone boundaries.

    real* mass_percent = new real[nzones-1];
    if (mass_percent == NULL) {
	cerr << "compute_general_mass_radii: "
	     << "not enough memory left for mass_percent\n";
	return;
    }

    int n = 0;
    if (bf == NULL) {
      for_all_daughters(dyn, b, bb)
	if (bb->get_mass()>=m_cutoff && bb->get_mass()<=M_cutoff)
	  n++;
    }
    else {
	for_all_daughters(dyn, b, bb)
	  if (bb->get_mass()>=m_cutoff && bb->get_mass()<=M_cutoff)
	    if ((*bf)(bb))
	      n++;
    }

    rm_pair_ptr rm_table = new rm_pair[n];

    if (rm_table == NULL) {
	cerr << "compute_general_mass_radii: "
	     << "not enough memory left for rm_table\n" << flush;
	return;
    }

    // Use the geometric center if the density center is unknown.

    //vector dc_pos = 0;
    //if (find_qmatch(b->get_dyn_story(), "density_center_pos")) 
    //dc_pos = getvq(b->get_dyn_story(), "density_center_pos");

    // Set up an array of (radius, mass) pairs.  Also find the total
    // mass of all nodes under consideration.

    real total_mass = 0;

    int i = 0;
    for_all_daughters(dyn, b, bi) {
      if (bf == NULL || (*bf)(bi)) {
	if (bi->get_mass()>=m_cutoff && bi->get_mass()<=M_cutoff) {
	  total_mass += bi->get_mass();
	  rm_table[i].radius = abs(bi->get_pos() - dc_pos);
	  rm_table[i].mass = bi->get_mass();
	  i++;
//	    PRL(i);
	}
      }
    }	

    // Sort the array by radius.

    qsort((void *)rm_table, (size_t)n, sizeof(rm_pair), compare_radii);

    // Determine Lagrangian radii.

    //    cerr << "Determine Lagrangian radii" << endl;
    int k;
    for (k = 0; k < nzones-1; k++) {
        if (!nonlin) 
	    mass_percent[k] = ((k + 1) / (real)nzones) * total_mass;
	else
	    mass_percent[k] = nonlin_masses[k] * total_mass;
    }

    //    real *rlagr = new real[nzones-1];
    real cumulative_mass = 0.0;
    i = 0;


    for (k = 0; k < nzones-1; k++) {

      while (cumulative_mass < mass_percent[k]) {
	cumulative_mass += rm_table[i++].mass;
      }

      r_lagr[k] = rm_table[i-1].radius;
    }


    // Place the data in the root dyn story.

#if 0
    if (bf == NULL)
	putiq(b->get_dyn_story(), "boolfn", 0);
    else
	putiq(b->get_dyn_story(), "boolfn", 1);

    putiq(b->get_dyn_story(), "n_nodes", n);

    putiq(b->get_dyn_story(), "n_lagr", nzones-1);
    putra(b->get_dyn_story(), "r_lagr", rlagr, nzones-1);
#endif

    delete mass_percent;
    delete rm_table;
//    delete rlagr;
    
    for (k = 0; k < nzones-1; k ++) 
      cerr << " " << b->get_starbase()->conv_r_dyn_to_star(r_lagr[k])
	             * Rsun_pc;
    cerr << " [pc]" << endl;
    
}

// Convenient synonyms:

local void  compute_mass_radii_quartiles(dyn * b, vector dc_pos,
					 real m_cutoff, real M_cutoff)
{
  int nl = 4;
  real *r_lagr = new real[nl-1];
  for (int i=0; i<nl-1; i++)
    r_lagr[i] = 0;

    compute_general_mass_radii(b, dc_pos, nl, r_lagr, m_cutoff, M_cutoff);
}

local void  compute_mass_radii_percentiles(dyn * b, vector dc_pos,
					   real m_cutoff, real M_cutoff)
{
  int nl = 10;
  real *r_lagr = new real[nl-1];
  for (int i=0; i<nl-1; i++)
    r_lagr[i] = 0;
  
    compute_general_mass_radii(b, dc_pos, nl, r_lagr, m_cutoff, M_cutoff);
}

local bool single_fn(dyn * b)
{

  //  int test = getiq(b->get_log_story(), "mass_doubled");
  //  PRL(test);
  if (getiq(b->get_log_story(), "mass_doubled") == 0 )
    return true;
  else
    return false;
}

local bool double_fn(dyn * b)
{
  if (getiq(b->get_log_story(), "mass_doubled") == 1 )
    return true;
  else
    return false;
}  

local void print_lagrangian_radii(dyn* b, vector dc_pos, int which_lagr,
				  real m_cutoff, real M_cutoff,
				  int which = 0)
{
    bool nonlin = false;

    int nl=0, indent=0;
    if (which_lagr == 0) {
	nl = 4;
	indent = 15;
	PRI(4); cerr << "quartiles: ";
    } else if (which_lagr == 1) {
	nl = 10;
	indent = 20;
	PRI(4); cerr << "10-percentiles: ";
    } else {
	nl = 10;
	indent = 26;
	PRI(4); cerr << "selected percentiles: ";
	nonlin = true;
    }
    
    //cerr << "before testnode call" << endl;
    //    if (which != 2 || 3) compute_general_mass_radii(b, nl, nonlin, testnode);

    real* r_lagr = new real[nl-1];
    if (r_lagr==NULL) {
      cerr << "print_lagrangian_radii not enough memory to allocate "
	   << "r_lagr[" << nl-1<< "]" << endl;
      return;
    }
    for (int i=0; i<nl-1; i++)
      r_lagr[i] = 0;
    
    if (which == 0) compute_general_mass_radii(b, dc_pos, nl, r_lagr,
					       m_cutoff, M_cutoff, nonlin); 

    //    cerr << "before single_fn call" << endl;
    //    PRL(which);
    if (which == 2) compute_general_mass_radii(b, dc_pos, nl, r_lagr,
					       m_cutoff, M_cutoff, nonlin, single_fn);

    //     cerr << "before double_fn call" << endl;
    //    PRL(which);
    if (which == 3) compute_general_mass_radii(b, dc_pos, nl, r_lagr,
					       m_cutoff, M_cutoff, nonlin, double_fn);

    //    PRC(dc_pos);PRC(m_cutoff);PRC(M_cutoff);PRL(which);

    //for (int k = 0; k < nl; k ++) 
    // cerr << " " << r_lagr[k];
    //cerr <<endl;

#if 0
    if (find_qmatch(b->get_dyn_story(), "n_lagr")) {
        int n_lagr = getiq(b->get_dyn_story(), "n_lagr");
	getra(b->get_dyn_story(), "r_lagr", r_lagr, n_lagr);

	for (int k = 0; k < n_lagr; k += 5) {
	  if (k > 0) {
	    cerr << endl;
	    for (int kk = 0; kk < indent; kk++) cerr << " ";
	  }
	  for (int i = k; i < k+5 && i < n_lagr; i++)
	    cerr << " " << r_lagr[i];
	}
	cerr << endl;
    }
#endif
    
    delete r_lagr;
}

/*-----------------------------------------------------------------------------
 *  compute_density  --  Get the density for all particles.
 *               note:
 *                    the neighbor_dist_sq[q] array will contain the distances
 *                    to the q-th nearest neighbor (in the form of the squares 
 *                    thereof); therefore  neighbor_dist_sq[0] = 0.0  since
 *                    this indicates the distance of a particle to itself.
 *                    Including this trivial case simplifies the algorithm
 *                    below by making it more homogeneous.
 *-----------------------------------------------------------------------------
 *              Litt.:
 *                    Stefano Casertano and Piet Hut: Astroph.J. 298,80 (1985).
 *                    To use their recipe, k >= 2 is required.
 *-----------------------------------------------------------------------------
 */
local void  compute_density(dyn * b,   	/* pointer to an Nbody system */
			    int k,	/* use kth nearest neighbor */
			    int n,
			    real local_density[],
			    real m_cutoff,
			    real M_cutoff)
{
    int  q, r, is = 0;
    real *neighbor_dist_sq;
    real  volume;
    real  delr_sq;
    
    if (k >= n)                           /* no k-th nearest neighbor exists */
        cerr << "compute_density: k = " << k << " >= n = " << n << endl;
    if (k <= 1)
        cerr << "compute_density: k = " << k << " <= 1" << endl;

    neighbor_dist_sq = new real[k+1];

    //  for_all_leaves(dyn, b, d)		// The density is now defined
    for_all_daughters(dyn, b, d)	// ONLY for top-level nodes.
      if (d->get_mass()>=m_cutoff && d->get_mass()<=M_cutoff)
        {
	neighbor_dist_sq[0] = 0.0;
	for (q = 1; q <= k; q++)
	    neighbor_dist_sq[q] = VERY_LARGE_NUMBER;

	//      for_all_leaves(dyn, b, dd)
	for_all_daughters(dyn, b, dd)
	  if (dd->get_mass()>=m_cutoff && dd->get_mass()<=M_cutoff)
	    {
	    if (d == dd)
	        continue;

	    delr_sq = square(something_relative_to_root(d, &dyn::get_pos)
			     - something_relative_to_root(dd, &dyn::get_pos));

	    if (delr_sq < neighbor_dist_sq[k])
	        for (q = k-1; q >= 0; q--)
	            {
		    if (delr_sq > neighbor_dist_sq[q])
		        {
			for (r = k; r > q+1; r--)
			    neighbor_dist_sq[r] = neighbor_dist_sq[r-1];
			neighbor_dist_sq[q+1] = delr_sq;
		        break;
			}
		    }
	    }
	    
	volume = (4.0/3.0) * PI * pow(neighbor_dist_sq[k], 1.5);

        real density =  (k - 1) / volume;           /* Ap. J. 298, 80 (1985) */
	story * s = d->get_dyn_story();
	//        putiq(s, "density_k_level", k);
	//        putrq(s, "density", density);
	local_density[is] = density;
	is++;
	}

    delete neighbor_dist_sq;
}



/*-----------------------------------------------------------------------------
 *  compute_mean_cod -- Returns the position and velocity of the density
 *		        center of an N-body system.
 *-----------------------------------------------------------------------------
 */
local void compute_mean_cod(dyn *b,
 			    real local_density[],
			    vector& pos, vector& vel, real m_cutoff,
			    real M_cutoff)
{
    real total_weight = 0;
    pos = 0;
    vel = 0;
    int is=0;
    
    for_all_daughters(dyn, b, d) 
      if (d->get_mass()>=m_cutoff && d->get_mass()<=M_cutoff) {
	    
	//	real density = getrq(d->get_dyn_story(), "density");
	real density = local_density[is++];
	real dens2 = density * density;        // weight factor

        total_weight += dens2;
	pos += dens2 * d->get_pos();
	vel += dens2 * d->get_vel();
    }	

    pos /= total_weight;
    vel /= total_weight;
}

/*===========================================================================*/

local void compute_core_parameters(dyn* b, int k, vector& center,
				   real& rcore, int& ncore, real& mcore,
				   real & vcore, int n,
				   real local_density[], 
				   real m_cutoff, real M_cutoff)
{
    vector vel;
    int is=0;

    real rcore2 = 0;
    if (rcore<=0) {
      compute_density(b, k, n, local_density, m_cutoff, M_cutoff);
      compute_mean_cod(b, local_density, center, vel, m_cutoff, M_cutoff);

      real total_weight = 0;
      real sum = 0;
      for_all_daughters(dyn, b, bi) 
	if (bi->get_mass()>=m_cutoff && bi->get_mass()<=M_cutoff) {
	  real density = local_density[is++];
	  real dens2 = density * density;        // weight factor
	  //real density = getrq(bi->get_dyn_story(), "density");
	  //real dens2 = density * density;              // weight factor

	  total_weight += dens2;
	  sum += dens2 * square(bi->get_pos() - center);
	}

      rcore2 = sum/total_weight;
      rcore = sqrt(rcore2);
    }
    else
      rcore2 = square(rcore);

    ncore = 0;
    mcore = 0;
    vcore = 0;

    for_all_daughters(dyn, b, bj)
      if (bj->get_mass()>=m_cutoff && bj->get_mass()<=M_cutoff) 
	if (square(bj->get_pos() - center) <= rcore2) {
	    ncore++;
	    mcore += bj->get_mass();
	    vcore += bj->get_vel()*bj->get_vel();
	}
    mcore = mcore/(1.0*ncore);
    vcore = sqrt(vcore)/(1.0*ncore);
}

local void print_core_parameters(dyn* b, vector& density_center, real& rcore,
				 int & ncore, real & vcore,
				 int n, real local_density[],
				 real m_cutoff, real M_cutoff,
				 bool verbose)
  {
    real mcore;

    compute_core_parameters(b, 6, density_center, rcore, ncore, mcore,
			    vcore, n, local_density,
			    m_cutoff, M_cutoff);

    real r_density_center = sqrt(density_center*density_center);

    vcore = conv_v_dyn_to_star(vcore,
			       b->get_starbase()->conv_r_star_to_dyn(1),
			       b->get_starbase()->conv_t_star_to_dyn(1));
    if (verbose) {

      real M_lower = b->get_starbase()->conv_m_dyn_to_star(m_cutoff);
      real M_upper = b->get_starbase()->conv_m_dyn_to_star(M_cutoff);

      cerr << "n = " << n << ", M (lower/upper) = ( "
	                  << M_lower << ",  " << M_upper << " ) [Msun]\n";
      cerr << "     Core: (N, R, m, v) = ( " << ncore << ", "
	   << Rsun_pc * b->get_starbase()->conv_r_dyn_to_star(rcore)
	   << ", "  
	   << b->get_starbase()->conv_m_dyn_to_star(mcore) <<", "
	   << vcore
	   << " ) [solar/km/s]\n";
      cerr << "     R_to_density_center = "
           << Rsun_pc *
	      b->get_starbase()->conv_r_dyn_to_star(r_density_center)
	   << " [pc]\n";
      

      /*
      //place the data in the root dyn-story.
      putiq(b->get_dyn_story(), "n_cutoff", n);
      putrq(b->get_dyn_story(), "m_cutoff", m_cutoff);
      putrq(b->get_dyn_story(), "M_cutoff", M_cutoff);
      putvq(b->get_dyn_story(), "density_center_pos", density_center);
      putrq(b->get_dyn_story(), "core_radius", rcore);
      putiq(b->get_dyn_story(), "n_core", ncore);
      putrq(b->get_dyn_story(), "m_core", mcore);
      */
    }
    else {
              PRC(n); PRC(m_cutoff); PRL(M_cutoff);
      PRI(4); PRC(rcore); PRC(ncore); PRL(mcore);
      PRI(4); PRL(r_density_center);
    }
	
}

local real system_energy(dyn* b, real m_cutoff, real M_cutoff)
{
    real kin = 0;
    for_all_leaves(dyn, b, bj)
      if (bj->get_mass()>=m_cutoff && bj->get_mass()<=M_cutoff)
	kin += bj->get_mass()
	     * square(something_relative_to_root(bj, &dyn::get_vel));
    kin *= 0.5;

    real pot = 0.0;
    for_all_leaves(dyn, b, bi) 
      if (bi->get_mass()>=m_cutoff && bi->get_mass()<=M_cutoff) {
	real dpot = 0.0;
	for_all_leaves(dyn, b, bj) 
	  if (bj->get_mass()>=m_cutoff && bj->get_mass()<=M_cutoff) {
	    if (bj == bi) break;
	    vector dx = something_relative_to_root(bi, &dyn::get_pos)
			  - something_relative_to_root(bj, &dyn::get_pos);
	    dpot += bj->get_mass() / abs(dx);
	}
	pot -= bi->get_mass() * dpot;
    }

    return kin + pot;
}

local void get_energies(dyn * n, real eps2, 
		  real& kinetic_energy, real& potential_energy,
		  real m_cutoff, real M_cutoff)
{
    dyn  * ni, *nj;

    kinetic_energy = 0;
    for (ALL_i) 
      if (ni->get_mass()>=m_cutoff && ni->get_mass()<=M_cutoff) {
      
	vector v = ni->get_vel();
        kinetic_energy += ni->get_mass() * v * v;
    }
    kinetic_energy *= 0.5;

    potential_energy = 0;
    for (ALL_i) 
      if (ni->get_mass()>=m_cutoff && ni->get_mass()<=M_cutoff) {
        real dphi = 0;
	vector xi = ni->get_pos();
	for (j_ABOVE_i) 
	  if (nj->get_mass()>=m_cutoff && nj->get_mass()<=M_cutoff) {
	    vector xij = nj->get_pos() - xi;
	    dphi += nj->get_mass()/sqrt(xij*xij + eps2);
	}
	potential_energy -= ni->get_mass() * dphi;
    }
}


local void print_energies(dyn* b, real& kT, real m_cutoff, real M_cutoff, bool verbose)
{
    // Energies (top-level nodes):

    real kinetic_energy, potential_energy;
    get_energies(b, 0.0, kinetic_energy, potential_energy, m_cutoff, M_cutoff);

    real total_energy = kinetic_energy + potential_energy;
    real virial_ratio = -kinetic_energy / potential_energy;
    kT = -total_energy / (1.5*b->n_daughters());

    if (verbose) {
      PRI(4); cerr << "top-level nodes:\n";
      PRI(8); PRC(kinetic_energy); PRL(potential_energy);
      PRI(8); PRC(total_energy); PRC(kT); PRL(virial_ratio);

      PRI(4); cerr << "total energy (full system) = "
		   << system_energy(b, m_cutoff, M_cutoff) << endl;
    }
    else
      cerr <<"\t"<< kinetic_energy <<" "<< potential_energy
	   <<" "<< total_energy <<" "<< kT
	   <<" "<< virial_ratio <<" "<< system_energy(b, m_cutoff, M_cutoff) << endl;
}

local void print_relaxation_time(dyn* b, real m_cutoff, real M_cutoff, bool verbose)
{
    // Print the relaxation time

    real potential_energy, kinetic_energy;
    real r_virial, t_relax;
    real total_mass = 0;

    int n=0;
    for_all_leaves(dyn, b, bj) 
      if (bj->get_mass()>=m_cutoff && bj->get_mass()<=M_cutoff) {
	n++;
	total_mass += bj->get_mass();
    }

    get_energies(b, 0.0, kinetic_energy, potential_energy, m_cutoff, M_cutoff);

    r_virial = -0.5 * total_mass * total_mass / potential_energy;

    t_relax = 9.62e-2 * sqrt(r_virial * r_virial * r_virial / total_mass)
              * n / log10(0.4 * n);

    if (verbose) {
      PRI(4); cerr << "r_virial = "
		   << Rsun_pc *
		      b->get_starbase()->conv_r_dyn_to_star(r_virial) << endl;
      PRI(4); cerr << "t_relax = "
		   <<  b->get_starbase()->conv_t_dyn_to_star(t_relax)  << endl;
    }
    else
      cerr <<"\t"<< r_virial <<" "<< t_relax << endl;
}

local void get_composition(dyn* b,
			   real m_cutoff, real M_cutoff,
			   int& n, real& mmean,
			   real& rmean,
			   real& vdisp) {
  
  n=0;
  real T_eff, L_eff, R_eff;
  real m_total=0, v2_disp=0, r_total=0;
  
  for_all_daughters(dyn, b, bi)
    if (bi->get_mass()>=m_cutoff && bi->get_mass()<=M_cutoff) {
      n++;
      m_total += bi->get_mass();
      v2_disp += bi->get_vel()*bi->get_vel();
      //R_eff = getrq(bi->get_star_story(), "R_eff");
      //if (R_eff<=0) {

      T_eff = getrq(bi->get_star_story(), "T_eff");
      L_eff = getrq(bi->get_star_story(), "L_eff");
      R_eff = sqrt(Starlab::max(0.1, (1130.*L_eff)/pow(T_eff, 4)));
	    
      r_total += R_eff;
    }

  if (n>0) {
    mmean = b->get_starbase()->conv_m_dyn_to_star(m_total)/(1.0*n);
    rmean = r_total/(1.0*n);

    vdisp = sqrt(v2_disp)/(1.0*n);
    vdisp = conv_v_dyn_to_star(vdisp,
			       b->get_starbase()->conv_r_star_to_dyn(1),
			       b->get_starbase()->conv_t_star_to_dyn(1));
  }
}

local real encounter_rate(int im, real ni, real mi, real ri, real vi, 
			  real rci,
			  int jm, real nj, real mj, real rj, real vj, 
			  real rcj) 
{
	  
    real to_volume = 4*PI/3.;

    //PRC(im);PRC(ni);PRC(mi);PRC(ri);PRC(vi);PRL(rci);
    //PRC(jm);PRC(nj);PRC(mj);PRC(rj);PRC(vj);PRL(rcj);
  
    real n_rhoi = 1./(to_volume * pow(rci, 3));
    real n_rhoj = 1./(to_volume * pow(rcj, 3));
    real r_tot  = ri + rj;
    real m_tot  = mi + mj;
    real v_rel  = sqrt(pow(vi, 2) + pow(vj, 2));
    real rcore = Starlab::min(rci, rcj);

    real geometric = 6.31e-15 * v_rel * pow(r_tot, 2);
    real grav_foc = 3.61e-9 * m_tot * r_tot / v_rel;

    real rate = 0;
	
    if (im==jm) {

      if (ni>1)
	rate = 0.5 * ((ni-1) * n_rhoi)
                   *  (nj * n_rhoj)
                   * pow(rcore, 3)
                   * (grav_foc + geometric);
    }
    else
      rate = (ni * n_rhoi)
           * (nj * n_rhoj)
           * pow(rcore, 3)
	   * (grav_foc + geometric);
    
#if 0
    if(rate>0)
      cerr << jm <<" "<< im <<" "
	   << mj <<" "<< mi <<" "
	   << nj <<" "<< ni <<" "
	   << rate << endl;
#endif

    return rate;
}



//#else

main(int argc, char **argv)
{
    bool binaries = true, verbose = true, out = false,
         calc_e = true, Mcut=true, mcut=true;
    int nzones = 10, nstep = 1;
    bool auto_detect = false;
    real mmin = 0.1, m_min = 0.1;
    real mmax = 1000, m_max = 1000;
    real dt = 0;

    extern char *poptarg;
    int c;
    char* param_string = "N:n:MmeoT:v";	// Note: "v" removed because only the
					// "verbose" option currently works.

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {

	    case 'e': calc_e = !calc_e;
		      break;
	    case 'n': nzones = atoi(poptarg);
		      break;
	    case 'N': nstep = atoi(poptarg);
		      break;
	    case 'M': Mcut = !Mcut;
		      break;
	    case 'm': mcut = !mcut;
		      break;
	    case 'a': auto_detect = !auto_detect;
		      break;
	    case 'T': dt = atof(poptarg);
		      break;
	    case 'o': out = !out;
		      break;
	    case 'v': verbose = !verbose;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Loop over input until no more data remain.

    dyn *b;

    real to_volume = 4*PI/3.;
    real n_rho1, n_rhoi, n_rhoj;
    real r_tot, m_tot, v_rel;      
    int n, ncore;
    real mass_int;
    vector density_center;
    real rcore=0, vcore=0;
    real kT;
    real M_cutoff, m_cutoff;
    real m_total;

    //    real dmass  = log10(mmax-mmin)/(1.*nzones);
        real dmass  = (mmax-mmin)/(1.*nzones);

    real mmean, rmean, vdisp;
    //     real *core_radius = (real *) calloc(nzones+1, 8);
    //     int  *nstars = (int *) calloc(nzones+1, 4);
    //     real *m_mean = (real *) calloc(nzones+1, 8);
    //     real *v_disp = (real *) calloc(nzones+1, 8);
    //     real *r_mean = (real *) calloc(nzones+1, 8);
    
    real *core_radius = new real[nzones+1];
    real  *ncstars = new real[nzones+1];
    real  *nstars = new real[nzones+1];
    real *m_mean = new real[nzones+1];
    real *v_disp = new real[nzones+1];
    real *vc_disp = new real[nzones+1];
    real *r_mean = new real[nzones+1];

    real  *new_nstars = new real[nzones+1];
    real  *old_nstars = new real[nzones+1];
    real  *loss = new real[2*nzones+1];
    real  *gain = new real[2*nzones+1];
    real   new_mass;
    int    km;
	
    real geometric, grav_foc, rate, total_rate;
    
    while (b = get_dyn(cin)) {

      total_rate = rcore = vcore=0;
      
      for (int i=0; i<=nzones; i++) {
	core_radius[i] = m_mean[i] = v_disp[i] = vc_disp[i] = r_mean[i] = 0;
	ncstars[i] = nstars[i] = 0;
	old_nstars[i] = new_nstars[i] = 0;
      }
      for (int i=0; i<=2*nzones; i++)
	loss[i]=gain[i]=0;

      b->get_starbase()->print_stellar_evolution_scaling(cerr);
      cerr << "Time = " << b->get_system_time() << endl;
      
      // Count stars and find mass extremes.
      n=0;
      if (auto_detect) {
	m_min = VERY_LARGE_NUMBER;
	m_max = -m_min;
	for_all_daughters(dyn, b, bi) {
	  n++;
	  m_min = Starlab::min(m_min, bi->get_mass());
	  m_max = Starlab::max(m_max, bi->get_mass());
	}
      }
      else {
	for_all_daughters(dyn, b, bi) 
	  n++;
	m_min = b->get_starbase()->conv_m_star_to_dyn(m_min);
	m_max = b->get_starbase()->conv_m_star_to_dyn(m_max);
      }
      //      mass_int = (log10(m_max)-log10(m_min))/(1.0*nzones-1);
      mass_int = (m_max-m_min)/(1.0*nzones-1);

      //rcore = getrq(b->get_dyn_story(), "core_radius");
      //density_center = getrq(b->get_dyn_story(), "density_center_pos");
      
      real* local_density = new real[n+1];
      if(local_density==NULL) {
	cerr << "Not enough memory to allocate local_density[" << n << "].\n";
	exit(1);
      }
      for(int i=0; i<n; i++)
	local_density[i] = 0;

      print_core_parameters(b, density_center, rcore, ncore, vcore,
			    n, local_density, m_min, m_max, verbose);
      print_energies(b, kT, m_min, m_max, verbose);
      print_relaxation_time(b, m_min, m_max, verbose);
      print_lagrangian_radii(b, density_center, 0, m_min, m_max, 0);

      core_radius[0] = Rsun_pc * b->get_starbase()->conv_r_dyn_to_star(rcore);
      nstars[0] = n;
      ncstars[0] = ncore;
      vc_disp[0] = vcore;

      get_composition(b, m_min, m_max, n, mmean, rmean, vdisp);
      m_mean[0] = mmean;
      v_disp[0] = vdisp;
      r_mean[0] = rmean;
      
      //PRC(nstars[0]);PRC(m_min);PRC(m_max);PRL(mass_int);
      //PRC(rcore);PRC(m_mean[0]);PRC(r_mean[0]);PRC(v_disp[0]);

      for (int i=0; i<nzones; i++) {

	//	m_cutoff = pow(10., log10(m_min)+i*mass_int);
	//	M_cutoff = pow(10., log10(m_min)+(i+1)*mass_int);
	m_cutoff = (m_min)+i*mass_int;
	M_cutoff = (m_min)+(i+1)*mass_int;

	if (!Mcut) M_cutoff = m_max;
	if (!mcut) m_cutoff = m_min;

	for(int is=0; is<nstars[0]; is++)
	  local_density[is] = 0;

	mmean = vdisp = rmean = 0;
	get_composition(b, m_cutoff, M_cutoff, n, mmean, rmean, vdisp); 
	nstars[i+1] = n;
	//	m_mean[i+1] = mmean;
	m_mean[i+1] = M_cutoff;
	v_disp[i+1] = vdisp;
	r_mean[i+1] = rmean;
	
	if (n > 6 &&
	    !(n==nstars[i] && (!Mcut || !mcut))) {

	  rcore = 0;
	  print_core_parameters(b, density_center, rcore, ncore, vcore,
				n, local_density, m_cutoff, M_cutoff, verbose);
	  core_radius[i+1] = Rsun_pc * b->get_starbase()
	                                ->conv_r_dyn_to_star(rcore);
	  ncstars[i+1] = ncore;
	  vc_disp[i+1] = vcore;

	  print_energies(b, kT, m_cutoff, M_cutoff, verbose);
	  print_relaxation_time(b, m_cutoff, M_cutoff, verbose);
	  print_lagrangian_radii(b, density_center, 0, m_cutoff, M_cutoff, 0);

	  //PRC(n);PRC(core_radius[0]);PRL(core_radius[i+1]);

	  n_rho1 = 1./(to_volume * pow(core_radius[0], 3));
	  n_rhoi = ncstars[i+1]/(to_volume * pow(core_radius[i+1], 3));
	  r_tot  = r_mean[0] + r_mean[i+1];
	  m_tot  = m_mean[0] + m_mean[i+1];
	  v_rel  = sqrt(pow(vc_disp[0], 2) + pow(vc_disp[i+1], 2));

	  //PRC(vc_disp[0]);PRC(vc_disp[i+1]);PRL(v_rel);
	  
	  geometric = 6.31e-9 * v_rel * pow(r_tot, 2);
	  grav_foc = 3.61e-3 * m_tot * r_tot / v_rel;
	  rate = (n_rho1/1000.) * (n_rhoi/1000)
	                 * pow(core_radius[0], 3)
	                 * (grav_foc + geometric);

	  if (verbose) {
	    cerr << "     Stars: (N, r, m, v) = ( "
		 << nstars[i+1] <<", "
		 << r_mean[i+1] <<", "
		 << m_mean[i+1] <<", "
		 << v_disp[i+1] <<" ) [solar/km/s]\n";
	    cerr << "     rate = "
		 << rate << " [per Myr]\n"
	         << "          (geo, gf) = ( "
		 << 100 * geometric/(geometric+grav_foc)
		 << ", " << 100 * grav_foc/(geometric+grav_foc)
		 << " ) [\%]\n\n";
	  }
	  else {
	    PRI(4); PRC(m_mean[i+1]);PRC(v_disp[i+1]);PRL(r_mean[i+1]);
	    PRI(4);PRC(geometric); PRL(grav_foc);
	    PRI(4);PRL(rate);
	    cerr << endl;

	  }

	}
      }
      
      // Now print the encounter probability in a mass-mass plane.
      cerr << log10(b->get_starbase()
		     ->conv_m_dyn_to_star(m_min)) <<" "
	   << log10(b->get_starbase()
		     ->conv_m_dyn_to_star(m_max)) << " "
	   << nzones << endl;
      cerr << log10(b->get_starbase()
		     ->conv_m_dyn_to_star(m_min)) <<" "
	   << log10(b->get_starbase()
		     ->conv_m_dyn_to_star(m_max)) << " "
	   << nzones << endl;
      
      for (int im=1; im<nzones+1; im++) {
	for (int jm=im; jm<nzones+1; jm++) {

	  //PRC(im);PRC(ncstars[im]);PRC(m_mean[im]);PRL(core_radius[im]);
	  //PRC(jm);PRC(ncstars[jm]);PRC(m_mean[jm]);PRL(core_radius[jm]);

	  rate = 0;
	  
	  if (core_radius[im]>0) {
	    if(core_radius[jm]>0) 
	      
	    rate = encounter_rate(im,  ncstars[im], m_mean[im], r_mean[im], 
			   vc_disp[im], core_radius[im], 
			   jm, ncstars[jm], m_mean[jm], r_mean[jm],
			   vc_disp[jm], core_radius[jm]);
	    else 
	      
	    rate = encounter_rate(im, ncstars[im], m_mean[im], r_mean[im], 
			   vc_disp[im], core_radius[im], 
			   jm, nstars[jm]*core_radius[0], 
			   m_mean[jm], r_mean[jm],
			   vc_disp[0], core_radius[0]);
	  }
	  else if(core_radius[jm]>0) 
	      
	    rate = encounter_rate(im, nstars[im]*core_radius[0], 
			   m_mean[im], r_mean[im], 
			   vc_disp[0], core_radius[0], 
			   jm, ncstars[jm], m_mean[jm], r_mean[jm],
			   vc_disp[jm], core_radius[jm]);	  

	  total_rate += rate;
	  }

	}
    }
    //rmtree(b);
      cerr << 0 << " " << 0 << " "
	   << 0 << " " << 0 << " "
	   << 0 << " " << 0 << " "
	   << 0 << endl;
      cerr << "     Total encounter rate = "
	   << total_rate << " [per Myr]" << endl;
      


      //Fill in the Gaps.
      real mc_tot = 0;
      real new_mc_tot = 0;
      for (km=1; km<nzones+1; km++) {
	//PRC(km);PRC(ncstars[km]);PRC(m_mean[km]);PRL(core_radius[km]);
	
	if (core_radius[km]<=0) {
	    core_radius[km] = core_radius[0];
	    new_nstars[km] = pow(1.0*nstars[km], 1./3.);
	    //m_mean[km] = pow(10, (km-1)*dmass) -1 + mmin;
	    r_mean[km] = pow(m_mean[km], 0.7);
	    vc_disp[km] = vc_disp[0];
	}
	else 
	  new_nstars[km] = ncstars[km];
	
	m_mean[km] = km*dmass;
	nstars[km] -= new_nstars[km];
	mc_tot += new_nstars[km] * m_mean[km];

	//PRC(km);PRC(new_nstars[km]);PRC(m_mean[km]);PRL(core_radius[km]);
      }
      for (km=1; km<nzones+1; km++) {
	old_nstars[km] = new_nstars[km];
      }

      real nstar_loss = 0;
      for(int ii=1; ii<nstep; ii++) {
      for (int im=1; im<nzones+1; im++) {
	for (int jm=im; jm<nzones+1; jm++) {

	  rate = encounter_rate(im,  new_nstars[im], m_mean[im], r_mean[im], 
				vc_disp[im], core_radius[im], 
				jm, new_nstars[jm], m_mean[jm], r_mean[jm],
				vc_disp[jm], core_radius[jm]);

      
	  //PRC(im);PRC(jm);PRC(m_mean[im]);PRC(m_mean[jm]);
	  //PRC(new_nstars[im]);PRL(rate);
	    new_mass = m_mean[im] + m_mean[jm];
	    //	    km = int((log10(new_mass)+log10(mmin))/dmass);
	    //      km = int((log10(new_mass) - log10(mmin))/dmass);
	    //km = int(1 + log(new_mass - mmin + 1)/dmass);
	    //km = 1 + int(new_mass - mmin/dmass);
	    km = Starlab::min(2*nzones, im + jm);
	    //m_mean[km] = pow(10, (km-1)*dmass) -1 + mmin;
	    
	    real ncoll = Starlab::max(0., Starlab::min(rate * dt,
				 Starlab::min(new_nstars[im]-loss[im],
				     new_nstars[jm]-loss[im])));

	    //PRC(new_mass);PRC(km);PRC(rate);PRL(ncoll);
	    loss[im] += ncoll;
	    loss[jm] += ncoll;
	    gain[km] += ncoll;
	}
      }

      real total_loss = 0;
      real total_gain = 0;
      for (int im=1; im<nzones+1; im++) {
	total_loss += loss[im] * m_mean[im];
	total_gain += gain[im] * m_mean[im];
	new_nstars[im] = new_nstars[im] - loss[im] + gain[im];
	nstar_loss = nstar_loss + gain[im] - loss[im];
	
      }
      for (int im=0; im<=2*nzones; im++)
	loss[im]=gain[im]=0;

      cerr << "Total loss = " << total_loss
	   << " total gain = " << total_gain << endl;
      }

      cerr << "Old and new mass function: "<<endl;
      for (km=1; km<nzones+1; km++) {
	new_mc_tot += new_nstars[km] * m_mean[km];
	if(old_nstars[km]<=1.e-7)
	  old_nstars[km]=0;
	if(new_nstars[km]<=1.e-7)
	  new_nstars[km]=0;
	if(old_nstars[km]>0 || new_nstars[km]>0)
	  cerr << km <<" "<< m_mean[km] <<" "
	       << old_nstars[km] << " " << new_nstars[km] << endl;
      }

      cerr << " N star lost = " << nstar_loss << endl;
      cerr << "Total mass goes from: "<< mc_tot << " to " << new_mc_tot <<endl;
      //	  new_nstars[km] = nstars[km] + rate * dt;

    // Now evolvefor one timestep and compute the
    // new mass spectrum

}

#endif
