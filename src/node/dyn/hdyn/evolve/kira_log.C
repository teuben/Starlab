
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

// Functions associated with log output, snapshots, etc.
//
// Externally visible functions:
//
//	void print_dstar_params
//	bool print_dstar_stats
//	void get_energies_with_tidal
//	void refine_cluster_mass
//	void print_statistics
//	void update_cpu_counters
//	void log_output
//	void snap_output

#include "hdyn.h"
#include "star/dstar_to_kira.h"

typedef  struct {
    real  dt;
    real  dt_a;
    hdyn*  b;
} dt_pair, *dt_pair_ptr;

local int compare_dt(const void * pi, const void * pj)	  // increasing dt
{
    if (((dt_pair_ptr) pi)->dt < ((dt_pair_ptr) pj)->dt)
        return -1;
    else if (((dt_pair_ptr)pi)->dt > ((dt_pair_ptr)pj)->dt)
        return +1;
    else
        return 0;
}

typedef  struct {
    real  count;
    hdyn*  b;
} count_pair, *count_pair_ptr;

local int compare_steps(const void * pi, const void * pj) // decreasing count
{
    if (((count_pair_ptr) pi)->count > ((count_pair_ptr) pj)->count)
        return -1;
    else if (((count_pair_ptr)pi)->count < ((count_pair_ptr)pj)->count)
        return +1;
    else
        return 0;
}

int get_effective_block(real dt)
{
    // For general timestep dt, find the effective block it is in, by
    // determining the smallest k for which the step can be written as
    //
    //		dt = n * 2^{-k}
    //
    // for integral n.

    int k = 0;
    real dtb = 1;

    while (fmod(dt, dtb) != 0 && k < 50) {
	k++;
	dtb /= 2;
    }

    return k;
}

#define NTOP 10
#define NBOT 10

local void print_timestep_stats(hdyn* b)
{
    // Time step statistics.

    real dt_mean = 0, dt_min = VERY_LARGE_NUMBER, dt_max = -dt_min;

    real n_dt = 0;		// actually these are integer counters
    real n_dt_p = 0;
    real n_dt_u = 0;
    real total_steps = 0;
    real total_steps_a = 0;
    real total_steps_p = 0;

    // Timestep histograms:

    int count[30],		// all nodes, perturbed steps
        count_a[30],		// all nodes, actual steps
        count_p[30],		// all perturbed nodes
        count_u[30],		// all unperturbed nodes, unperturbed steps
        count_ue[30],		// all unperturbed nodes, effective blocks
        count_next[30];		// next_time distribution

    int i;
    for (i = 0; i < 30; i++)
	count[i] = count_a[i]
	    	 = count_p[i]
		 = count_u[i]
		 = count_ue[i]
		 = count_next[i]
		 = 0;

    // Count the total number of particles/nodes.

    int n = 0;
    for_all_nodes(hdyn, b, bb)  n++;
    n--;

    dt_pair_ptr list_dt  = new dt_pair[n];
    dt_pair_ptr list_dtb = new dt_pair[n];

    // Total timestep counters:

    count_pair_ptr steps_inst  = new count_pair[n];
    count_pair_ptr steps_delta = new count_pair[n];
    count_pair_ptr f_dir_delta = new count_pair[n];
    count_pair_ptr f_ind_delta = new count_pair[n];

    n = 0;
    int np = 0, nd = 0, ni = 0;

    real delta_steps_tot = 0;
    real delta_f_dir_tot = 0;
    real delta_f_ind_tot = 0;

    bool print_dtb = false;

    for_all_nodes(hdyn, b, bb)
	if (bb != b) {

	    real dt = bb->get_timestep();

	    // This dt is the timestep of a top-level or perturbed binary
	    // node, and the perturbed timestep of an unperturbed binary.

	    // Actual timestep:

	    real dt_a = dt;
	    if (bb->get_kepler())
		dt_a = bb->get_unperturbed_timestep();

	    // Don't include younger binary sisters in *any* statistic.

	    if (bb->is_top_level_node() || bb->get_elder_sister() == NULL) {

		// Basic stats (all steps):

		dt_mean += dt_a;
		dt_min = min(dt_min, dt_a);
		dt_max = max(dt_max, dt_a);

		n_dt += 1;
		real steps = rint(1/dt);
		real steps_a = rint(1/dt_a);
		total_steps += steps;
		total_steps_a += steps_a;

		// Update perturbed_steps counters.

		if (!bb->get_kepler()) {

		    list_dt[np].dt = dt;
		    list_dt[np].dt_a = dt;
		    list_dt[np].b = bb;
		    steps_inst[np].count = steps;
		    steps_inst[np].b = bb;
		    np++;
		}

		// 1. All perturbed steps, including unperturbed
		//    binaries: count.

		i = 0;
		real dtbin = 1;
		while (dt < dtbin) {
		    i++;
		    dtbin *= 0.5;
		}
		if (i > 29) i = 29;
		count[i]++;

		if (!bb->get_kepler()) {

		    // 2. Perturbed nodes only: count_p.

		    total_steps_p += steps;
		    n_dt_p += 1;
		    count_p[i]++;

		}

		// 3. All actual steps, including unperturbed
		//    binaries: count_a.

		i = 0;
		dtbin = 1;
		while (dt_a < dtbin) {
		    i++;
		    dtbin *= 0.5;
		}
		if (i > 29) i = 29;
		count_a[i]++;

		if (bb->get_kepler()) {

		    // 4a. Unperturbed nodes only, unperturbed steps: count_u.

		    n_dt_u += 1;
		    count_u[i]++;

		    // 4b. Unperturbed nodes only, effective blocks: count_ue.

		    int kb = get_effective_block(dt_a);
		    if (kb > 29) kb = 29;
		    count_ue[kb]++;
		}

		// 5. Next_time, all nodes, actual steps: count_next.

		real dt_next = bb->get_next_time() - bb->get_system_time();
		int kb = get_effective_block(dt_next);
		real dtb = pow(2.0, -kb);

		list_dtb[n].dt = -kb;
		list_dtb[n].dt_a = dt_a;
		list_dtb[n].b = bb;

		if (bb->get_kepler() || dt_a != dtb) {

		    if (!print_dtb) {
			cerr << endl << "  Actual and block timesteps:"
			     << endl;
			print_dtb = true;
		    }

		    if (bb->get_kepler())
			cerr << "    unperturbed ";
		    else
			cerr << "    ***mismatch ";

		    cerr << bb->format_label() << ":  ";
		    PRC(kb), PRC(dt_a); PRL(dtb);
		}

		if (kb > 29) kb = 29;
		count_next[kb]++;		// next_time

		// Total timestep/force counters:

		// Instantaneous:

		steps_delta[n].count = bb->get_steps();
		if (find_qmatch(bb->get_dyn_story(), "prev_steps"))
		    steps_delta[n].count -= getrq(bb->get_dyn_story(),
						  "prev_steps");
		putrq(bb->get_dyn_story(), "prev_steps", bb->get_steps());
		steps_delta[n].b = bb;
		delta_steps_tot += steps_delta[n].count;

		n++;

		// Direct (top-level/GRAPE):

		f_dir_delta[nd].count = bb->get_direct_force();
		if (find_qmatch(bb->get_dyn_story(), "prev_f_dir"))
		    f_dir_delta[nd].count -= getrq(bb->get_dyn_story(),
						   "prev_f_dir");
		putrq(bb->get_dyn_story(), "prev_f_dir",
		      bb->get_direct_force());
		f_dir_delta[nd].b = bb;
		delta_f_dir_tot += f_dir_delta[nd].count;
		nd++;

		// Indirect (low-level/front-end):

		f_ind_delta[ni].count = bb->get_indirect_force();
		if (find_qmatch(bb->get_dyn_story(), "prev_f_ind"))
		    f_ind_delta[ni].count -= getrq(bb->get_dyn_story(),
						   "prev_f_ind");
		putrq(bb->get_dyn_story(), "prev_f_ind",
		      bb->get_indirect_force());
		f_ind_delta[ni].b = bb;
		delta_f_ind_tot += f_ind_delta[ni].count;
		ni++;
	    }
	}

    //
    //----------------------------------------------------------------------
    //
    // Print out the results:

    // First print the timestep histograms...

    cerr << endl << "  Timesteps (younger binary components excluded,"
	 << endl << "             mean = " << dt_mean/max((int)rint(n_dt), 1)
	 << "  min = " << dt_min << "  max = " << dt_max << "):"
	 << endl;

    if (n_dt_u > 0)
	cerr << endl << "  all nodes (perturbed steps, bin #1 = 1, dt /= 2):"
	     << endl << "        ";
    else
	cerr << endl << "  all nodes (actual steps, bin #1 = 1, dt /= 2):"
	     << endl << "        ";

    for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count[i]);
    cerr << endl << "        ";
    for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count[i]);
    cerr << endl << "        ";
    for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count[i]);

    cerr << endl << "  total = " << n_dt
	 << ",  weighted timestep sum = "
	 << total_steps/n_dt << endl;

    if (n_dt_u > 0) {

	cerr << endl << "  all nodes (actual steps):"
	     << endl << "        ";
	for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count_a[i]);
	cerr << endl << "        ";
	for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count_a[i]);
	cerr << endl << "        ";
	for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count_a[i]);

	cerr << endl << "  total = " << n_dt
	     << ",  weighted timestep sum = "
	     << total_steps_a/n_dt << endl;

	if (n_dt_p > 0) {
	    cerr << endl << "  excluding unperturbed systems:"
		 << endl << "        ";
	    for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count_p[i]);
	    cerr << endl << "        ";
	    for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count_p[i]);
	    cerr << endl << "        ";
	    for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count_p[i]);

	    cerr << endl << "  total = " << n_dt_p
		 << ",  weighted timestep sum = "
		 << total_steps_p/n_dt_p
		 << endl;
	}

	cerr << endl << "  unperturbed systems only:"
	     << endl << "        ";
	for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count_u[i]);
	cerr << endl << "        ";
	for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count_u[i]);
	cerr << endl << "        ";
	for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count_u[i]);

	cerr << endl << "  total = " << n_dt_u << endl;

	cerr << endl << "  unperturbed systems only (effective blocks):"
	     << endl << "        ";
	for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count_ue[i]);
	cerr << endl << "        ";
	for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count_ue[i]);
	cerr << endl << "        ";
	for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count_ue[i]);
	cerr << endl;
    }

    cerr << endl
	 << "  next_step distribution (all nodes):"
	 << endl << "        ";
    for (i =  0; i < 10; i++) fprintf(stderr, " %5d", count_next[i]);
    cerr << endl << "        ";
    for (i = 10; i < 20; i++) fprintf(stderr, " %5d", count_next[i]);
    cerr << endl << "        ";
    for (i = 20; i < 30; i++) fprintf(stderr, " %5d", count_next[i]);
    cerr << endl;    

    // ...then the list of shortest steps and force counters.

    if (n > NBOT+NTOP) {

	// Print out NBOT shortest perturbed time steps.

	cerr << endl << "  " << NBOT << " shortest perturbed time steps:"
	     << endl << endl;

	qsort((void *)list_dt, (size_t)np, sizeof(dt_pair), compare_dt);

	for (i = 0; i < NBOT; i++)
	    fprintf(stderr, "    %3d.     %-15s  %.4e\n",
		    i+1,
		    list_dt[i].b->format_label(),
		    list_dt[i].dt);

	// Print out NBOT shortest effective block steps.

	cerr << endl << "  " << NBOT << " shortest effective block steps:"
	     << endl << endl;

	qsort((void *)list_dtb, (size_t)n, sizeof(dt_pair), compare_dt);

	for (i = 0; i < NBOT; i++)
	    fprintf(stderr, "    %3d.     %-15s  %3d   %.4e   %.4e\n",
		    i+1,
		    list_dtb[i].b->format_label(),
		    -(int)list_dtb[i].dt,	   // "dt" here is really -kb
		    pow(2.0, list_dtb[i].dt),
		    list_dtb[i].dt_a);

	// Print out top NTOP particles in (1) instantaneous cost (steps),
	//				   (2) integrated cost since the
	//				       last output.
	//				   (3) direct ("GRAPE") forces
	//				       since the last output
	//				   (4) indirect ("front-end") forces
	//				       since the last output.

	// Sort the arrays in order of decreasing count:

	qsort((void *)steps_inst,  (size_t)np, sizeof(count_pair),
	      compare_steps);
	qsort((void *)steps_delta, (size_t)n,  sizeof(count_pair),
	      compare_steps);
	qsort((void *)f_dir_delta, (size_t)nd, sizeof(count_pair),
	      compare_steps);
	qsort((void *)f_ind_delta, (size_t)ni, sizeof(count_pair),
	      compare_steps);

	cerr << endl << "  Top " << NTOP
	     << " nodes by instantaneous (excl. unperturbed)"
	     << " and delta (all) steps:"
	     << endl << endl;

	fprintf(stderr, "             %-26s      %-30s\n",
		"--- instantaneous ---",
		"------- delta -------");
	fprintf(stderr, "            %-26s      %-30s\n\n",
		"(normalized to dt = 1) ",
		"(since last log output)");

	char tmp[1024];

	for (i = 0; i < NTOP; i++) {

	    if (i < np)
		sprintf(tmp, "%3d.     %-10s %10.0f",
			i+1,
			steps_inst[i].b->format_label(),
			steps_inst[i].count);
	    else
		strcpy(tmp, " ");
	    fprintf(stderr, "    %-35s", tmp);

	    if (i < n)
		sprintf(tmp, "%-10s %10.0f",
			steps_delta[i].b->format_label(),
			steps_delta[i].count);
	    else
		strcpy(tmp, " ");
	    fprintf(stderr, "      %-30s\n", tmp);
	}
	cerr << endl << "  total delta_steps = " << delta_steps_tot << endl;

	cerr << endl << "  Top " << NTOP
	     << " nodes by delta(force calculations):"
	     << endl << endl;

	fprintf(stderr, "             %-26s      %-30s\n",
		"------- direct ------",
		"------ indirect -----");

	for (i = 0; i < NTOP; i++) {

	    if (i < nd)
		sprintf(tmp, "%3d.     %-10s %10.0f",
			i+1,
			f_dir_delta[i].b->format_label(),
			f_dir_delta[i].count);
	    else
		strcpy(tmp, " ");
	    fprintf(stderr, "    %-35s", tmp);

	    if (i < ni && f_ind_delta[i].count > 0)
		sprintf(tmp, "%-10s %10.0f",
			f_ind_delta[i].b->format_label(),
			f_ind_delta[i].count);
	    else
		strcpy(tmp, " ");
	    fprintf(stderr, "      %-30s\n", tmp);
	}
	cerr << endl << "  total delta_f_dir = " << delta_f_dir_tot
	     << "  delta_f_ind = " << delta_f_dir_tot
	     << endl;
    }

    delete [] list_dt;
    delete [] steps_inst;
    delete [] steps_delta;
    delete [] f_dir_delta;
    delete [] f_ind_delta;
}

void print_dstar_params(dyn* b)			// b is binary center of mass
{
    if (((hdyn*)b)->get_use_dstar())
	print_binary_dstars(b);
}

bool print_dstar_stats(dyn* b, bool mass_function,    // b is binary 
		       vector center, bool verbose)   // center of mass
{
    if (((hdyn*)b)->get_use_dstar()) {
	dstar_stats(b, mass_function, center, verbose);

	return true;
    }

    return false;
}

void get_energies_with_tidal(dyn* b, real eps2, 
			     real &potential, real &kinetic, real &total,
			     bool cm)
{
    // Coerce calculate_energies_with_tidal into the same calling
    // sequence as dyn::calculate_energies, for use by sys_stats...
    // Discard eps2 (--> 0).

    calculate_energies_with_tidal((hdyn*)b, kinetic, potential, total, cm);
}

#define G_TOL 1.e-6

local real solve3(real a, real b, real c)
{
    // Iteratively solve
    //
    //		a z^3 + b z + c  =  0
    //
    // where we know a > 0, b > 0, c < 0.

    // Proceed by Newton-Raphson iteration.  Since a > 0, b > 0, so
    // there are are no turning points to complicate the logic.

    real z = 2*max(pow(abs(c/a), 1.0/3), sqrt(abs(b/a)));
    real g = c + z * (b + a*z*z);
    real gpr = b + 3*a*z*z;		// always > 0

    while (abs(g) > G_TOL) {
	z = z - g/gpr;
	g = c + z * (b + a*z*z);
	gpr = b + 3*a*z*z;
    }

    return z;
}

#define M_TOL 1.e-4
#define M_ITER_MAX 20

void refine_cluster_mass(hdyn *b)
{
    if (b->get_tidal_field() == 0) return;

    // Self-consistently determine the total mass within the outermost
    // closed zero-velocity surface.  Use a point-mass approximation for
    // the cluster potential and iterate until the actual mass within
    // the surface agrees with the mass used to generate the surface.
    // The total potential is
    //
    //		phi  =  -GM/r + (alpha1 x^2 + alpha3 z^2) / 2
    //
    // where we measure everything relative to lagr_pos (which is
    // the density center, if known, and the modified center of mass
    // otherwise).  The Jacobi radius for this field is
    //
    //		r_J  =  (-GM/alpha1)^{1/3}
    //
    // and the potential of the last closed surface is
    //
    //		phi_J  =  1.5 alpha1 r_J^2.

    vector center(0);
    if (find_qmatch(b->get_dyn_story(), "lagr_pos"))
	center = getvq(b->get_dyn_story(), "lagr_pos");

    real M_inside = total_mass(b), M = -1;
    real r_J, r_x2, r_y2, r_z2, r_max_inside;
    int  N_inside, iter = 0;

    real M_J, M_x, M_y, M_z;
    int  N_J, N_x, N_y, N_z;

    cerr << endl << "  refine_cluster_mass: getting M by iteration"
	 << "; total system mass = " << M_inside
	 << endl;

    while (iter++ < M_ITER_MAX
	   && M_inside > 0
	   && abs(M_inside/M - 1) > M_TOL) {

	M = M_inside;
	M_inside = 0;
	M_J = M_x = M_y = M_z = 0;
	N_J = N_x = N_y = N_z = 0;

	r_J = pow(-M/b->get_alpha1(), 0.3333333);

	real phi_J = 1.5 * b->get_alpha1() * r_J * r_J;
	real r_max = 0;

	r_x2 = square(r_J);		// zero-velocity surface crosses x-axis
	r_y2 = square(-M/phi_J);	// zero-velocity surface crosses y-axis

	// zero-velocity surface crosses z-axis where
	//
	//	-GM/z + 0.5 alpha3 z^2  =  phi_J
	// so
	//	0.5 alpha3 z^3 -phi_J z -GM  =  0

	r_z2 = square(solve3(0.5*b->get_alpha3(), -phi_J, -M));

	N_inside = 0;
	r_max_inside = 0;

	for_all_daughters(hdyn, b, bb) {

	    vector dx = bb->get_pos() - center;

	    real x = dx[0];
	    real y = dx[1];
	    real z = dx[2];
	    real r = abs(dx);

	    if (r < r_J) {

		N_J++;
		M_J += bb->get_mass();

		if (r == 0 
		    || -M/r + 0.5 * (b->get_alpha1()*x*x + b->get_alpha3()*z*z)
			    < phi_J) {
		    N_inside++;
		    M_inside += bb->get_mass();
		    r_max_inside = max(r, r_max_inside);
		}

	    }
	    r_max = max(r, r_max);

	    // Count projected masses and numbers.

	    real xy = x*x + y*y;
	    real xz = x*x + z*z;
	    real yz = y*y + z*z;

	    if (xy < r_x2) {		// z projection
		M_x += bb->get_mass();
		N_x++;
	    }
	    if (xz < r_x2) {		// y projection
		M_x += bb->get_mass();
		N_x++;
	    }

	    if (yz < r_y2) {		// x projection
		M_y += bb->get_mass();
		N_y++;
	    }
	    if (xy < r_y2) {		// z projection
		M_y += bb->get_mass();
		N_y++;
	    }

	    if (yz < r_z2) {		// x projection
		M_z += bb->get_mass();
		N_z++;
	    }
	    if (xz < r_z2) {		// y projection
		M_z += bb->get_mass();
		N_z++;
	    }
	}

	M_x /= 2;
	N_x /= 2;
	M_y /= 2;
	N_y /= 2;
	M_z /= 2;
	N_z /= 2;

#if 0
	fprintf(stderr, "    %2d  ", iter);
	PRC(M); PRC(r_J); PRL(r_max);
	PRI(8); PRC(N_inside); PRC(M_inside); PRL(r_max_inside);
#endif

    }

    if (iter >= M_ITER_MAX)
	cerr << "    (too many iterations)" << endl;

    cerr << "  within last zero-velocity surface:"
	 << "  M = " << M_inside << "  N = " << N_inside
	 << endl
	 << "                                    "
	 << "  r_J = " << r_J << "  r_max = " << r_max_inside
	 << endl
	 << "  within r_J:                       "
	 << "  M = " << M_J << "  N = " << N_J
	 << endl
	 << "  within r_J (projected):           "
	 << "  M = " << M_x << "  N = " << N_x
	 << endl
	 << "  within r_y (projected):           "
	 << "  M = " << M_y << "  N = " << N_y << "  r_y = " << sqrt(r_y2)
	 << endl
	 << "  within r_z (projected):           "
	 << "  M = " << M_z << "  N = " << N_z << "  r_z = " << sqrt(r_z2)
	 << endl;
}

void print_statistics(hdyn* b,
		      int long_binary_output)		// default = 2
{
    // General system statistics (omit time output
    // and O(N^2) operations):

    sys_stats(b,
	      0.5,			// energy cutoff
	      true,			// verbose output
	      (long_binary_output > 0),	// print binaries
	      (long_binary_output > 1),	// full binary_output
	      2,			// which lagr
	      false,			// don't print time
	      false,			// don't compute energy
	      false,			// don't allow n^2 ops
	      get_energies_with_tidal,	// energy calculation function
	      print_dstar_params,	// allow access to dstar_to_kira
	      print_dstar_stats);	// from dyn sys_stats function

    // Statistics on time steps and bottlenecks:

    print_timestep_stats(b);
    print_sort_counters();
}

local int n_bound(hdyn* b, vector& cod_vel)
{
    int nb = 0;
    for_all_daughters(hdyn, b, bb)
	if (0.5*square(bb->get_vel()-cod_vel) + bb->get_pot() < 0)
	    nb += bb->n_leaves();
    return nb;
}

local real m_bound(hdyn* b, vector& cod_vel)
{
    real mb = 0;
    for_all_daughters(hdyn, b, bb)
	if (0.5*square(bb->get_vel()-cod_vel) + bb->get_pot() < 0)
	    mb += bb->get_mass();
    return mb;
}

local void get_n_and_m_bound(hdyn* b, vector& cod_vel, int& nb, real& mb)
{
    nb = 0;
    mb = 0;
    for_all_daughters(hdyn, b, bb)
	if (0.5*square(bb->get_vel()-cod_vel) + bb->get_pot() < 0) {
	    nb += bb->n_leaves();
	    mb += bb->get_mass();
    }
}



void update_cpu_counters(hdyn * b)
{
    real last_cpu = b->get_kira_counters()->cpu_time;
    b->get_kira_counters()->cpu_time = cpu_time();
    b->get_kira_counters()->total_cpu_time += (cpu_time() - last_cpu);
    write_counters_to_log(b);
}

void log_output(hdyn * b, real count, real steps,
		kira_counters* kc_prev,
		int long_binary_output)	// default = 2
{
    // Standard log output (to stderr, note).

    vector cod_pos, cod_vel;
    compute_com(b, cod_pos, cod_vel);

//#if 0

    // Suppressed for now while we are debugging other parts of
    // the GRAPE-6 interface (Steve, 7/00):

    compute_densities(b, cod_pos, cod_vel);	// does nothing unless
						// GRAPE is present
//#endif

    int n_bound;
    real m_bound;
    get_n_and_m_bound(b, cod_vel, n_bound, m_bound);
    
    cerr << endl;
    for (int k = 0 ; k < 40; k++) cerr << '-';
    cerr << endl << endl
	 << "Time = " << b->get_system_time()
	 << "  N = "    << b->n_leaves() << " (" << n_bound << " bound)"
	 << "  mass = " << total_mass(b) << " (" << m_bound << " bound)"
	 << endl
	 << "          block steps = " << count
	 << "  total steps = " << steps
	 << "  steps per block = " << steps/max(1.0, count)
	 << endl;

    if (b->get_use_sstar())
	print_sstar_time_scales(b);

    int p = cerr.precision(INT_PRECISION);
    print_recalculated_energies(b, true, true);
    cerr.precision(p);

    // Note: on return from print_recalculated_energies, hdyn::pot is
    // properly set, and includes the tidal field, if any.  The earlier
    // "bound" quantities may be suspect...

    if (steps > 0) {
	update_cpu_counters(b);
	print_counters(b->get_kira_counters(), kc_prev);
    }

    print_statistics(b, long_binary_output);

    refine_cluster_mass(b);		// won't be called from sys_stats;
					// make the output appear in the
					// same place as in the standalone
					// version

    cerr << flush;
}

void snap_output(hdyn * b, real steps, int& snaps,
		 bool reg_snap, bool last_snap,
		 char * snap_save_file,
		 xreal t, xreal ttmp, real t_end,
		 real& t_snap, real dt_snap, bool verbose)
{
    // Save a snapshot.

    update_cpu_counters(b);
    b->set_time(t);

    if (last_snap) {

	// Write (overwrite) snap_save_file.

	ofstream file(snap_save_file);

	if (file) {
	    put_node(file, *b, b->get_kira_options()->print_xreal);
	    file.close();
	    cerr << "Snapshot saved in file "
		 << snap_save_file
		 << " (save CPU = "
		 << cpu_time() - b->get_kira_counters()->cpu_time
		 << ")\n";
	} else {
	    if (snap_save_file == NULL)
		cerr << "No save file specified.\n";
	    else
		cerr << "Error opening save file "
		     << snap_save_file << "\n";
	}
    }

    if (reg_snap) {

	if (verbose) cerr << endl;
	    cerr << "Snapshot output #" << ++snaps
		 << " at steps = " << steps
		 << ", time = " << b->get_system_time() << endl;
	if (verbose) cerr << endl;

	put_node(cout, *b, b->get_kira_options()->print_xreal);
	cout << flush;

	// if (ttmp > t_end) exit(0);	// do this in kira, not here...
    }
}
