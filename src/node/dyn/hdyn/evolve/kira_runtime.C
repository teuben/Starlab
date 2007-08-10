
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// Implement run-time corrections/adjustments to kira.
//
// Externally visible functions:
//
//	void clean_up_files
//	bool check_file
//	bool check_runtime

#include "hdyn.h"
#include "kira_defaults.h"
#include <unistd.h>

local void clean_up_file(char* name)
{
    char tmp[80];

    // NOTE: this procedure may fail if the system is running
    // close to its memory limit.  Files can be opened and read,
    // and appropriate action taken, but the system call to delete
    // may not work.  This may represent a serious problem for
    // options that produce a lot of output. (Steve, 11/99)
    //
    // Fix is to use the "unlink" system command. (Steve, 6/00)

    // sprintf(tmp, "rm -f %s", name);
    // system(tmp);

    cerr << endl;
    if (unlink(name) == 0)
	cerr << "***** deleted file " << name << endl;
    else
	cerr << "***** error deleting file " << name << endl;

    ifstream file(name);
    if (file) {
	cerr << ">>>>> clean_up_file: warning: file " << name
	     << " NOT deleted" << endl;
	file.close();
    }
}

#define N_FILES 14

void clean_up_files()
{
    // List files to be explicitly cleaned up:

    char * files[N_FILES] = {"ALLOW_CHECK",
			     "ABORT",
			     "DIAG",
			     "DUMP",
			     "ENERGY",
			     "KILL_MULTIPLES",
			     "LOG",
			     "MULTIPLES",
			     "OPTIONS",
			     "PARAMS",
			     "PP",
			     "QT",
			     "SCHED",
			     "STOP"};

    // Clean up...

    for (int i = 0; i < N_FILES; i++) {
	ifstream file(files[i]);
	if (file) {
	    file.close();
	    clean_up_file(files[i]);
	}
    }
}

bool check_file(char* name,
		bool del)	// default = true
{
    ifstream file(name);
    if (file) {
	file.close();
	if (del) clean_up_file(name);
	return true;
    } else
	return false;
}

local void print_multiples(hdyn * b)
{
    cerr << endl << "***** In print_multiples:\n";
    int n_mult = 0;

    for_all_daughters(hdyn, b, bb)
	if (!bb->is_leaf())
	    if (bb->n_leaves() > 2) {
		n_mult++;
		cerr << "\nFound multiple system " << bb->format_label()
		     << endl;
		real eb = print_structure_recursive(bb);
		cerr << "Energy = " << eb << endl;
	    }

    if (n_mult == 0) cerr << "    (no multiples found)\n";

    cerr << endl;
}

local void force_step(hdyn * b, real t)
{
    real dt = b->get_timestep();
    int i = 0;
    while (b->get_next_time() > t && i++ < 10)
	b->set_timestep((dt /= 2));
    b->set_timestep((dt *= 2));
}

local void force_unperturbed_step(hdyn * b)
{
    kepler * kep = b->get_kepler();
    if (kep) {
	real dt_unp = b->get_unperturbed_timestep();
	while (b->get_next_time() > b->get_parent()->get_next_time())
	    b->set_unperturbed_timestep((dt_unp -= kep->get_period()));
	b->set_unperturbed_timestep((dt_unp += kep->get_period()));
    } else {
	cerr << "warning: " << b->format_label()
	     << " has no kepler pointer\n";
	force_step(b, b->get_parent()->get_next_time());
    }
}

local void adjust_vel(hdyn * b)
{
    real v = abs(b->get_vel());
    real r = abs(b->get_pos());

    b->set_vel(-v * b->get_pos() / r);	// Preserve v

    hdyn * s =  b->get_binary_sister();
    real mass_fac = -b->get_mass()/s->get_mass();

    s->set_vel(b->get_vel() * mass_fac);

    // Recompute acc and jerk (acc unchanged)...

    vec a_2b, j_2b;
    real p_2b;
    real d_min_sister = VERY_LARGE_NUMBER;
    hdyn * p_nn_sister = NULL;
    real eps2 = 0.0;
    b->calculate_partial_acc_and_jerk(b->get_parent(), b->get_parent(), b,
				      a_2b, j_2b, p_2b,
				      d_min_sister, p_nn_sister,
				      false, NULL, b);

//  b->set_acc_and_jerk_and_pot(mass_fac*a_2b, mass_fac*j_2b, -mass_fac*p_2b);

}

local bool kill_multiple(hdyn * b)
{
    // For now, adjust eccentricity of top-level node b only...

    hdyn * od = b->get_oldest_daughter();
    kepler * kep = od->get_kepler();

    if (kep) {

	if (kep->get_eccentricity() < 1) {

	    // Unperturbed system.  Reduce timesteps to advance
	    // system in scheduling, then modify kepler structure.

	    // Top-level CM step...

	    force_step(b, b->get_system_time());

	    // Top-level binary step...

	    force_unperturbed_step(od);

	    // Lower levels, if any...

	    for_all_nodes(hdyn, b, bb) {
		if (bb->get_parent() != b && bb->get_elder_sister() == NULL) {

		    if (bb->get_kepler() == NULL) {

			cerr << "warning: " << bb->format_label()
			     << " has no kepler pointer\n";

			force_step(bb, bb->get_parent()->get_next_time());

		    } else

			force_unperturbed_step(bb);
		}
	    }

	    // Modify top-level orbit....

	    kep->set_eccentricity(0.99999);

	    // Probably need to do more (reinitialize) here...	<<<<<<

	    return true;
	}
	
    } else {

	// Force CM and top-level steps only...

	force_step(b, b->get_system_time());
	force_step(b->get_oldest_daughter(), b->get_system_time());

	// Modify relative velocity of components.

	adjust_vel(b->get_oldest_daughter());

	return true;

    }

    return false;
}

local bool kill_multiples(hdyn * b)
{

    //---------------------------------------------------------------//
    cerr << endl << "***** kill_multiples: not yet implemented...\n\n";
    return false;
    //---------------------------------------------------------------//

    bool killed = false;
    for_all_daughters(hdyn, b, bb)
	if (!bb->is_leaf())
           if (bb->n_leaves() > 2) {
		cerr << "\nFound multiple system " << bb->format_label()
		     << endl;
		print_structure_recursive(bb);
		if (kill_multiple(bb)) {
		    cerr << "...killed\n";
		    killed = true;
		}
	    }

    return killed;
}



// Functions to allow modifications to run-time parameters and options.
// Note the many similarities between modify_params, modify_options,
// and modify_diag -- possible merger of the three?

local void modify_params(hdyn * b, char * name,
			 real& t_end, real& new_dt_log, real& new_dt_snap,
			 int& long_binary_output, char* new_snap_save_file)
{
    // Allow runtime modification of parameters that can be set on
    // the kira command line.

    ifstream file(name);
    if (file) {

	cerr << endl << "***** reading parameter changes from file " << name
	     << "\n\n";

	// File should contain lines of the form:
	//
	//		-keyword value
	//
	// where the keywords a, b, C, d, D, g, k, K, l, N, O, and t are
	// permitted (meanings are the same as if specified on the command
	// line, except that "-t" increments t_end).
	//
	// Parameter changes take place IMMEDIATELY.

	char *s, line[128];

	while (file.get(line, 128, '\n')) {

	    // NOTE: If get() found the terminating '\n', it left it
	    // in the input stream as the first unread character (a C++
	    // feature, not a bug).  To avoid an infinite loop, read
	    // that character now.  Partial reads of lines probably
	    // mean an error, so we simply read characters until the
	    // terminator is found.

	    char c;
	    while (file.get(c) && c != '\n') ;

	    // Parse the line.

	    if (s = strstr(line, "-a")) {

		real tmp = atof(s+2);
		if (tmp > 0 && tmp < 1) {
		    cerr << "setting eta = " << tmp << endl;
		    b->set_eta(tmp);
		    putrq(b->get_log_story(), "eta", tmp);
		}

	    } else if (s = strstr(line, "-b")) {

		int tmp = atoi(s+2);
		if (tmp > 0) {
		    cerr << "setting long_binary_output = " << tmp << endl;
		    long_binary_output = tmp;
		}

	    } else if (s = strstr(line, "-C")) {

		real tmp = atof(s+2);
		if (tmp > 0 && tmp <= 1800) {	// arbitrary upper limit...
		    cerr << "setting grape_max_cpu = " << tmp << endl;
		    b->get_kira_options()->grape_max_cpu = tmp;
		}

	    } else if (s = strstr(line, "-d")) {

		real tmp = atof(s+2);
		if (tmp > 0) {
		    cerr << "setting dt_log = " << tmp << endl;
		    new_dt_log = tmp;
		}

	    } else if (s = strstr(line, "-D")) {

		real tmp = atof(s+2);
		if (tmp > 0) {
		    cerr << "setting dt_snap = " << tmp << endl;
		    new_dt_snap = tmp;
		}

	    } else if (s = strstr(line, "-l")) {

		// This version is a simple toggle, with minimal checks...

		bool f = b->get_p_friction();
		f = !f;
		if (!b->get_plummer()) f = false;
		b->set_p_friction(f);

	    } else if (s = strstr(line, "-g")) {

		real tmp = atof(s+2);
		if (tmp > 0 && tmp < 100) {
		    cerr << "setting lag_factor = " << tmp << " (squared)\n";
		    b->set_lag_factor(tmp*tmp);
		    putrq(b->get_log_story(), "lag_factor", tmp);
		}

	    } else if (s = strstr(line, "-k")) {

		real tmp = atof(s+2);
		if (tmp > 0 && tmp < 1) {
		    cerr << "setting gamma = " << tmp << endl;
		    b->set_gamma2(tmp*tmp);
		    putrq(b->get_log_story(), "gamma", tmp);
		}

	    } else if (s = strstr(line, "-K")) {

		int tmp = atoi(s+2);
		if (tmp > 0) {
		    int kappa_max = (int) pow(2.0, tmp);
		    b->set_max_slow_factor(kappa_max);
		    cerr << "setting maximum slowdown factor = "
			 << kappa_max << endl;
		    putiq(b->get_log_story(), "log_max_slow", tmp);
		}

	    } else if (s = strstr(line, "-l")) {

		real tmp = atof(s+2);
		if (tmp > 0 && tmp < 100) {
		    cerr << "setting d_min_fac = " << tmp << endl;
		    b->set_d_min_fac(tmp);
		    putrq(b->get_log_story(), "d_min_fac", tmp);
		}

	    } else if (s = strstr(line, "-N")) {

		int tmp = atoi(s+2);
		if (tmp > 0) {
		    cerr << "setting n_check = " << tmp << endl;
		    b->get_kira_diag()->n_check_heartbeat = tmp;
		}

	    } else if (s = strstr(line, "-O")) {

		strcpy(new_snap_save_file, s+3);
		if (new_snap_save_file[0] != '\0') {
		    cerr << "setting snap_save_file = "
			 << new_snap_save_file << endl;
		}

	    } else if (s = strstr(line, "-t")) {

		real tmp = atof(s+2);
		t_end += tmp;
		cerr << "setting t_end = " << t_end << endl;

	    } else if (s = strstr(line, "-u")) {

		int tmp = atoi(s+2);
		if (tmp == 0)
		  set_allow_multiples(b, false);
		else
		  set_allow_multiples(b, true);
		print_unperturbed_options(b);

	    } else if (s = strstr(line, "-U")) {

		toggle_unperturbed(b, 0);
		print_unperturbed_options(b);

	    } else if (s = strstr(line, "-x")) {

		b->get_kira_options()->print_xreal
		    = !b->get_kira_options()->print_xreal;
		cerr << "setting print_xreal = "
		     << b->get_kira_options()->print_xreal << endl;

	    }
	}

	file.close();
	clean_up_file(name);

	cerr << endl;
    }
}

local int get_value(char *s, int default_value)		// overloaded!
{
    // Look for an "=" in the line.

    char *ss;
    if (ss = strstr(s, "="))
	return atoi(ss+1);
    else
	return default_value;
}

local real get_value(char *s, real default_value)
{
    char *ss;
    if (ss = strstr(s, "="))
	return atof(ss+1);
    else
	return default_value;
}

local bool get_value(char *s, real& v1, real& v2)
{
    char *ss;
    if (ss = strstr(s, "="))
	if (sscanf(ss+1, "%f %f", &v1, &v2) == 2)
	    return true;

    return false;
}

local char *get_string(char *s)
{
    char *ss;
    if (ss = strstr(s, "=")) {
	ss++;
	while (*ss <= ' ') ss++;
	char *s1 = ss;			// first nonblank character after "="
	char *s2 = s1;
	while (*s2 > ' ') s2++;		// s2 = next blank character after s1

	// Create permanent data...

	int i;
	char *string = new char[s2 - s1 + 1];
	for (i = 0; i < s2-s1; i++) string[i] = *(s1+i);
	string[i] = '\0';

	return string;
    }

    return NULL;
}


local void modify_options(hdyn * b, char * name, bool del = true)
{
    // Modify entries in the hdyn::kira_options class.

    // ***************************************************************
    // ****  Note that options/diag names need to be unique       ****
    // ****  to accommodate use of combined .kira file.           ****
    // ****  Abbreviated names are unhelpful -- reconsider?	  ****
    // ****  Future option names should be long and descriptive.  ****
    // ***************************************************************

    int count = 0;

    ifstream file(name);
    if (file) {

	char *s, line[128];

	while (file.get(line, 128, '\n')) {

	    char c;
	    while (file.get(c) && c != '\n') ;	// see note in modify_params

	    // Use "keyword [= value]" format (if a value is required),
	    // with reasonable abbreviations and alternatives allowed...

	    // Allow blank/comment lines:

	    if (line[0] == '\n' || line[0] == '\0' || line[0] == '#')
		continue;

	    bool print = true;

	    // First option is to change the general Starlab "precision"
	    // setting that controls format of real output from put_node().

	    if (s = strstr(line, "precision"))

		adjust_starlab_precision(get_value(s, -1));  // default = reset
	
	    else if (s = strstr(line, "print"))

		b->get_kira_options()->print();
	
	    else if (s = strstr(line, "reset")) {

		// Restore factory settings...

		delete b->get_kira_options();
		b->set_kira_options(new kira_options);
	

		// *** print_xreal is taken care of in modify_params() ***


	    } else if ((s = strstr(line, "perturber_criterion"))
		       || (s = strstr(line, "perturber_c")))

		b->get_kira_options()->perturber_criterion
		    = get_value(s, DEFAULT_PERTURBER_CRITERION);

	    else if ((s = strstr(line, "optimize_scheduling"))
		       || (s = strstr(line, "optimize_s")))

		b->get_kira_options()->optimize_scheduling
		    = get_value(s, DEFAULT_OPTIMIZE_SCHEDULING);

	    else if ((s = strstr(line, "optimize_block"))
		     || (s = strstr(line, "optimize_b")))

		b->get_kira_options()->optimize_block
		    = get_value(s, DEFAULT_OPTIMIZE_BLOCK);


	    //  ***  allow_unperturbed and allow_multiples are  ***
	    //  ***  taken care of in modify_params()...        ***


	    else if ((s = strstr(line, "min_unpert_steps"))
		     || (s = strstr(line, "min_unpert")))

		b->get_kira_options()->min_unpert_steps
		    = get_value(s, DEFAULT_MIN_UNPERT_STEPS);

	    else if ((s = strstr(line, "full_merge_tol_for_close"))
		     || (s = strstr(line, "full_merge_tol_"))
		     || (s = strstr(line, "full_merge_close")))

		b->get_kira_options()->full_merge_tol_for_close_binary
		    = get_value(s, DEFAULT_FULL_MERGE_TOL_FOR_CLOSE_BINARY);

	    else if ((s = strstr(line, "full_merge_tolerance"))
		     || (s = strstr(line, "full_merge")))

		b->get_kira_options()->full_merge_tolerance
		    = get_value(s, DEFAULT_FULL_MERGE_TOLERANCE);

	    else if ((s = strstr(line, "relax_factor"))
		     || (s = strstr(line, "relax")))

		b->get_kira_options()->relax_factor
		    = get_value(s, DEFAULT_RELAX_FACTOR);

	    else if ((s = strstr(line, "partial_merge_factor"))
		     || (s = strstr(line, "partial_merge")))

		b->get_kira_options()->partial_merge_factor
		    = get_value(s, DEFAULT_PARTIAL_MERGE_FACTOR);

	    else if ((s = strstr(line, "multiple_merge_tolerance"))
		     || (s = strstr(line, "multiple_merge")))

		b->get_kira_options()->multiple_merge_tolerance
		    = get_value(s, DEFAULT_MULTIPLE_MERGE_TOLERANCE);

	    else if ((s = strstr(line, "unconditional_stable_fac"))
		     || (s = strstr(line, "uncond_stable"))
		     || (s = strstr(line, "stable")))

		b->get_kira_options()->unconditional_stable_fac
		    = get_value(s, DEFAULT_UNCONDITIONAL_STABLE_FAC);

	    else if ((s = strstr(line, "partial_stable_fac"))
		     || (s = strstr(line, "partial_stable")))

		b->get_kira_options()->partial_stable_fac
		    = get_value(s, DEFAULT_PARTIAL_STABLE_FAC);

	    else if ((s = strstr(line, "aarseth_stable_fac"))
		     || (s = strstr(line, "aarseth_stable"))
		     || (s = strstr(line, "aarseth_fac")))

		b->get_kira_options()->aarseth_stable_fac
		    = get_value(s, DEFAULT_AARSETH_STABLE_FAC);

	    else if ((s = strstr(line, "use_aarseth_criterion"))
		     || (s = strstr(line, "use_aarseth"))
		     || (s = strstr(line, "aarseth")))

		b->get_kira_options()->use_aarseth_criterion
		    = get_value(s, DEFAULT_USE_AARSETH_CRITERION);

	    else if ((s = strstr(line, "close_criterion"))
		     || (s = strstr(line, "close")))

		b->get_kira_options()->close_criterion
		    = get_value(s, DEFAULT_CLOSE_CRITERION);

	    else if ((s = strstr(line, "allow_keplstep"))
		     || (s = strstr(line, "allow_kepl")))

		b->get_kira_options()->allow_keplstep
		    = get_value(s, DEFAULT_ALLOW_KEPLSTEP);

	    else if ((s = strstr(line, "use_old_correct_acc_and_jerk"))
		     || (s = strstr(line, "use_old_corr")))

		b->get_kira_options()->use_old_correct_acc_and_jerk
		    = get_value(s, DEFAULT_USE_OLD_CORRECT_ACC_AND_JERK);

	    else if ((s = strstr(line, "grape_check_count"))
		     || (s = strstr(line, "grape_check")))

		b->get_kira_options()->grape_check_count
		    = get_value(s, DEFAULT_GRAPE_CHECK_COUNT);

	    else if ((s = strstr(line, "grape_coll_freq"))
		     || (s = strstr(line, "grape_coll")))

		b->get_kira_options()->grape_coll_freq
		    = get_value(s, DEFAULT_GRAPE_COLL_FREQ);

	    else if ((s = strstr(line, "grape_pert_freq"))
		     || (s = strstr(line, "grape_pert")))

		b->get_kira_options()->grape_pert_freq
		    = get_value(s, DEFAULT_GRAPE_PERT_FREQ);

	    else if ((s = strstr(line, "grape_max_cpu"))
		     || (s = strstr(line, "grape_max")))

		b->get_kira_options()->grape_max_cpu
		    = get_value(s, DEFAULT_GRAPE_MAX_CPU);

	    else if ((s = strstr(line, "use_perturbed_list"))
		     || (s = strstr(line, "use_pert")))

		b->get_kira_options()->use_perturbed_list
		    = get_value(s, DEFAULT_USE_PERTURBED_LIST);

	    else if ((s = strstr(line, "max_slow_factor"))
		     || (s = strstr(line, "max_slow_fac")))

		b->set_max_slow_factor(get_value(s, DEFAULT_MAX_SLOW_FACTOR));

	    else if ((s = strstr(line, "max_slow_perturbation"))
		     || (s = strstr(line, "max_slow_pert")))

		b->set_max_slow_perturbation(get_value(s,
					     DEFAULT_MAX_SLOW_PERTURBATION));
	    else

		print = false;

	    if (print) {
		if (count++ == 0) cerr << endl;
		cerr << "read runtime option " << line << endl;
	    }
	}

	file.close();
	if (del) clean_up_file(name);
    }
}


local void modify_diag(hdyn * b, char * name, bool del = true)
{
    // Modify entries in the hdyn::kira_diag class.

    // *************************************************************
    // ****  Note that options/diag names need to be unique     ****
    // ****  to accommodate use of combined .kira file.         ****
    // ****  Abbreviated names are unhelpful -- reconsider?	****
    // ****  Future diag names should be long and descriptive.  ****
    // *************************************************************

    int count = 0;

    ifstream file(name);
    if (file) {

	char *s, line[128];

	while (file.get(line, 128, '\n')) {

	    char c;
	    while (file.get(c) && c != '\n') ;	// see note in modify_params

	    // Use "keyword [= value]" format (if a value is required),
	    // with reasonable abbreviations and alternatives allowed...
	    //
	    // Convention: for true/false switches, 0 = false, 1 = true
	    //		   for levels,   -1 turns off the output,
	    //			       >= 0 sets the value

	    // Allow blank/comment lines:

	    if (line[0] == '\n' || line[0] == '\0' || line[0] == '#')
		continue;

	    bool print = true;

	    if (s = strstr(line, "print"))

		b->get_kira_diag()->print();
	
	    else if (s = strstr(line, "name")) {

		if (s = get_string(line))
		    b->get_kira_diag()->set_name(s);
		else
		    b->get_kira_diag()->clear_name();

	    } else if (s = strstr(line, "dt")) {

		real dt = get_value(s, 0.0);
		if (dt > 0)
		    b->get_kira_diag()->set_range(b->get_system_time(),
						  b->get_system_time() + dt);

	    } else if (s = strstr(line, "range")) {

		real t1, t2;
		if (get_value(line, t1, t2))
		    if (t2 > t1 && t2 > b->get_system_time())
			b->get_kira_diag()->set_range(t1, t2);

	    } else if (s = strstr(line, "reset")) {

		// Restore factory settings...

		delete b->get_kira_diag();
		b->set_kira_diag(new kira_diag);
	
	    } else if (s = strstr(line, "kira_main"))

		b->get_kira_diag()->kira_main = get_value(s, 1);

	    else if (s = strstr(line, "n_check_heartbeat"))

		b->get_kira_diag()->n_check_heartbeat
		    = get_value(s, 50000);

	    else if (s = strstr(line, "check_heartbeat"))

		b->get_kira_diag()->check_heartbeat
		    = get_value(s, 1);

	    else if (s = strstr(line, "n_check_runtime"))

		b->get_kira_diag()->n_check_runtime = get_value(s, 2500);

	    else if ((s = strstr(line, "unpert_function_id"))
		     || (s = strstr(line, "unpert_id")))

		b->get_kira_diag()->unpert_function_id
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_start_unperturbed"))
		     || (s = strstr(line, "report_start_unpert"))
		     || (s = strstr(line, "unpert_start")))

		b->get_kira_diag()->report_start_unperturbed
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_continue_unperturbed"))
		     || (s = strstr(line, "report_cont_unpert"))
		     || (s = strstr(line, "unpert_cont")))

		b->get_kira_diag()->report_continue_unperturbed
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_end_unperturbed"))
		     || (s = strstr(line, "report_end_unpert"))
		     || (s = strstr(line, "unpert_end")))

		b->get_kira_diag()->report_end_unperturbed
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_pericenter_reflection"))
		     || (s = strstr(line, "report_peri_unpert"))
		     || (s = strstr(line, "report_peri"))
		     || (s = strstr(line, "unpert_peri")))

		b->get_kira_diag()->report_pericenter_reflection
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_multiple"))
		     || (s = strstr(line, "report_mult_unpert"))
		     || (s = strstr(line, "unpert_mult")))

		b->get_kira_diag()->report_multiple
		    = get_value(s, 1);

	    else if ((s = strstr(line, "end_unpert_report_level"))
		     || (s = strstr(line, "end_unpert_level")))

		b->get_kira_diag()->end_unpert_report_level
		    = get_value(s, 0);

	    else if ((s = strstr(line, "unpert_report_level"))
		     || (s = strstr(line, "unpert_level")))

		b->get_kira_diag()->unpert_report_level
		    = get_value(s, 0);

	    else if ((s = strstr(line, "multiple_report_level"))
		     || (s = strstr(line, "mult_report_level")))

		b->get_kira_diag()->multiple_report_level
		    = get_value(s, 0);

	    else if ((s = strstr(line, "tree_level"))) {

		b->get_kira_diag()->tree_level
		    = get_value(s, 0);
		if (b->get_kira_diag()->tree_level < 0)
		    b->get_kira_diag()->tree = false;

	    } else if ((s = strstr(line, "tree"))
		       || (s = strstr(line, "diag_tree")))

		b->get_kira_diag()->tree
		    = get_value(s, 1);

	    else if ((s = strstr(line, "ev_function_id")))

		b->get_kira_diag()->ev_function_id
		    = get_value(s, 1);

	    else if ((s = strstr(line, "grape_level"))) {

		b->get_kira_diag()->grape_level
		    = get_value(s, 0);
		if (b->get_kira_diag()->grape_level < 0)
                    b->get_kira_diag()->grape = false;

	    } else if ((s = strstr(line, "grape_ev")))

		b->get_kira_diag()->grape
		    = get_value(s, 1);

	    else if ((s = strstr(line, "timestep_check")))

		b->get_kira_diag()->timestep_check
		    = get_value(s, 1);

	    else if (s = strstr(line, "diag_correct"))

		b->get_kira_diag()->correct
		    = get_value(s, 1);

	    else if (s = strstr(line, "kira_ev"))

		b->get_kira_diag()->kira_ev = get_value(s, 1);

	    else if ((s = strstr(line, "slow_perturbed"))
		     || (s = strstr(line, "slow_pert")))

		b->get_kira_diag()->slow_perturbed
		    = get_value(s, 1);

	    else if ((s = strstr(line, "slow_level"))) {

		b->get_kira_diag()->slow_level
		    = get_value(s, 0);
		if (b->get_kira_diag()->slow_level < 0)
		    b->get_kira_diag()->slow = false;

	    } else if (s = strstr(line, "slow_diag"))

		b->get_kira_diag()->slow
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_stellar_evolution"))
		       || (s = strstr(line, "report_stellar_e"))
		     || (s = strstr(line, "stellar_e"))) {

		b->get_kira_diag()->report_stellar_evolution
		    = get_value(s, 1);

//		PRL(b->get_kira_diag()->report_stellar_evolution);

	    } else if ((s = strstr(line, "report_stellar_mass_loss"))
		       || (s = strstr(line, "report_stellar_m"))
		       || (s = strstr(line, "stellar_m")))

		b->get_kira_diag()->report_stellar_mass_loss
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_binary_mass_loss"))
		       || (s = strstr(line, "report_b"))
		       || (s = strstr(line, "binary")))

		b->get_kira_diag()->report_binary_mass_loss
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_adjust_perturbed"))
		       || (s = strstr(line, "report_adj"))
		       || (s = strstr(line, "adjust_pert")))

		b->get_kira_diag()->report_adjust_perturbed_list
		    = get_value(s, 1);

	    else if ((s = strstr(line, "report_correct_perturber_list"))
		       || (s = strstr(line, "report_corr"))
		       || (s = strstr(line, "correct_pert")))

		b->get_kira_diag()->report_correct_perturber_list
		    = get_value(s, 1);

	    else if (s = strstr(line, "diag_ev"))    // poor choice of keyword!

		b->get_kira_diag()->ev
		    = get_value(s, 1);

	    else

		print = false;

	    if (print) {
		if (count++ == 0) cerr << endl;
		cerr << "read runtime diag option " << line << endl;
	    } else if (strstr(line, "=")) {
		if (count++ == 0) cerr << endl;
		cerr << "ignored unknown diag option \"" << line << "\""
		     << endl;
	    }
	}

	file.close();
	if (del) clean_up_file(name);
    }
}


local void pp(hdyn * b, char * name)
{
    ifstream file(name);
    if (file) {

	cerr << endl << "***** Running pp3 on particle(s) listed in file "
	     << name << ":\n\n";

	// File should contain valid particle names, one per line.

	char *s, line[128];

	while (file.get(line, 128, '\n')) {

	    // NOTE: If get() found the terminating '\n', it left it
	    // in the input stream as the first unread character (a C++
	    // feature, not a bug).  To avoid an infinite loop, read
	    // that character now.  Partial reads of lines probably
	    // mean an error, so we simply read characters until the
	    // terminator is found.

	    char c;
	    while (file.get(c) && c != '\n') ;

	    // Add a trailing '\0' to the input line at the first
	    // non-printing location.

	    bool ok = false;
	    for (int i = 0; i < 128; i++) {
		if (line[i] < ' ') {
		    line[i] = '\0';
		    ok = true;
		    break;
		}
	    }

	    // (Silently ignore unterminated lines.)

	    if (ok) {

		// Parse the line.

		hdyn* bb = (hdyn*) node_with_name(line, b->get_root());

		if (!bb)
		    cerr << "Unable to find node " << line << endl;
		else
		    pp3(bb);
	    }
	}

	file.close();
	clean_up_file(name);

	cerr << endl;
    }
}

local void qt(hdyn * b, char * name)
{
    ifstream file(name);
    if (file) {
	cerr << endl << "kira_runtime: launching Qt browser at time "
	     << b->get_system_time();
	char display[1024];
	int ret;
	if (file.get(display, 1023, '\n')) {
	    cerr << " with X display " << display << endl;
	    ret = Qt_pp3(b, display);
	} else {
	    cerr << " with default X display" << endl;
	    ret = Qt_pp3(b);
	}
	file.close();
	clean_up_file(name);
	cerr << "...continuing with return value " << ret << endl;
    }
}

local void dump_to_file(hdyn* b, char* name)
{
    ifstream file(name);
    bool written = false;
    if (file) {
	char dumpfile[128];
	if (file.get(dumpfile, 128, '\n')) {
	    ofstream dump(dumpfile);
	    if (dump) {
		put_node(b, dump, b->get_kira_options()->print_xreal);
		dump.close();
		cerr << "Data written to file " << dumpfile << endl;
		written = true;
	    }
	}
	file.close();
	clean_up_file(name);
    }
    if (!written)
	cerr << "No data written" << endl;
}


// Check_kira_init:  read and process configuration files if present,
//		     but do *not* delete them.

void check_kira_init(hdyn *b)
{
    // Look in standard places for nonstandard options.  Merge the old
    // old "INIT_DIAG" and "INIT_OPTIONS" settings into a single .kira
    // file (change made by Steve, 1/02).

    char file[1024];

    // Look first in $HOME/.kira ...
    // Note: may not exist if (e.g.) we are running kira in a cgi script.

    if (getenv("HOME")) {

	strcpy(file, getenv("HOME"));
	strcat(file, "/.kira");

	if (check_file(file, false)) {
	    cerr << endl
		 << "reading init file " << getenv("HOME") << "/.kira"
		 << endl << flush;
	    modify_diag(b, file, false);
	    modify_options(b, file, false);
	}
    }

    // ... then allow local setup in ./.kira.

    strcpy(file, "./.kira");

    if (check_file(file, false)) {
	cerr << endl
	     << "reading init file " << "./.kira"
	     << endl << flush;
	modify_diag(b, file, false);
	modify_options(b, file, false);
    }
}


// Check_kira_runtime:  read and process reconfiguration files if present,
//			then delete them.

bool check_kira_runtime(hdyn* b,
			real& t_end, real& new_dt_log, real& new_dt_snap,
			int& long_binary_output, char* new_snap_save_file,
			bool& tree_changed)
{
    static bool allow_check = true;
    bool status = false;

    if (check_file("ALLOW_CHECK")) {
	allow_check = !allow_check;
	cerr << endl
	     << "**** setting allow_check = " << allow_check
	     << endl;
    }

    if (!allow_check) return false;

    // Note the order in which we check things -- placing ABORT
    // at end allows us to get other output before quitting.

    if (check_file("PARAMS", false))
	modify_params(b, "PARAMS", t_end,
		      new_dt_log, new_dt_snap,
		      long_binary_output, new_snap_save_file);

    if (check_file("DIAG", false))
	modify_diag(b, "DIAG");

    if (check_file("OPTIONS", false))
	modify_options(b, "OPTIONS");

    if (check_file("PP", false))
	pp(b, "PP");

    if (check_file("QT", false))
	qt(b, "QT");

    if (check_file("DUMP", false)) {
	cerr << endl
	     << "***** (Unsynchronized system dump forced by user)\n\n";
	cerr << "Time = " << b->get_system_time()
	     << "  N = "    << b->n_leaves()
	     << "  mass = " << total_mass(b) << endl;
	dump_to_file(b, "dump");
	cerr << endl;
    }

    if (check_file("ENERGY")) {
	cerr << endl << "***** (Energy output forced by user)\n\n";
	cerr << "Time = " << b->get_system_time()
	     << "  N = "    << b->n_leaves()
	     << "  mass = " << total_mass(b) << endl;

	// Note that, as presently coded, the computation of the energy
	// also recomputes acc and jerk for all nodes, and therefore
	// affects the future evolution, making runs non-reproducible.

	// pp3(b);
	print_recalculated_energies(b);
	// pp3(b);
	cerr << endl;
    }

    if (check_file("LOG")) {
	cerr << endl << "***** (Log output forced by user)\n\n";
	cerr << "Time = " << b->get_system_time()
	     << "  N = "    << b->n_leaves()
	     << "  mass = " << total_mass(b) << endl;

	// print_recalculated_energies(b);	// omit because this can
						// change the evolution...
						// -- doesn't affect evolution
						// any more, but omit anyway

	print_counters(b->get_kira_counters());

	set_block_check(false);			// avoid annoying
	print_statistics(b);			// "mismatch" messages
	set_block_check(true);

	cerr << endl;
    }

    if (check_file("MULTIPLES")) {
	cerr << endl << "***** (Multiple output forced by user)\n\n";
	cerr << "Time = " << b->get_system_time()
	     << "  N = "    << b->n_leaves()
	     << "  mass = " << total_mass(b) << endl;
	print_multiples(b);
	cerr << endl;
    }

    if (check_file("KILL_MULTIPLES"))
	tree_changed |= kill_multiples(b);

    if (check_file("SCHED")) {
	cerr << endl << "***** (Scheduling output forced by user)\n\n";
	cerr << "Time = " << b->get_system_time()
	     << "  N = "    << b->n_leaves()
	     << "  mass = " << total_mass(b) << endl;
	print_counters(b->get_kira_counters());

	cerr << "Particles just advanced:" << endl;
	for_all_nodes(hdyn, b, bb)
	    if (bb->get_time() == bb->get_system_time())
		cerr << "    " << bb->format_label()
		     << "  timestep = " << bb->get_timestep()
		     << endl;

	cerr << endl;
    }

    // Check this last...

    if (check_file("ABORT")) {

	cerr << endl << "***** Abort forced by user at time "
	     << b->get_system_time() << "\n\n" << flush;
	clean_up_files();

	status = true;
    }

    return status;	// "true" ==> quit kira on return
}
