
// hydro_leapfrog.C

// Replaced all dyn references by dyn pointers... (Steve, 5/03).

#include "dyn.h"
#include "hydro.h"

#ifdef TOOLBOX

local void accumulate_acceleration(dyn *bj,           // target body
				   dyn *bi,           // body whose force is required
				   real eps_squared,  // softening length squared
				   dyn *& bm)         // if non-null, merger has
						      // happened between bi and bm
{
    bm = NULL;
    if (bj->get_oldest_daughter()) {
	for_all_daughters(dyn, bj, bb) {
	    accumulate_acceleration(bb, bi, eps_squared, bm);
    	    if (bm) return;            // merger has happened, go back to merge
	}                              // and recalculate all accelerations
    } else {
	if (bi != bj) {
	    vector d_pos = bi->get_pos() - bj->get_pos();
	    real d_pos_squared = d_pos * d_pos;
	    real sum_of_stellar_r_eff =
		    get_effective_radius(bi) + get_effective_radius(bj);
	    if (sum_of_stellar_r_eff * sum_of_stellar_r_eff > d_pos_squared) {
		bm = bj;
		return;                // merger has happened, go back to merge
	    }                          // and recalculate all accelerations

	    real soft_d_pos_squared = d_pos_squared + eps_squared;

	    real inverse_d_pos_cubed =
		1 / (soft_d_pos_squared * sqrt(soft_d_pos_squared));

	    bi->inc_acc(-inverse_d_pos_cubed * bj->get_mass() * d_pos);
	}
    }
}


local void calculate_acceleration(dyn *bj,
				  dyn *b,            // n-body system pointer
				  real eps_squared,  // softening length squared
				  dyn *& bm1,        // if non-null, merger has
				  dyn *& bm2)        // happened between bm1 and bm2
{
    bm1 = bm2 = NULL;
    if (b->get_oldest_daughter()) {
	for_all_daughters(dyn, b, bb) {
	    calculate_acceleration(bj, bb, eps_squared, bm1, bm2);
	    if (bm1) return;
	}
    } else {
	b->clear_acc();
	dyn *bm;                             // points to potential merger
	accumulate_acceleration(bj, b, eps_squared, bm);
	if (bm) {
	    bm1 = b;
	    bm2 = bm;
	}
    }
}

local void predict_step(dyn *b,		// n-body system pointer
			real dt)	// timestep
{
    if (b->get_oldest_daughter()) {
	for_all_daughters(dyn, b, bb)
	    predict_step(bb, dt);
    } else {
	b->inc_vel(0.5 * dt * b->get_acc());
	b->inc_pos(dt * b->get_vel());
    }
}

local void old_correct_step(dyn *b,	// n-body system pointer
			    real dt)	// timestep
{

    // Full recursive implementation.
    
    if (b->get_oldest_daughter()) {
	old_correct_step(b->get_oldest_daughter(), dt);
    } else {
	b->inc_vel( 0.5 * dt * b->get_acc() );
    }
    if (b->get_younger_sister())
	old_correct_step(b->get_younger_sister(), dt);
}

local void correct_step(dyn *b,		// n-body system pointer
			real dt)	// timestep
{

//  *** Non full-recursive implementation...
    
    if (b->get_oldest_daughter()) {
	for_all_daughters(dyn, b, bb)
	    correct_step(bb,dt);
    } else {
	b->inc_vel(0.5 * dt * b->get_acc());
    }
}

// delete_node_from_general_tree
// delete node b. Do not check if the parent has more than 2 remaining
// daughters or not.

local void delete_node_from_general_tree(dyn *b)
{
    if (b) {
	cerr << "Warning: delete_node_from_general_tree, b is NULL\n";
	return;
    }

    dyn *parent = b->get_parent();
    
    // check if b is head without parent or sisters

    if (parent == NULL) {
	cerr << "Warning: delete_node_from_general_tree, b has no parent\n";
	return;
    }

    b->set_parent(NULL);

    dyn *elder_sister = b->get_elder_sister();
    dyn *younger_sister = b->get_younger_sister();

    if (parent->get_oldest_daughter() == b) {
	parent->set_oldest_daughter(younger_sister);
    }

    if (elder_sister) {
	elder_sister->set_younger_sister(younger_sister);
    }
    if (younger_sister) {
	younger_sister->set_elder_sister(elder_sister);
    }
}

// add_node
// insert b into the tree as the oldest_daughter of the
// parent. The ex-oldest node of the parent becomes the
// younger sister of b.

local void add_node(dyn *b, dyn *parent)
{
    if (b == NULL) {
	cerr << "Warning: add_node, b is NULL\n";
	return;
    }

    if (parent == NULL) {
	cerr << "Warning: add_node, parent is NULL\n";
	return;
    }

    // Set the pointers of ex-oldest-daughter of parent.

    dyn *ex = parent->get_oldest_daughter();
    if (ex)
	ex->set_elder_sister(b);

    parent->set_oldest_daughter(b);

    b->set_elder_sister(NULL);
    b->set_younger_sister(ex);
    b->set_parent(parent);
}

local void merge(dyn *bm1,         // stellar radius proportional to mass
		 dyn *bm2,         // i.e.  R_bm = R_bm1 + R_bm2
		 real t)           // time
{
    dyn *bm = new dyn();

    cerr << "entering merge(" << bm1->get_index() << ", " << bm2->get_index()
	 << ") with r1 = " << bm1->get_pos() << "\n                     "
         << " and r2 = " << bm2->get_pos()   << "\n                     "
	 << "   at t = " << t << endl;

    real  m1 = bm1->get_mass();
    real  m2 = bm2->get_mass();
    real  m_tot = m1 + m2;
    
    bm->set_mass(m_tot);
    bm->set_pos((m1*bm1->get_pos() + m2*bm2->get_pos())/m_tot);
    bm->set_vel((m1*bm1->get_vel() + m2*bm2->get_vel())/m_tot);

    real  new_r_eff = get_effective_radius(bm1) + get_effective_radius(bm2);
    
    addhydro(bm, new_r_eff);
    
    dyn *root = bm1->get_parent(); 

    delete_node_from_general_tree(bm1);
    delete_node_from_general_tree(bm2);
  
    add_node(bm, root);
}

local void calculate_acc_and_merge(dyn *bj,
				   dyn *b,             // n-body system pointer
				   real eps_squared,   // softening length squared
				   real t)             // time
{
    dyn *bm1;           // if non-null, merger has
    dyn *bm2;           // happened between bm1 and bm2

    calculate_acceleration(bj, b, eps_squared, bm1, bm2);

    while (bm1) {
	merge(bm1, bm2, t);
        calculate_acceleration(bj, b, eps_squared, bm1, bm2); // recalculate acc
    }
}

local void step(real &t,        // time
		dyn  *b,        // dyn array                   
		real dt,        // time step of the integration 
		real eps)       // softening length             
{
    t += dt;
    
    predict_step(b,dt);
    calculate_acc_and_merge(b, b, eps*eps, t);
    correct_step(b, dt);
}

local void evolve(real &t,        // time                         
		  dyn *b,         // dyn array                   
		  real delta_t,   // time span of the integration 
		  real dt,        // time step of the integration 
		  real dt_out,    // output time interval
		  real dt_snap,   // snapshot output interval
		  real eps)       // softening length             
{
    real t_end = t + delta_t;      // final time, at the end of the integration
    real t_out = t + dt_out;       // time of next diagnostic output
    real t_snap = t + dt_snap;     // time of next snapshot;
    int  steps = 0;
 
    calculate_acc_and_merge(b, b, eps*eps, t);
    
    while (t < t_end){
	step(t, b, dt, eps);
	steps++;

//      Check for (trivial) output to cerr.

	if (t >= t_out) {
	    cerr << "Time = " << t << "  steps = " << steps << endl;
	    t_out += dt_out;
	}

//      Output a snapshot at the scheduled time or at the end of the run.

	if (t >= t_snap || t >= t_end) {
	    cerr << "time = " << t << endl;
	    put_node(b);
	    cout << flush;
	    t_snap += dt_snap;
	}
    }
}

main(int argc, char **argv)
{
    dyn*  b;             // pointer to the nbody system
    real  t = 0;         // time

    real  delta_t = 1;   // time span of the integration
    real  dt = 0.02;     // time step of the integration
    real  dt_out;        // diagnostic output time interval
    real  dt_snap;       // snap output interval
    real  eps = 0.05;    // softening length 	       	   
    char  *comment;

    extern char *poptarg;
    int  pgetopt(int, char **, char *), c;

    bool  a_flag = FALSE;
    bool  c_flag = FALSE;
    bool  d_flag = FALSE;
    bool  D_flag = FALSE;
    bool  e_flag = FALSE;
    bool  q_flag = FALSE;
    bool  t_flag = FALSE;
    
    while ((c = pgetopt(argc, argv, "a:c:d:D:e:qt:")) != -1)
	switch(c) {
	    case 'a': a_flag = TRUE;
		      dt = atof(poptarg);
		      break;
	    case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
	    case 'd': d_flag = TRUE;
		      dt_out = atof(poptarg);
		      break;
	    case 'D': D_flag = TRUE;
		      dt_snap = atof(poptarg);
		      break;
	    case 'e': e_flag = TRUE;
		      eps = atof(poptarg);
		      break;
	    case 'q': q_flag = TRUE;
		      break;
	    case 't': t_flag = TRUE;
		      delta_t = atof(poptarg);
		      break;
            case '?': cerr <<
		      "usage: hydro_leapfrog -t # -a # -e # -D # -d # " <<
		      "[-c \"..\"]\n" <<
		      "for t (time span), a (time step length), " <<
		      "d (output interval) D (snapshot output interval) " <<
		      "and e (softening length)\n";
		      exit(1);
	}            

    if (!q_flag) {

	// Check input arguments and echo defaults.

	if (!t_flag) cerr << "default delta_t = " << delta_t << "\n";
	if (!a_flag) cerr << "default dt = " << dt << "\n";
	if (!d_flag) cerr << "default dt_out = " << dt_out << "\n";
	if (!e_flag) cerr << "default eps = " << eps << "\n";
    }

    if (!D_flag) dt_snap = delta_t;

    b = get_dyn(cin, new_hydro);
    
    if (c_flag == TRUE) b->log_comment(comment);
    b->log_history(argc, argv);

    evolve(t, b, delta_t, dt, dt_out, dt_snap, eps);
}

#endif
