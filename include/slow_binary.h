
/// @file slow_binary.h  Slowdown of lightly perturbed binaries.
//			 (Based on "slow-KS" treatment of Mikkola & Aarseth.)

/// \a slow_binary:  Class handling slowdown of lightly perturbed binaries.

class slow_binary {

    private:

	int  kappa;		///< Slowdown factor (power of 2).
	xreal t_init;		///< Time when this kappa value began.
	xreal t_apo;		///< Time of last "apocenter" (convenience).
	real tau;		///< Slow time.
	real tau_pred;		///< Slow t_pred.
	real dtau;		///< Slow timestep: true timestep = dtau*kappa.
	vec acc_p;		///< Perturbation to CM acceleration.
	vec jerk_p;		///< Perturbation to CM jerk.
	vec old_acc_p;		///< Previous acc_p.
	vec old_jerk_p;		///< Previous jerk_p.
	bool stop;		///< Flag to force stop next time around.

    public:

	slow_binary(int k = 1) {
	    kappa = k;
	    t_init = t_apo = tau = tau_pred = dtau = 0;
	    acc_p = jerk_p = old_acc_p = old_jerk_p = 0;
	    stop = false;
	}

	// Convert slow time to real time.

	xreal tau_to_time(real tt = 0) {
	    if (tt == 0)
		return t_init + kappa * tau;
	    else
		return t_init + kappa * tt;
	}

	/// Convert time to slow time.

	real time_to_tau(xreal t)	{return ((real)(t - t_init)) / kappa;}

	void set_kappa(int k)		{kappa = k;}
	int  get_kappa()		{return kappa;}

	void set_t_init(xreal t)	{t_init = t;}
	xreal get_t_init()		{return t_init;}

	void set_t_apo(xreal t)		{t_apo = t;}
	xreal get_t_apo()		{return t_apo;}

	void set_tau(real t)		{tau = t;}
	void inc_tau(real dt)		{tau += dt;}
	real get_tau()			{return tau;}

	void set_tau_pred(real t)	{tau_pred = t;}
	void clear_tau_pred()		{tau_pred = -VERY_LARGE_NUMBER;}
	void init_tau_pred()		{tau_pred = tau;}

	real get_tau_pred()		{return tau_pred;}
	
	void set_dtau(real dt)		{dtau = dt;}
	real get_dtau()			{return dtau;}

	void set_acc_p(vec a)		{acc_p = a;}
	vec get_acc_p()			{return acc_p;}

	void set_jerk_p(vec j)		{jerk_p = j;}
	vec get_jerk_p()		{return jerk_p;}

	void set_old_acc_p(vec a)	{old_acc_p = a;}
	vec get_old_acc_p()		{return old_acc_p;}

	void set_old_jerk_p(vec j)	{old_jerk_p = j;}
	vec get_old_jerk_p()		{return old_jerk_p;}

	void store_old_force() {
	    old_acc_p = acc_p;
	    old_jerk_p = jerk_p;
	}

	void set_stop(bool v = true)	{stop = v;}
	bool get_stop()			{return stop;}

};

class _dyn_;			// to permit the_node pointer below...

/// \a slow_perturbed:  Class handling for slow binaries perturbed by this node.

class slow_perturbed {

    private:

	_dyn_ *the_node;	///< CM node (hdyn, really) of a slow binary
	int kappa;		///< expected slowdown factor
	vec acc_p;		///< perturbation to CM acceleration
	vec jerk_p;		///< perturbation to CM jerk
	vec old_acc_p;		///< previous acc_p
	vec old_jerk_p;		///< previous jerk_p
	slow_perturbed *sp;	///< linked list pointer

    public:

	slow_perturbed() {
	    the_node = NULL;
	    kappa = 1;
	    acc_p = jerk_p = old_acc_p = old_jerk_p = 0;
	    sp = NULL;
	}

	~slow_perturbed() {
	    if (sp) delete sp;		// recursive
	}

	void set_node(_dyn_ *n)		{the_node = n;}
	_dyn_ *get_node()		{return the_node;}

	void set_kappa(int k)		{kappa = k;}
	int  get_kappa()		{return kappa;}

	void set_acc_p(vec a)		{acc_p = a;}
	vec get_acc_p()			{return acc_p;}

	void set_jerk_p(vec j)		{jerk_p = j;}
	vec get_jerk_p()		{return jerk_p;}

	void set_old_acc_p(vec a)	{old_acc_p = a;}
	vec get_old_acc_p()		{return old_acc_p;}

	void set_old_jerk_p(vec j)	{old_jerk_p = j;}
	vec get_old_jerk_p()		{return old_jerk_p;}

	void store_old_force() {
	    old_acc_p = acc_p;
	    old_jerk_p = jerk_p;
	    if (sp) sp->store_old_force();	// recursive!
	}

	void set_next(slow_perturbed *s)	{sp = s;}
	slow_perturbed *get_next()		{return sp;}
};
