local void print_binary_diagnostics(real t, tdyn *b, tdyn *bb, pdyn *curr)
{
    PRC(t); PRL(b->format_label());
    PRL(bb->get_time());
    PRL(b->get_time());

    tdyn *n = b->get_next();
    if (n) PRL(n->get_time());

    real fac = 1 + bb->get_mass()
		/ bb->get_binary_sister()->get_mass();
    real M = bb->get_parent()->get_mass();

    real rbb = abs(fac*bb->get_pos());
    real vbb = square(fac*bb->get_vel());
    real ebb = 0.5*vbb - M / rbb;

    real rb = abs(fac*b->get_pos());
    real vb = square(fac*b->get_vel());
    real eb = 0.5*vb - M / rb;

    real rc = abs(fac*curr->get_pos());
    real vc = square(fac*curr->get_vel());
    real ec = 0.5*vc - M / rc;

    PRI(4); PRC(rbb); PRC(vbb); PRL(ebb);
    PRI(4); PRC(rb); PRC(vbb); PRL(eb);
    PRI(4); PRC(rc); PRC(vc); PRL(ec);

    if (n) {
	real rn = abs(fac*n->get_pos());
	real vn = square(fac*n->get_vel());
	real en = 0.5*vn - M / rn;
	PRI(4); PRC(rn); PRC(vn); PRL(en);
    }

    PRI(4); PRL(-M/(2*ebb));

#ifndef NEW_INTERP

    // Old code:

    PRL(interpolate_pos(b, b->get_time(), bb));
    PRL(interpolate_pos(b, t, bb));
    if (n)
	PRL(interpolate_pos(b, n->get_time(), bb));

#else

    // Equivalent new code uses set_interpolated_pos():

    vec interp_pos_time;
    set_interpolated_pos(b, b->get_time(), interp_pos_time, bb);
    PRL(interp_pos_time);
    vec interp_pos_t;
    set_interpolated_pos(b, t, interp_pos_t, bb);
    PRL(interp_pos_t);
    if (n) {
	vec interp_pos_ntime;
	set_interpolated_pos(b, n->get_time(), interp_pos_ntime, bb);
	PRL(interp_pos_ntime);
    }

#endif

    real dt = t - (real)b->get_time();
    real dtn = 0;
    if (n) dtn = n->get_time() - b->get_time();
    cerr << "interpolation..." << endl;
    PRC(dt); PRL(dtn);
    PRL(b->get_pos());
    PRL(dt*b->get_vel());
    PRL(dt*dt*b->get_acc());
    PRL(dt*dt*dt*b->get_jerk());
}
