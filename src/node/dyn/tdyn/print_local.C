local void print_events(segment *s, real t)
{
    PRI(4); cerr << "events for s = " << s << ":" << endl;

    tdyn *bn = s->get_first_event();
    PRI(4); cerr << "base node "; PRC(bn); PRC(bn->get_time());
    PRL(bn->format_label());

    tdyn *b = bn;
    while (b->get_time() < t) {
	PRI(8); PRC(b); PRL(b->get_time());
	if (b->get_next()) b = b->get_next();
    }
    PRI(8); PRC(b); PRC(b->get_time());
    if (b->get_next()) {
	PRL(b->get_next()->get_time());
    } else {
	cerr << "next = NULL" << endl << endl;
    }
}

local void print_details(worldbundle *wb, tdyn *p, real t)
{
    // Re-locate node p at time t and print out relevent information
    // on the local worldline/segment/event structure.

    cerr << "details..." << endl;

    worldline *w = wb->find_worldline(p);
    segment *s = w->get_first_segment();
    segment *sprev = NULL;

    PRI(4); PRL(p->format_label());
    PRI(4); PRC(w); PRL(s);

    while (s->get_t_end() < t) {
	sprev = s;
	s = s->get_next();
    }

    PRI(4); PRL(s);

    PRI(4); PRC(t); PRC(s->get_t_start()); PRL(s->get_t_end());
    if (sprev) {
	PRI(4); PRC(sprev->get_t_start()); PRL(sprev->get_t_end());
	PRI(4); PRL(sprev->get_first_event());
    }

    print_events(s, t);
}
