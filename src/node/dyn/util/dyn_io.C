
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// dyn_io:  test Starlab dyn class I/O functions.
////
////          Define scan_dyn_story and print_dyn_story for
////          the dyn class.
////
//// Options: none

#include "dyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize all static dyn data here...

xreal dyn::system_time          = 0.0;
real dyn::real_system_time      = 0.0;
bool dyn::use_sstar	        = false;

bool dyn::ignore_internal	= false;

unsigned int  dyn::external_field	= 0;

int  dyn::tidal_type	= 0;
real dyn::omega	= 0;
real dyn::omega_sq	= 0;
real dyn::alpha1	= 0;
real dyn::alpha3	= 0;
vector dyn::tidal_center = vector(0,0,0);

real dyn::p_mass = 0;
real dyn::p_scale_sq = 0;
vector dyn::p_center = vector(0,0,0);

real dyn::pl_coeff = 0;
real dyn::pl_scale_sq = 0;
real dyn::pl_exponent = 0;
vector dyn::pl_center = vector(0,0,0);
real dyn::pl_cutoff = 0;
real dyn::pl_cutoff_sq = 0;
real dyn::pl_mass = 0;
real dyn::pl_softening_sq = 0;

FILE* dyn::ifp = 0;
FILE* dyn::ofp = 0;

bool dyn::col_output = false;

void dyn::print_static(ostream& s)		// default = cerr
{
    node::print_static(s);

    s << "system_time = " << system_time << endl;
    s << "real_system_time = " << real_system_time << endl;

    s << "use_sstar = " << use_sstar << endl;
    s << "ignore_internal = " << ignore_internal << endl;

    s << "external_field = " << external_field << endl;
    s << "tidal_type = " << tidal_type << endl;

    s << "omega = " << omega << endl;
    s << "omega_sq = " << omega_sq << endl;
    s << "alpha1 = " << alpha1 << endl;
    s << "alpha3 = " << alpha3 << endl;
    s << "tidal_center = " << tidal_center << endl;

    s << "p_mass = " << p_mass << endl;
    s << "p_scale_sq = " << p_scale_sq << endl;
    s << "p_center = " << p_center << endl;

    s << "pl_coeff = " << pl_coeff << endl;
    s << "pl_scale_sq = " << pl_scale_sq << endl;
    s << "pl_exponent = " << pl_exponent << endl;
    s << "pl_center = " << pl_center << endl;

    s << "pl_cutoff = " << pl_cutoff << endl;
    s << "pl_cutoff_sq = " << pl_cutoff_sq << endl;
    s << "pl_mass = " << pl_mass << endl;
    s << "pl_softening_sq = " << pl_softening_sq << endl;
}

static bool read_xreal = false;

istream & dyn::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s,input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

    	if (!strcmp("real_system_time", keyword)) {

	    read_xreal = true;
	    last_real = true;

	    // We don't know anything about parent nodes yet, so it is
	    // not easy to know if we are the root node.  Rule: if we
	    // find real_system_time, assume that we should read an
	    // xreal as system_time.  Otherwise, read it as real.
	    // The tortuous logic is to keep the determination of
	    // which choice we should make completely local.
	    //
	    // Unfortunately, this logic must be duplicated in all
	    // other *dyn::scan_dyn_story functions (see _dyn_io.C,
	    // hdyn_io.C, sdyn3_io.C)...

	} else if (!strcmp("system_time", keyword)) {

	    // Check input format before reading.

	    if (!last_real) read_xreal = false;

	    if (read_xreal) {

		//cerr << "dyn::scan_dyn_story: input "
		//     << "time data type is xreal"
		//     << endl;

		// The following should set real_system_time too...

		set_system_time(get_xreal_from_input_line(input_line));

	    } else {

		//cerr << "dyn::scan_dyn_story: input "
		//     << "time data type is real"
		//     << endl;

		real_system_time = system_time = strtod(val, NULL);
	    }
	} else {

	    last_real = false;

	    if (!strcmp("m", keyword))
		mass = strtod(val, NULL);
	    else if (!strcmp("r", keyword))
		set_vector_from_input_line(pos, input_line);
	    else if (!strcmp("v", keyword))
		set_vector_from_input_line(vel, input_line);
	    else
		add_story_line(dyn_story, input_line);
	}
    }

    return s;
}

ostream& dyn::print_dyn_story(ostream& s,
			      bool print_xreal,		// default = true
			      int short_output)	// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    // Print system time first (root node only).

    if (!parent) {

#ifdef USE_XREAL
	if (print_xreal) {

	    // Note (Steve 5/00): system_time is now xreal and hard to read.
	    // For convenience, also print out a real version of the time.
	    // By printing out real_system_time first, we set a flag that
	    // allows scan_dyn_story to read xreal, rather than real, input.

	    put_real_number(s, "  real_system_time  =  ", (real)system_time);
	    put_real_number(s, "  system_time  =  ", system_time);

	} else

	    put_real_number(s, "  system_time  =  ", (real)system_time);
#else

	put_real_number(s, "  system_time  =  ", system_time);

#endif
    }

    // Mass is now printed by node::print_dyn_story().

    node::print_dyn_story(s, print_xreal, short_output);

    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);

    return s;
}

// gets called by get_dyn when input format is columns of numbers.
dyn* get_col(istream& s,
	     npfp the_npfp,		// note: return value is always node*
	     hbpfp the_hbpfp,
	     sbpfp the_sbpfp,
	     bool use_stories)
{

  static int t;
  dyn* root = (dyn*)the_npfp(the_hbpfp, the_sbpfp, use_stories);
  root->set_index(0);
  // Write time to both the root node (static class data) and to
  // the root node (dyn story), because different user functions
  // may expect either placement.  (To be cleaned up...)
  root->set_system_time(t++);
  // Is time supposed to be incremented here?
  if (use_stories) putrq(root->get_dyn_story(), "t", t);

  static unsigned lineno = 0;
  bool first_dyn = true;
  dyn* bo;

  dyn::set_col_output(true);

  while (true) {
    char line[256];
    if (dyn::get_ifp()) {
      if (!fgets(line, sizeof line, dyn::get_ifp())) break;
    } else {
      if (s.getline(line, sizeof line - 1)) strcat(line, "\n");
      else break;
    }
    switch (++lineno, line[0]) {
    case '\n': return root;
    case '#': continue;
    case ';':
      if (use_stories) {	      // deal with Log (;) and Dyn (;;) stories
	line[strlen(line)-1] = '\0';
	if (line[1] == ';') {
	  if (strlen(&line[2])) {
	      story *st = root->get_dyn_story();
	      if (st)
		add_story_line(st, line+2);
	  }
	} else
	  if (strlen(&line[1])) root->log_comment(&line[1]);
      }
      continue;
    }
    int id;
    double m, x[3], v[3];
    if (sscanf(line, "%i %lg %lg %lg %lg %lg %lg %lg\n", &id, &m, &x[0], &x[1],
               &x[2], &v[0], &v[1], &v[2]) != 8)
      cerr << "Malformed input on line #" << lineno << '\n', exit(1);
    dyn* const b = (dyn*)the_npfp(the_hbpfp, the_sbpfp, use_stories);
    b->set_index(id);
    b->set_mass(m);
    b->set_pos(vector(x[0], x[1], x[2]));
    b->set_vel(vector(v[0], v[1], v[2]));
    b->set_parent(root);
    if (first_dyn) root->set_oldest_daughter(b), first_dyn = false;
    else  b->set_elder_sister(bo), bo->set_younger_sister(b);
    bo = b;
  }

  if (!first_dyn) return root;
  else { delete root; return NULL; }
}

// this I/O problem is getting to be a real pain in the (my) ass

void put_col(dyn* root, ostream& s) {

  story *st;		// note elegant replicated code below!

  if (dyn::get_ofp()) {

    // Special treatment to preserve the root Log and Dyn stories.

    st = root->get_log_story();
    if (st) put_simple_story_contents(dyn::get_ofp(), *st, ";");

    st = root->get_dyn_story();
    if (st) put_simple_story_contents(dyn::get_ofp(), *st, ";;");

    for_all_daughters(dyn, root, i)
      fprintf(dyn::get_ofp(), "%i %g %g %g %g %g %g %g\n", i->get_index(),
	      i->get_mass(), i->get_pos()[0], i->get_pos()[1], i->get_pos()[2],
	      i->get_vel()[0], i->get_vel()[1], i->get_vel()[2]);
    putc('\n', dyn::get_ofp());

  } else {

    // Special treatment to preserve the root Log and Dyn stories.

    st = root->get_log_story();
    if (st) put_simple_story_contents(s, *st, ";");

    st = root->get_dyn_story();
    if (st) put_simple_story_contents(s, *st, ";;");

    for_all_daughters(dyn, root, i)
      s << i->get_index() << ' ' << i->get_mass() << ' ' << i->get_pos() << ' '
	<< i->get_vel() << '\n';
    s << '\n';
  }

}


#else
main(int argc, char** argv)
{
    dyn *b;
    check_help();

    while (b = get_dyn()) {
	cout << "TESTING put_dyn:" << endl;
        put_dyn(b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	delete b;
    }
    cerr << "Normal exit\n";
}
#endif
