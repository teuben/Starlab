
#include "../evolve/sigma.h"
#ifndef TOOLBOX

#define DEFAULT_MASS_RATIO 1
#define DEFAULT_ECC 0
#define DEFAULT_SMA 1
#define SMA_REDUCE 10
#define DEFAULT_R1 0
#define DEFAULT_R2 0
#define TIDAL_TOL_FACTOR 1e-6


local void pp(sdyn* b, ostream & s, int level = 0) {

    s.precision(4);

    for (int i = 0; i < 2*level; i++) s << " ";

    b->pretty_print_node(s);
    s << " \t"<< b->get_mass() << " \t"
      << b->get_pos() << "   "
      << b->get_vel() <<endl;

    for (sdyn * daughter = b->get_oldest_daughter();
	 daughter != NULL;
	 daughter = daughter->get_younger_sister())
	pp(daughter, s, level + 1);	
}

// initialize_root: Set up the outermost (scattering) trajectory.
//                  The plane and orientation of the trajectory define the
//                  coordinate system used.  Assume a target semi-major
//                  axis of 1 in determining the initial separation.

local void initialize_root(sdyn* root, real v_inf, 
			   real rho_sq_min, real rho_sq_max, real peri,
			   real r_init,
			   real projectile_mass, real projectile_radius) {

  //  cerr << "Root: " << endl;
    //PRC(v_inf);PRC(rho_sq_min);PRC(rho_sq_max);PRC(peri);PRC(r_init);
    //    PRC(projectile_mass);PRL(projectile_radius);
    // Explicitly state target mass and radius:

    // NOTE that this will give unit radius to the target star in
    // a two-body encounter!

    real target_mass = 1;
    real target_sma = 1;

    real m_total = target_mass + projectile_mass;
    root->set_mass(m_total);
    root->set_pos(vector(0,0,0));
    root->set_vel(vector(0,0,0));

    sdyn* target = root->get_oldest_daughter();
    target->set_mass(target_mass);
    target->set_radius(target_sma);

    sdyn* projectile = target->get_younger_sister();
    projectile->set_mass(projectile_mass);
    projectile->set_radius(projectile_radius);

    // Establish the orbital elements of the scattering:

    v_inf *= sqrt(0.5*m_total);
    // Exactly one of rho_sq_max and peri should be < 0.

    if (rho_sq_max*peri > 0) err_exit("Inconsistent initial conditions.");

    if(rho_sq_min<0)
      rho_sq_min = 0;
    
    real rho_max = -1;
    if(rho_sq_max>0)
      rho_max = sqrt(rho_sq_max);
    //    PRL(v_inf);
    //    PRL(peri);

    if (peri > 0) {
	if (v_inf > 0)
	    rho_max = peri * sqrt(1 + 2*m_total/(peri*v_inf*v_inf));
	else
	    rho_max = peri;
    }
    // Adjusted by (SPZ: 11 Oct 2000)
    // randomize rho \propto \sqrt{rho_max^2}
    //    PRL(rho_sq_min);PRL(rho_max);
    real rho = sqrt(randinter(rho_sq_min, pow(rho_max, 2)));
    //    PRL(rho);

    // Note: Unit of velocity is v_g, not v_crit (identical for equal masses)

    real energy = .5 * v_inf * v_inf;
    real ang_mom = rho * v_inf;

    //    PRC(energy);PRC(ang_mom);PRC(v_inf);PRL(rho);
    //    PRL(ang_mom);
 
    real ecc = (energy == 0 ? 1	: sqrt( 1 + 2 * energy
				              * pow(ang_mom/m_total, 2)));

    // Special case: if v_inf = 0, assume rho is periastron.

    real virial_ratio = rho;

    kepler k;

    // Don't use mean_anomaly = 0 here (see scatter3.C):

    //PRC(m_total);PRC(energy);PRC(ecc);PRL(virial_ratio);
    make_standard_kepler(k, 0, m_total, energy, ecc, virial_ratio, -0.001, 1);

    // Radius for "unperturbed" inner binary:

    real r_unp = (target_sma + projectile_radius)
	             * pow(TIDAL_TOL_FACTOR / m_total, -1/3.0);

    real r_start = min(r_init, r_unp);
    if (k.get_separation() < r_start)
        k.return_to_radius(r_start);
    else
        k.advance_to_radius(r_start);

    // Initialize the dynamics.

    root->set_time(k.get_time());
    target->set_time(k.get_time());
    projectile->set_time(k.get_time());
    
    target->set_pos(-projectile_mass * k.get_rel_pos() / m_total);
    target->set_vel(-projectile_mass * k.get_rel_vel() / m_total);

    projectile->set_pos(target_mass * k.get_rel_pos() / m_total);
    projectile->set_vel(target_mass * k.get_rel_vel() / m_total);

//    cerr << "Initialized root...\n";
//    pp(root, cerr);

//    cerr << "t: " << target->get_index()
//	 <<" " << target->get_name() << endl;
//    cerr << "p: " << projectile->get_index()
//	 << " "   << projectile->get_name() << endl;
}

// split_particle: Split the specified node into a binary with the specified
//                 parameters.  All unspecified orbit elements are chosen
//                 randomly.  Newly created nodes have names "n1" and "n2",
//                 where "n" is the name of the node being split.

local void split_particle(sdyn* current, real ecc, real sma, int planar,
			  real mass_ratio, real r1, real r2) {

    if (current->get_oldest_daughter() != NULL)
	err_exit("Can't split a binary node!");

    // Update the binary tree structure:

    sdyn* d1 = new sdyn;
    sdyn* d2 = new sdyn;

    current->set_oldest_daughter(d1);

    d1->set_parent(current);
    d2->set_parent(current);

    d1->set_younger_sister(d2);
    d2->set_elder_sister(d1);

    // Set new masses and radii:

    real m_total = current->get_mass();
    real m1 = m_total / (1 + mass_ratio);
    real m2 = m1 * mass_ratio;

    // By convention, first component has larger mass.

    if (m1 < m2) {
	real temp = m1;
	m1 = m2;
	m2 = temp;
    }

    d1->set_mass(m1);
    d2->set_mass(m2);
    
    //    PRC(m1);PRC(r1);PRC(m2);PRL(r2);
    d1->set_radius(r1);
    d2->set_radius(r2);
    
    // Sets parent radius (center of mass) to zero
    d2->get_parent()->set_radius(0);

    // Set internal orbital elements:

    kepler k;

    if (sma < 0) {
	sma = DEFAULT_SMA;
	for (int i = 1; i < strlen(current->get_name()); i++)
	    sma /= SMA_REDUCE;
    }

    if (ecc < 0 || ecc >= 1) {

      // factor 2 for MIN_PERI_FACTOR see sdyn/util/make_tree.C
      real peri_min = 2*max(d1->get_radius(), d2->get_radius());
      real e_max = 1 - peri_min/sma;

      if(e_max<0) 
	err_exit("initial eccentricity out of bounds in split_particle(...)"); 

      // Thermal distribution in eccentricity
      if(STABILITY)
	ecc = sqrt(randinter(0, e_max*e_max));	
      else
	ecc = sqrt(randinter(0,1));	
    }

    //    cerr << "Splitting " << current->get_name();
    //    cerr << ",  planar flag = " << planar << endl;
    //    cerr << "ecc, sma = " << ecc << " " << sma << endl;

    real peri = 1; // Default value (unimportant unless ecc = 1).
    if (ecc == 1) peri = 0;

    // For now, binary phase is random.

    real mean_anomaly = randinter(-PI, PI);

    make_standard_kepler(k, 0, m_total, -0.5 * m_total / sma, ecc,
			 peri, mean_anomaly);

    set_random_orientation(k, planar);

    d1->set_time(current->get_time());
    d2->set_time(current->get_time());

    d1->set_pos(-m2 * k.get_rel_pos() / m_total);
    d1->set_vel(-m2 * k.get_rel_vel() / m_total);

    d2->set_pos(m1 * k.get_rel_pos() / m_total);
    d2->set_vel(m1 * k.get_rel_vel() / m_total);

    // Naming convention:

    d1->set_name(current->get_name());
    d2->set_name(current->get_name());
    strcat(d1->get_name(), "1");
    strcat(d2->get_name(), "2");

//    cerr << "Pretty-print " << current->get_name() << endl;
//    pp(current, cerr);
}

// locate_label_and_set_defaults: Return a pointer to the node with the
//                                specified name, creating the binary tree
//                                above it if necessary and establishing
//                                default parameters at each level.

// Note that the node name is strictly defined and the number of characters in
// the name is the level in the binary tree: "t", "t1", "t11", "t112", etc.

#include <string.h>
local sdyn* locate_label_and_set_defaults(sdyn* root, char* label, int planar) {

    sdyn* b = root;

    for (int i = 0; i < strlen(label); i++) {

	char* name = strdup(label);
	name[i+1] = '\0';		// Truncate name to current level

	sdyn* node = NULL;

//	cerr << "Looking for " << name<<endl;

	while (node == NULL) {
	    for_all_daughters(sdyn, b, bb)
		if (strcmp(bb->get_name(), name) == 0) node = bb;
	    if (node == NULL) {
		if (b->get_oldest_daughter() != NULL) return NULL;

//		cerr << "Not found, making default daughters\n";

		split_particle(b, DEFAULT_ECC, -1, planar,
			       DEFAULT_MASS_RATIO, DEFAULT_R1, DEFAULT_R2);
	    } 
//	    else
//		cerr << "Found " << node->get_name()<<endl;
	}
	b = node;
    }
    return b;
}

// first_leaf: return a pointer to the first leaf descended from the
//             specified node.

// functions removed by SPZ&AG on 15Okt2002
#if 0
sdyn* first_leaf(sdyn* b)
{
    while (b->get_oldest_daughter()) b = b->get_oldest_daughter();
    return b;
}

// next_leaf: return a pointer to the next leaf after the present node.

sdyn* next_leaf(sdyn* b)
{
    if (b == NULL)

	return NULL;

    else if (b->get_oldest_daughter())

	return first_leaf(b);

    else // Ascend the tree looking for a younger sister.

	while (b != NULL && b->get_younger_sister() == NULL)
	    b = b->get_parent();

    return (b == NULL ? (sdyn*)NULL : first_leaf(b->get_younger_sister()));
}

#endif

sdyn* mkscat(int argc, char **argv, sigma_input &input) {

    // Establish standard configuration and names:

    sdyn* root = mksdyn(2);	// Top-level (unbound) scattering orbit.
    root->set_name("r");
    root->set_index(-1);	// Undo indexing effect of mksdyn...

    sdyn* target = root->get_oldest_daughter();
    target->set_name("t");
    target->set_index(-1);

    sdyn* projectile = target->get_younger_sister();
    projectile->set_name("p");
    projectile->set_index(-1);

    sdyn* current = root;	// Node currently under consideration.

    // Establish defaults for the top level:

    real projectile_mass = 1;	// Note convention: projectile and target have
    real projectile_radius = 0;	// 		    EQUAL masses by default.
    real v_inf = 1;
    real rho = 0;
    real peri = -1;
    real initial_separation = VERY_LARGE_NUMBER;

    real mass_ratio = DEFAULT_MASS_RATIO;   // Parameters to use when the
    real ecc = DEFAULT_ECC;		    // next binary is built
    real sma = DEFAULT_SMA;
    real r1 = DEFAULT_R1;
    real r2 = DEFAULT_R2;

    int  planar = 0;
    bool debug = FALSE;

    // First scan the argument list to get the actual random seed used:

    int  seed = 0;
    int i = 0;
    int random_seed;
    while (++i < argc)
	if (argv[i][0] == '-' && argv[i][1] == 's') {
	  seed = atoi(argv[++i]);
	  if(seed!=0) 
	    random_seed = srandinter(seed);
	  else
	    cerr << "Seed not re-initialized" << endl;
	  seed = -1;
	}
    if(seed>0) 
      random_seed = srandinter(seed);

    // Now parse the rest of the command line and initialize the system.

    i = 0;
    while (++i < argc) if (argv[i][0] == '-')
	switch (argv[i][1]) {

	    // Top-level parameters:

	    case 'M': if (current != root) cerr <<
		        "Too late to initialize projectile mass!\n";
		      projectile_mass = atof(argv[++i]);
		      input.pmass = projectile_mass;
		      break;
	    case 'R': projectile_radius = atof(argv[++i]);
		      break;
	    case 'r': switch(argv[i][2]) {
			  case '\0':    if (current != root) cerr <<
			      "Too late to initialize impact parameter!\n";
			                if(input.rho_sq_max<0 && input.peri<0) {
					  rho = atof(argv[++i]);
					  peri = -1;
					  input.rho_sq_max = rho*rho;
					  input.peri = -1;
			                }
					break;
			  case 'm':    if (current != root) cerr <<
			      "Too late to initialize pericenter!\n";
			                if(input.rho_sq_max<0 && input.peri<0) {
					  peri = atof(argv[++i]);
					  rho = -1;
					  input.peri = peri;
					  input.rho_sq_max = rho;
					}
					break;
			  case '1':	r1 = atof(argv[++i]);
					break;
			  case '2':	r2 = atof(argv[++i]);
					break;
			  default:      cerr << "Incorrect 'r' flag ignored\n";
					break;
		      }
		      break;
	    case 'S': initial_separation = atof(argv[++i]);
		      break;
	    case 'v': if (current != root) cerr <<
		        "Too late to initialize velocity at infinity!\n";
	              if(input.rho_sq_max<0 && input.peri<0) {
			v_inf = atof(argv[++i]);
			input.v_inf = v_inf;
		      }
		      break;

	    // Begin accumulating parameters on a new node:

	    case 'p':
	    case 't': if (current == root)
			  initialize_root(root, input.v_inf, 
					  input.rho_sq_min, input.rho_sq_max, 
					  input.peri,
					  initial_separation,
					  projectile_mass, projectile_radius);
		      else
			  split_particle(current, ecc, sma, planar,
					 mass_ratio, r1, r2);

		      current = locate_label_and_set_defaults(root,
							      &argv[i][1],
							      planar);
		      if (current == NULL) err_exit("Illegal node name.");

		      // Reset the defaults (note that planar isn't reset):

		      mass_ratio = DEFAULT_MASS_RATIO;
		      //		      ecc = DEFAULT_ECC;
		      //ecc = sqrt(randinter(0, 1));
		      ecc = -1;
		      sma = -1;	       	// Impossible value
		      r1 = DEFAULT_R1;
		      r2 = DEFAULT_R2;

		      break;

	    // Binary parameters:

	    case 'a': sma = atof(argv[++i]);
		      break;
	    case 'e': ecc = atof(argv[++i]);
		      break;
	    case 'q': mass_ratio = atof(argv[++i]);
		      break;

	    case 'P': switch(argv[i][2]) {
		          case '-':	planar = -1;
					break;
		          case '0':	planar = 0;
					break;
		          case '+':	planar = 1;
					break;
		          case '\0':	if (planar == 0)
			                    planar = 1;
					else
					  planar = 0;
					break;
			  default:      cerr << "Incorrect 'P' flag ignored\n";
					break;
		      }
			break;

	    case 'd': debug = 1 - debug;
		      break;
	    case 's': break;

            default:  cerr << "usage: mkscat [-d] [-s #]"
		           << " [-M #] [-r #] [-rm #] [-v #] [-R #] [-S #]"
		           << " [ -t/p... [-a #] [-e #] [-q #] [-P[+/-]"
			   << " [-r1 #] [-r2 #] ] [-t/p... ...]"
			   << endl;
		      exit(1);
	}

    // Initialize the last piece:

    if (current == root)
	initialize_root(root, input.v_inf, input.rho_sq_min, input.rho_sq_max, 
			input.peri, 
			initial_separation,
			projectile_mass, projectile_radius);
    else
	split_particle(current, ecc, sma, planar,
		       mass_ratio, r1, r2);

    // Finally, assign sequential indices to the particles, for the
    // convenience of xstarplot and other display/reduction programs.

    int id = 0;
    sdyn* b = first_leaf(root);
    while (b) {
	b->set_index(++id);
	b = next_leaf(b);
    }

//----------------------------------------------------------------------------

    root->log_history(argc, argv);

    char seedlog[80];
    sprintf(seedlog, "           random seed = %d",
	    random_seed);
    root->log_comment(seedlog);

    if (0) {
	cerr << "mkscat: pretty-print system (seed = "
	     << random_seed << "):\n";
	pp(root, cerr);
    }

    return root;
}

#define MAX_ARGS 256
#define TEXT_BUFFER_SIZE 1024

static char temp[TEXT_BUFFER_SIZE];	// argv will point into this array
static char* name = "parse_string";	// Any name will do...

// parse_string: Split a character string into words.  Return a list of
//               substrings, using the same conventions as argc and argv
//               in UNIX (0 --> "parse_string").

#if 0
void parse_string(char* s, int& argc, char* argv[]) {

    int length = strlen(s);
    if (length >= TEXT_BUFFER_SIZE) err_exit("String too long.");

    argv[0] = name;
    argc = 1;

    int first = 0, last = -1;  // First and last characters of the current word

    for (int i = 0; i <= length; i++) {
	temp[i] = s[i];
	if (temp[i] == ' ' || i == length) {
	    if (last >= 0) {				// End of word
		temp[i] = '\0';
		argv[argc++] = temp + first;		// *No* check on argc
	    }
	    last = -1;
	} else {
	    if (last < 0) first = i;			// Start of word
	    last = i;
	}
    }
}
#endif

sdyn* mkscat(char* s, sigma_input &input) {	
                                    // Character string interface to mkscat
    int argc;
    char* argv[MAX_ARGS];

    parse_string(s, argc, argv);
    if (argc > MAX_ARGS) err_exit("Too many arguments.");

    return mkscat(argc, argv, input);
}    

#endif
