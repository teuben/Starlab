//
// proto_star.C
//

#include "proto_star.h"
#include "dyn.h"
#include "kepler.h"

local void add_secondary(dyn* original, real mass_ratio)
{
    dyn* primary = new dyn;
    dyn* secondary = new dyn;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    primary->set_mass(original->get_mass());
    secondary->set_mass(mass_ratio*original->get_mass());
    original->inc_mass(secondary->get_mass());

    // Naming convention:

    if (original->get_name() == NULL)
	if (original->get_index() >= 0) {
	    char tmp[64];
	    sprintf(tmp, "%d", original->get_index());
	    original->set_name(tmp);
	}

    if (original->get_name() != NULL) {

	// Label components "a" and "b".

	primary->set_name(original->get_name());
	secondary->set_name(original->get_name());
	strcat(primary->get_name(), "a");
	strcat(secondary->get_name(), "b");

	// Make standard CM name.

	char tmp[256];
	sprintf(tmp, "(%s,%s)", primary->get_name(), secondary->get_name());
	original->set_name(tmp);
    }
}

#if 0
// To be worked out
local void set_kepler_orientation_from_aml(kepler *k, 
					   vec angular_momentum) {

    // normal unit vector is normlaized angular momentum vector
    vec n = angular_momentum/abs(angular_momentum);

    vec temp = vec(1, 0, 0);
    if (abs(n[0]) > .5) temp = vec(0, 1, 0);	// temp is not parallel to n
    if (n[2] < 0) temp = -temp;

    vec b = n ^ temp;
    b /= abs(b);
    vec a = b ^ n;
    if (n[2] < 0) a = -a;	// Force (a, b) to be (x, y) for n = +/-z
    
    // Construct *random* unit vectors l and t perpendicular to each
    // other and to n (psi = 0 ==> periastron along a):

    vec l = cos(psi)*a + sin(psi)*b;
    vec t = n ^ l;

    k.set_orientation(l, t, n);
}
#endif

local void add_dynamics(dyn* cm, real ecc, real energy) 
		//	vec angular_momentum)
{
    dyn* primary = cm->get_oldest_daughter();
    dyn* secondary = primary->get_younger_sister();

    real m_total = cm->get_mass();
    real m1 = primary->get_mass();
    real m2 = secondary->get_mass();

    PRC(m_total);PRC(m1);PRL(m2);

    // Set internal orbital elements:

    kepler k;

    real peri = 1; // Default value (unimportant unless ecc = 1).
    if (ecc == 1) peri = 0;

    // For now, binary phase is random.

    real mean_anomaly = randinter(-PI, PI);

    // Energies here are really binding energies (i.e. > 0 for a bound
    // orbit) per unit mass; kepler package expects energy < 0.

    energy = -energy;

    make_standard_kepler(k, 0, m_total, energy, ecc,
			 peri, mean_anomaly);

    //PRC(m_total);
    //PRC(energy);
    //PRL(ecc);

    //Orientation for protostar should not be random.
    //set_random_orientation(k);
    // ceate normal vector from angular momentum vector.
    
    //    set_kepler_orientation_from_aml(k, angular_momentum);
    int planar = 0; // ??
    set_random_orientation(k, planar);

    k.initialize_from_shape_and_phase();

    k.print_elements(cerr);

    // Set positions and velocities.

    primary->set_pos(-m2 * k.get_rel_pos() / m_total);
    primary->set_vel(-m2 * k.get_rel_vel() / m_total);

    cerr << primary->get_pos() << endl;
    cerr << primary->get_vel() << endl;

    secondary->set_pos(m1 * k.get_rel_pos() / m_total);
    secondary->set_vel(m1 * k.get_rel_vel() / m_total);

    cerr << secondary->get_pos() << endl;
    cerr << secondary->get_vel() << endl;
}

void proto_star::create_binary_from_proto_star() {

  //    star_transformation_story(Main_Sequence);
  //    new main_sequence(*this);
  //    return;

    real m_tot = core_mass;
    real q = randinter(0, 1);
    real mprim = m_tot/(1+q);
    real msec = mprim*q;

    //	angular_momentum = sqrt(2 * total_mass * periastron);

    real am = 1;
    //    real am = angular_momentum();

    real a_min   = pow(am, 2)/m_tot;
    real two_r   = 2*core_radius;
    real ecc_min = Starlab::max(0., 1 - pow(am/(two_r*m_tot), 2));
    real ecc_max = Starlab::min(1., 1 - two_r/a_min);

    //	angular_momentum = sqrt(abs(1 - eccentricity * eccentricity)
    //	                         * total_mass * semi_major_axis);

    //	periastron = semi_major_axis * abs(1 - eccentricity);
    PRC(am);PRC(q);PRC(mprim);PRC(msec);PRC(a_min);PRL(two_r);    
    PRC(ecc_min);PRL(ecc_max);
    
//    if ((!cnsts.parameters(proto_star_to_binary) ||
//	  get_use_dyn()) && (ecc_max<0 || ecc_min>=1)) {

       star_transformation_story(Main_Sequence);
       new main_sequence(*this);
       return;
//    }

    star_transformation_story(Double);

    //Was squared
    ecc_min = sqrt(ecc_min);

    real ecc = sqrt(randinter(ecc_min, ecc_max));
    real sma = pow(am, 2) / (m_tot*(1-pow(ecc, 2)));
    real energy = 0.5 * (mprim+msec) / sma;

    dump(cerr, false);
    PRC(ecc);PRC(sma);PRL(energy);

    relative_mass = mprim;
    get_node()->set_mass(conv_m_star_to_dyn(mprim));

    add_secondary(dynamic_cast(dyn*, get_node()), q);
    cerr << "Secondary added"<<endl;

    addstar(((dyn*)get_node())->get_parent(), 0., Main_Sequence);
    cerr << "star added"<<endl;

    add_dynamics(dynamic_cast(dyn*, get_node()->get_parent()), 
		 ecc, energy);

    cerr << "dynamics double created" << endl;

    for_all_leaves(dyn, ((dyn*)get_node())->get_root(), si) {
      cerr << si->format_label() << " ";
      PRL(si->get_starbase()->get_element_type());
      cerr << "pos: " << si->get_pos() << endl;
    }
}






