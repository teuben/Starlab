
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

#include "scatter3.h"

// set_orientation:  Set up the orientation of the outer orbit with
//		     respect to the inner binary [which always lies in
//		     the (x-y) plane].   See "sigma3.h" for details.

void set_orientation(kepler &k, phase3 &p)
{
    real mu = p.cos_theta;
    real sin_theta = sqrt(1 - mu * mu);

    // Construct the normal vector:

    vec n = vec(sin_theta*cos(p.phi), sin_theta*sin(p.phi), mu);

    // Construct unit vectors a and b perpendicular to n:

    vec temp = vec(1, 0, 0);
    if (abs(n[0]) > .5) temp = vec(0, 1, 0);	// temp is not parallel to n
    if (n[2] < 0) temp = -temp;

    vec b = n ^ temp;
    b /= abs(b);
    vec a = b ^ n;
    if (n[2] < 0) a = -a;	// Force (a, b) to be (x, y) for n = +/-z
    
    // Construct *random* unit vectors l and t perpendicular to each
    // other and to n (psi = 0 ==> periastron along a):

    vec l = cos(p.psi)*a + sin(p.psi)*b;
    vec t = n ^ l;

    k.set_orientation(l, t, n);
    k.initialize_from_shape_and_phase();
}

// set_up_dynamics:  Turn the system described by two kepler structures
//		     into an sdyn3 tree; return a pointer to the root
//		     node of the tree.

sdyn3 * set_up_dynamics(real m2,      // mass of secondary in binary
			real m3,      // mass of projectile (m1 + m2 = 1)
			kepler & k1,  // inner binary
			kepler & k3)  // outer binary
{
    // Create bodies for integration.

    sdyn3 * b = mksdyn3(3);            	    // pointer to the 3-body system

    // Initialize the dynamics.

    b->set_time(k3.get_time());
    for_all_daughters(sdyn3, b, bb)
	bb->set_time(k3.get_time());

    sdyn3 * b1 = b->get_oldest_daughter();
    sdyn3 * b2 = b1->get_younger_sister();
    sdyn3 * b3 = b2->get_younger_sister();

    b1->set_label(1);
    b2->set_label(2);
    b3->set_label(3);

    b->set_mass(1 + m3);
    b1->set_mass(1 - m2);
    b2->set_mass(m2);
    b3->set_mass(m3);

    kepler_pair_to_triple(k1, k3, b1, b2, b3);

    return b;
}

// sdyn3_to_system:  Save sdyn3 dynamical information in the specified array.

void sdyn3_to_system(sdyn3 * root, body * system)
{
    int n = 0;
    for_all_daughters(sdyn3, root, bb) {
	system[n].index = bb->get_index();
	system[n].mass  = bb->get_mass();
	for (int kp = 0; kp < 3; kp++) system[n].pos[kp] = bb->get_pos()[kp];
	for (int kv = 0; kv < 3; kv++) system[n].vel[kv] = bb->get_vel()[kv];
	n++;
    }
}
