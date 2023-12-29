
#include "dyn.h"

void shift_to_com(dyn *b)
{
    // Compute the center of mass.

    real total_mass = 0;
    vec com_pos  = 0;
    vec com_vel  = 0;

    for_all_daughters(dyn, b, d) {

	total_mass += d->get_mass();
	com_pos    += d->get_mass() * d->get_pos();
	com_vel    += d->get_mass() * d->get_vel();

    }	

    com_pos /= total_mass;
    com_vel /= total_mass;

    // Shift all positions and velocities.

    for_all_daughters(dyn, b, d) {

        d->inc_pos(-com_pos);
        d->inc_vel(-com_vel);

    }
}

main()
{
    dyn *b = get_dyn(cin);

    shift_to_com(b);

    put_dyn(b, cout);
}
