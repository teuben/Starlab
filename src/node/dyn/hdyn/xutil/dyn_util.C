
#include "hdyn.h"
#include "xstarplot.h"

void convert_relative_to_absolute(hdyn* b)
{
    if (b->get_parent()) b->inc_pos(b->get_parent()->get_pos());
    for_all_daughters(hdyn, b, bb) convert_relative_to_absolute(bb);
}

void accumulate_potential_energy(hdyn* bj, hdyn*bi,
				 real& epot, real& rmin, hdynptr& bmin)

// Determine the potential energy of bi with respect to bj, along
// with bi's nearest neighbor and minimum neighbor distance.

{
    if (bj->get_oldest_daughter() != (hdyn*)NULL)
	for (hdyn * bb = bj->get_oldest_daughter(); bb != (hdyn*)NULL;
	    bb = bb->get_younger_sister()) {
	    accumulate_potential_energy(bb, bi, epot, rmin, bmin);
	}
    else
	if (bi != bj) {
	    vec d_pos = bi->get_pos() - bj->get_pos();
	    real mi = bi->get_mass();
	    real mj = bj->get_mass();
	    real r = sqrt(d_pos * d_pos);
	    epot += -mi*mj/r;
	    if (r < rmin) {rmin = r; bmin = bj;}
	}
}

// compute_energies: calculate the appropriate color code for particle bi
//		     relative to "particle" bj (usually the root node).

void compute_energies(hdyn* bj, hdyn* bi, char& c)
{
    c = default_color;

    real   epot = 0, rmin = VERY_LARGE_NUMBER;
    hdynptr bmin = bi;

    accumulate_potential_energy(bj, bi, epot, rmin, bmin);

    real   mi = bi->get_mass();
    real   mj = bmin->get_mass();
    vec d_vel = bi->get_vel() - bmin->get_vel();
    vec d_pos = bi->get_pos() - bmin->get_pos();
    real   r = sqrt(d_pos * d_pos);
    real   e0 = (0.5 * d_vel * d_vel - (mi + mj)/r);

    if (e0 < 0) {
	real   epot1 = 0, rmin1 = VERY_LARGE_NUMBER;
	hdynptr bmin1 = bi;

	accumulate_potential_energy(bj, bmin, epot1, rmin1, bmin1);

	if (bi == bmin1) {
	    real e  = - mi*mj / r;
	    vec R_vel = (mi*bi->get_vel()+mj*bmin->get_vel())/(mi+mj);
	    real ekin = 0.5*(mi+mj)*R_vel*R_vel;

	    if (epot + epot1 - 2*e + ekin < 0) c = bound_binary;
	    else c = unbound_binary;

	} else c = bound_single;

    } else {
	vec vel = bi->get_vel();
	real ekin = 0.5*bi->get_mass()*vel*vel;
	
	if (ekin + epot > 0.0) c = unbound_single;
    }
}

int nearest_index(hdyn* b, float r, float s, int kx, int ky)
{
    real r2_min = VERY_LARGE_NUMBER;
    int index = -1;
    for_all_leaves(hdyn, b, bb) {
	real dx = bb->get_pos()[kx] - r;
	real dy = bb->get_pos()[ky] - s;
	real r2 = dx*dx + dy*dy;
	if (r2 < r2_min) {
	    r2_min = r2;
	    index = bb->get_index();
	}
    }
    return index;
}

int count_nodes(hdyn* b)
// Return the total number of nodes BELOW the node b.
{
    if(b->get_oldest_daughter() == NULL)
	return 0;
    else {
	int  n = 0;
	for (hdyn* daughter = b->get_oldest_daughter();
	     daughter != NULL; daughter = daughter->get_younger_sister())
	    n += 1 + count_nodes(daughter);
	return n;
    }
}
