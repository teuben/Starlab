
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\     
 //                                                       //            _\|/_
//=======================================================//              /|\

#include "node.h"

//see: Eggleton, Fitchet & Tout 1989, ApJ 347, 998

local real mf_Miller_Scalo(real m_lower, real m_upper)
{
    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.19*rnd
	    / (pow(1-rnd, 0.75) + 0.032*pow(1-rnd, 0.25));
    }
    while(m_lower>m || m>m_upper); 
    return m;
}

// see: de la Fuente Marcos, Aarseth, Kiseleva, Eggleton 
// 1997, in Docobo, Elipe, McAlister (eds.), Visual Double
// Stars: Formation, dynamics and Evolutionary Tracks, KAP: ASSL
// Series vol. 223,  165

local real mf_Scalo(real m_lower, real m_upper) {
    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.3*rnd
	    / pow(1-rnd, 0.55);
    }
    while(m_lower>m || m>m_upper); 
    return m;
}

real get_random_stellar_mass(real m_lower, real m_upper) 
{
    //return mf_Miller_Scalo(m_lower, m_upper);
    return mf_Scalo(m_lower, m_upper);
}

void main(int argc, char ** argv)
{
    int n = 1000;
    int count = 0, countNS = 0, countBH = 0;
    int visible = 0, OBvisible = 0;
    real mlimit = 20;

    if (argc > 1) n = atoi(argv[1]);
    if (argc > 2) mlimit = atof(argv[2]);
    srandinter(0);

    for (int i = 0; i < n; i++) {

	real m = get_random_stellar_mass(0.1, 100);

	real t = randinter(0, 12);	// Gyr in past
	real lifetime = 12;
	if (m > 1) {
	    lifetime *= pow(m, -2.5);
	    if (lifetime < 0.01) lifetime = 0.01;
	}
	
	count++;
	if (m > mlimit)
	    countBH++;
	else if (m > 8)
	    countNS++;

	if (lifetime >= t) {
	    visible++;
	    if (m > 5) OBvisible++;
	}

//	cerr << m << endl;

    }
    cerr << "BH prog. frac = " << countBH / ((real) count)
	 << "  NS prog. frac = " << countNS / ((real) count) << endl;
    PRC(n), PRC(visible), PRL(OBvisible);
}
