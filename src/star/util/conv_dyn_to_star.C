 
//// conv_star_to_dyn: converts N-body system in
////         N-body units (M=G=-4E=1) to 
////         general standard units [Msun, pc, km/s, Myear]      
////
////         NOTE: conversion units are taken from input snapshot.
////
//// Options:    none
////
////
//  SPZ @MIT, August 2001

#include "dyn.h"
#include "constants.h"

#ifdef TOOLBOX

//#else

#define  FALSE  0
#define  TRUE   1

main(int argc, char ** argv) {

    int option  = 1;  // conv to cgs

    dyn *root;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    root = get_dyn(cin);
    root->log_history(argc, argv);

    real P, K, E;
    get_top_level_energies(root, 0.0, P, K);
    real r_virial = -0.5*root->get_mass()*root->get_mass()
                  / P;
    PRL(r_virial);

    real mf = root->get_starbase()->conv_m_star_to_dyn(1);
    real rf = root->get_starbase()->conv_r_star_to_dyn(1);
    real tf = root->get_starbase()->conv_t_star_to_dyn(1);
    PRC(mf);PRC(rf);PRL(tf);
    real mtot = 1/mf;
    real rvir = cnsts.parameters(Rsun) / (rf * cnsts.parameters(PC));
    real tvir = 1/tf;
    real vel = rvir/tvir; // [pc/Myr]
    real km_sec = cnsts.physics(km_per_s) * cnsts.physics(Myear) 
                   / cnsts.parameters(PC);
    vel /= km_sec;
    PRL(km_sec);
    PRC(mtot);PRC(rvir);PRC(tvir);PRL(vel);
    real Um = 1/mtot;
    real Ul = r_virial/rvir;  // virial radius may be non-unity
    real Ut = 1/tvir;
    real Uv = 1/vel;

    // Now transform the snapshot
    for_all_leaves(dyn, root, bi) {
      bi->set_mass(bi->get_mass() / Um);
      bi->set_pos(bi->get_pos() / Ul);
      bi->set_vel(bi->get_vel() / Uv);
    }
    root->set_mass(mtot);

    cerr << "Reset the standard starbase units" << endl;
    root->get_starbase()->set_stellar_evolution_scaling(-1, -1, -1);
    root->get_starbase()->print_stellar_evolution_scaling(cerr);

    put_node(cout, *root);
}

#endif

// end of: conv_dyn_to_star.C
 


