 
//// conv_star_to_dyn: converts stellar system to N-body system
////         N-body units are M=G=-4E = 1
////
////         NOTE: conversion units are taken from input snapshot.
////
//// Options:    -e #  softening
////             -v    express units in terms of the virial radius and
////                   virial radius crossing time.
////
////                   convertsion from 
////                   general standard input units 
////                   [Msun, pc, km/s, Myear]      
////                   to dimenson less N-body units
////
//  SPZ @MIT, August 2001

#include "dyn.h"
#include "constants.h"

#ifdef TOOLBOX

#define  FALSE  0
#define  TRUE   1

main(int argc, char ** argv) {

    bool virial_units = false;
    real eps = 0;

    dyn *root;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "e:v";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {
	    case 'e': eps = atof(poptarg);
	              break;
	    case 'v': virial_units = true;
	              break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	}

    // real input snapshot
    root = get_dyn();
    root->log_history(argc, argv);
    
    real mtot = root->get_mass();
    PRL(mtot);

    // Calculate unit for G (called Newton_G)
    // Note here that the assumed units are:
    // mass in [Msun]
    // position in [pc]
    // velocity in [km/s]
    // time in [Myr]
    // Yes, these are strange units.
    real Uk = cnsts.parameters(Msun)*square(cnsts.physics(km_per_s));
    real Up = cnsts.physics(G)*square(cnsts.parameters(Msun))
            / cnsts.parameters(PC);
    real Newton_G = Up/Uk;

    real P, K, E;
    get_top_level_energies(root, eps*eps, P, K);
    P *= Newton_G;
    E = K + P;
    real Q = -K/P;

    cerr << "basis scaling units" << endl;
    cerr << "Unit for kinetic energy:   " << K << endl;
    cerr << "Unit for potential energy: " << P << endl;
    cerr << "Unit for total energy:     " << E << endl;
    cerr << "Unit for virial ratio:     " << Q << endl;
    cerr << "Unit for G:              1/" << 1/Newton_G << endl;

    cerr << "Virial units" << endl;
    real r_vir = -0.5 * Newton_G * pow(mtot, 2) / P;
    PRL(r_vir);
    real t_vir= 21.0 * sqrt(Q*pow(r_vir, 3) / mtot);
    PRL(t_vir);

    real Um = mtot;
    real Ul = Newton_G * pow(Um, 2)/(-4*E);
    real Ut = Newton_G * pow(Um, 5./2.) / pow(-4*E, 3./2.);
        
    if(virial_units || Ul<0) {
      cerr << "Use virial radius instead of distance unit" << endl;
      Ul = r_vir;
      cerr << "Use virial radius crossing time as time unit" << endl;
      Ut = t_vir;
    }

    real Uv = Ul/Ut;

    cerr << "Unit for M:                " << Um << endl;
    cerr << "Unit for R:                " << Ul << endl;
    cerr << "Unit for T:                " << Ut << endl;
    cerr << "Unit for v:                " << Uv << endl;

    // Now transform the snapshot
    for_all_leaves(dyn, root, bi) {
      bi->set_mass(bi->get_mass() / Um);
      bi->set_pos(bi->get_pos() / Ul);
      bi->set_vel(bi->get_vel() / Uv);
    }

    // Now initialize the units for stellar evolution
    // These units are curious too:
    // mass unit is 1/mtot
    real mf = Um;
    // distance unit is Rsun/pc
    real rf = Ul;
      // cnsts.parameters(Rsun)/(cnsts.parameters(PC)*Ul);
    // time unit is 1/Myear
    real tf = Ut; //cnsts.physics(Myear);
    PRC(mf);PRC(rf);PRL(tf);

    root->get_starbase()->set_stellar_evolution_scaling(mf, rf, tf);
    root->get_starbase()->print_stellar_evolution_scaling(cerr);

    root->set_mass(1);
    put_node(root);
}

#endif

// end of: conv_star_to_dyn.C
 

