
#include "scatter3.h"
#include "pvm_scatt.h"
#include "pvm3.h"

main(int argc, char **argv)
{
    initial_state3 init;
    make_standard_init(init);

    int  seed 	    = 0;    	// seed for random number generator
    int n_rand      = 0;        // number of times to invoke the generator
                                // before starting for real
    int  n_experiments = 1;     // default: only one run
    real dt_out     =       	// output time interval
	  VERY_LARGE_NUMBER;
    real dt_snap    =       	// output time interval
	  VERY_LARGE_NUMBER;

    real cpu_time_check = 3600;
    real snap_cube_size = 10;

    int planar_flag = 0;
    int psi_flag = 0;
    real psi = 0;

    int   b_flag = 0;
    bool  q_flag = FALSE;
    bool  Q_flag = FALSE;

    extern char *poptarg;
    int c;
    const char *param_string = "A:bc:C:d:D:e:g:L:m:M:n:N:o:pPqQr:R:s:S:U:v:x:y:z:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1)
	switch(c) {

	    case 'A': init.eta = atof(poptarg);
		      break;
	    case 'b': b_flag = 1 - b_flag;
		      break;
	    case 'c': cpu_time_check = 3600*atof(poptarg);// (Specify in hours)
		      break;
	    case 'C': snap_cube_size = atof(poptarg);
		      break;
	    case 'd': dt_out = atof(poptarg);
		      break;
	    case 'D': dt_snap = atof(poptarg);
		      break;
	    case 'e': init.ecc = atof(poptarg);
		      break;
	    case 'g': init.tidal_tol_factor = atof(poptarg);
		      break;
	    case 'L': init.r_init_min = atof(poptarg);
		      break;
	    case 'm': init.m2 = atof(poptarg);
		      break;
	    case 'M': init.m3 = atof(poptarg);
		      break;
	    case 'n': n_experiments = atoi(poptarg);
		      break;
	    case 'N': n_rand = atoi(poptarg);
		      break;
	    case 'o': psi = atof(poptarg);
		      psi_flag = 1;
		      break;
	    case 'p': planar_flag = 1;
		      break;
	    case 'P': planar_flag = -1;
		      break;
	    case 'q': q_flag = 1 - q_flag;
		      break;
	    case 'Q': Q_flag = 1 - Q_flag;
		      break;
	    case 'r': init.rho = atof(poptarg);
		      break;
	    case 'R': init.r_stop = atof(poptarg);
		      init.r_init_min = init.r_init_max = abs(init.r_stop);
		      break;
	    case 's': seed = atoi(poptarg);
		      break;
	    case 'S': init.r_stop = atof(poptarg);
		      break;
	    case 'U': init.r_init_max = atof(poptarg);
		      break;
	    case 'v': init.v_inf = atof(poptarg);
		      break;
	    case 'x': init.r1 = atof(poptarg);
		      break;
	    case 'y': init.r2 = atof(poptarg);
		      break;
	    case 'z': init.r3 = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}            

    if (Q_flag) q_flag = TRUE;

    if (init.m2 > 1) {
	cerr << "scatter3: init.m2 = " << init.m2 << " > 1" << endl;
	exit(1);
    }

    cpu_init();
    int random_seed = srandinter(seed, n_rand);

    for (int i = 0; i < n_experiments; i++) {

	if (n_experiments > 1) cerr << i+1 << ": ";

	cerr << "Random seed = " << get_initial_seed()
	     << "  n_rand = " << get_n_rand() << endl << flush;

	// Normally, we just want random angles and phases.
	// However, in the planar case, it may be desireable to
	// control the orientation of the outer orbit.  The
	// angle between the inner and outer orbital axes in
	// the planar case is init.phase.psi.  Specify psi
	// in *degrees*, with psi = 0 meaning that the outer
	// periastron lies along the positive x-axis.

	randomize_angles(init.phase);

	if (planar_flag == 1) {
	    init.phase.cos_theta = 1;	// Planar prograde
	    if (psi_flag) init.phase.psi = psi * M_PI / 180.0;
	} else if (planar_flag == -1) {
	    init.phase.cos_theta = -1;	// Planar retrograde
	    if (psi_flag) init.phase.psi = psi * M_PI / 180.0;
	}

	intermediate_state3 inter;
	final_state3 final;

	real cpu = cpu_time();

	// scatter3(init, inter, final, cpu_time_check,
	//	    dt_out, dt_snap, snap_cube_size);

	//-----------------------------------------------------------------

	// Perform scattering on a separate processor.

	char exe[100];

	strcpy(exe, getenv("STARLAB_PATH"));	// ASSUME that these
	if (exe[0] == '\0')			// environment variables
	    strcpy(exe, getenv("PWD"));		// are defined
	else
	    strcat(exe, "/bin");
	strcat(exe, "/scatter3_slave.pvm");

	int tid;
	if (pvm_spawn(exe, (char**)0, 33, ".", 1, &tid) != 1) {
	    cerr << "Unable to spawn slave task\n";
	    exit(1);
	}

	if (pvm_recv(tid, HANDSHAKE_MSG) < 0) {
	    cerr << "Error getting slave host ID\n";
	    exit(1);
	}

	char name[100];
	pvm_upkstr(name);

	cerr << "initialized slave node " << name << endl;
	cerr << "tid = " << tid << endl;

	pvm_notify(PvmTaskExit, FAILURE_MSG, 1, &tid);

	// Pack the data.

	pvm_initsend(PvmDataRaw);

	int n_rand = get_initial_seed();
	pvm_pkint(&n_rand, 1, 1);
	pack(init);

	pvm_pkdouble(&cpu_time_check, 1, 1);
	pvm_pkdouble(&dt_snap, 1, 1);
	pvm_pkdouble(&snap_cube_size, 1, 1);

	// Send the data as message SEND_DATA_MSG.

	pvm_send(tid, SEND_DATA_MSG);

	struct timeval timeout;
	timeout.tv_sec = 3600;		// Give up after an hour.
	int bufid, length, msgid, tid_r;

	bufid = pvm_trecv(-1, -1, &timeout);
	pvm_bufinfo(bufid, &length, &msgid, &tid_r);

	cerr << "returned msgid = " << msgid << endl;
	if (msgid != RETURN_DATA_MSG) exit(0);

	pvm_upkint(&n_rand, 1, 1);

	// Note: n_rand should = random_seed[task]

	unpack(init);
	unpack(inter);
	unpack(final);

	int dum;
	pvm_upkint(&dum, 1, 1);
	pvm_upkdouble(&cpu, 1, 1);

	//-----------------------------------------------------------------

	// cpu = cpu_time() - cpu;

	cerr << ":  ";
	print_scatter3_outcome(inter, final, cerr);

	if (Q_flag) print_scatter3_summary(inter, final, cpu, cerr);

	if (!q_flag) print_scatter3_report(init, inter, final,
					   cpu, b_flag, cerr);
    }
}
