
// scatter3_pvm.C: Perform three-body scattering experiments.
//
// Input is an initial state, output is intermediate and final states.
// All data are received and returned via PVM message-passing calls.
// 
// Note the function scatter3 is completely deterministic.
// No randomization is performed at this level.

#ifdef HAS_PVM

#include "scatter3.h"
#include "pvm_scatt.h"
#include "pvm3.h"
#include <unistd.h>

extern "C" int gethostname(char*, unsigned int);

#define FORCE_TIMEOUT	0
#define TIME_TIMEOUT	100	// (seconds)

void main()
{
    // Start by echoing our hostname back to the parent process.

    char name[100];
    gethostname(name, 100);

    pvm_initsend(PvmDataRaw);
    pvm_pkstr(name);
    pvm_send(pvm_parent(), HANDSHAKE_MSG);

    while (1) {

	// Pick up the next message from parent and act accordingly.

	int bufid = pvm_recv(pvm_parent(), -1);

	int length, msgid, source;
	int info = pvm_bufinfo(bufid, &length, &msgid, &source);

	if (info < 0 || msgid != SEND_DATA_MSG || source != pvm_parent()) {
	    pvm_exit();
	    exit(0);
	}

	// Unpack the data.

	int n_rand;
	pvm_upkint(&n_rand, 1, 1);

	initial_state3 init;
	unpack(init);

	real cpu_time_check, dt_out, dt_snap, snap_cube_size;
	pvm_upkdouble(&cpu_time_check, 1, 1);
	pvm_upkdouble(&dt_snap, 1, 1);
	pvm_upkdouble(&snap_cube_size, 1, 1);

	// cerr << "unpacked data: m2, m3 = " << init.m2 <<" "<< init.m3
	//	<< endl;
	// cerr << "starting scattering..." << endl << flush;

	// Perform the scattering.

	intermediate_state3 inter;
	final_state3 final;
	real cpu_scatter = cpu_time();

	int result = single_scatter(init, inter, final,
				    cpu_time_check, dt_snap, snap_cube_size);
	cpu_scatter = cpu_time() - cpu_scatter;

	// cerr << "...back" << endl;
	// cerr << "cpu time = " << cpu_scatter << endl << flush;

#if FORCE_TIMEOUT

	// Temporary only...

	if (cpu_scatter > TIME_TIMEOUT) {

	    // Force a timeout: waste a lot of time...

	    real x = 0;
	    for (int j = 0; j < 1000000; j++)
	        for (int i = 0; i < 1000000000; i++)
		    x += sqrt(0.1*i + j);

	    cpu_scatter = x;

	}

#endif

	// Return the results.

	pvm_initsend(PvmDataRaw);

	pvm_pkint(&n_rand, 1, 1);

	pack(init);
	pack(inter);
	pack(final);

	pvm_pkint(&result, 1, 1);
	pvm_pkdouble(&cpu_scatter, 1, 1);

	// Send the data as message RETURN_DATA_MSG.

	if (pvm_send(pvm_parent(), RETURN_DATA_MSG) < 0) {
	    cerr << "Error returning data from " << name << endl;
	    cerr << "CPU time = " << cpu_scatter << endl;
	}

    }
}

#else

#include "stdinc.h"

void main()
{
    err_exit("PVM not available");
}

#endif
