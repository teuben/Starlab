
#include "sigma3.h"
#include "pvm_scatt.h"

#ifdef HAS_PVM

// Note: ALL the PVM code for sigma3 is contained in this file.

#include "pvm3.h"

//----------------------------------------------------------------------------

// Externally visible functions (defined at end):

void initialize_processors(int ntask, int debug);
void terminate_processors(int debug);
int multiscatter3(scatter_profile & prof, sigma_out & out,
		  real rho_sq_min, real rho_sq_max, int rho_zone,
		  real dt_snap, real snap_cube_size,
		  real cpu_time_check, real cpu_init, real &cpu_save,
		  int& scatt_total, real& cpu_total, stat_fp acc_stats,
		  int debug, int scatter_summary_flag);

//----------------------------------------------------------------------------

#define RECV_TIMEOUT	3600		// time out on receive after 1 hour

// Global variables (note capitalization):

static char Exe[128];			// slave process name

#define NTASK_MAX	16		// maximum number of tasks to start

static int  Tid[NTASK_MAX];		// TID of each task
static int  Nscatt[NTASK_MAX];		// # scatterings performed by each task
static real Cpu[NTASK_MAX];		// total CPU time consumed by each task
static char Name[NTASK_MAX][128];	// host name of each task

static int  Random_seed[NTASK_MAX];	// last random seed used for each task
static initial_state3 Init[NTASK_MAX];	// last init sent to each task

static int  N_task_init = 0;		// task # runs from 0 to N_task_init-1

//----------------------------------------------------------------------------

local int task_num(int task_id)
{
    for (int task = 0; task < N_task_init; task++)
	if (Tid[task] == task_id) return task;
    return -1;
}

local char* host_name(int task_id)
{
    int task = task_num(task_id);
    if (task >= 0)
	return Name[task];
    else
	return NULL;
}

local void initialize_task(int task, int reinit, int io)
{
    if (reinit) {

      // Start a new task on a specific processor.

      if ( pvm_spawn(Exe, (char**)0, 1,
		     Name[task], 1, Tid+task) != 1) {
	  cerr << "Unable to re-spawn task #" << task << endl;
	  exit(1);
      }

    } else {

        // Start a brand new task on any available processor.
      
      if ( pvm_spawn(Exe, (char**)0, 33, ".", 1, Tid+task) != 1) {
	  cerr << "Unable to spawn task #" << task << endl;
	  exit(1);
      }

    }

    if (pvm_recv(Tid[task], HANDSHAKE_MSG) < 0) {
        cerr << "Error getting new ID for task #" << task << endl;
	exit(1);
    }

    pvm_upkstr(Name[task]);

    if (io) {
        cerr << (reinit ? "Rei" : "I");
        cerr << "nitialized PVM task #" << task
	     <<",  TID = " << Tid[task]
	     << ",  node = " << Name[task] << endl << flush;
    }

    // Seek notification on task exit (including fatal error).

    pvm_notify(PvmTaskExit, FAILURE_MSG, 1, Tid+task);

}

local void initialize_scattering(scatter_profile & prof,
				 real rho_sq_min, real rho_sq_max,
				 int tid, int scatter_summary,
				 real cpu_time_check,
				 real dt_snap, real snap_cube_size)
{
    int task = task_num(tid);
    single_scatter_init(prof, rho_sq_min, rho_sq_max,
			Init[task], Random_seed[task],
			scatter_summary, dt_snap, snap_cube_size);

    // Pack the data.

    pvm_initsend(PvmDataRaw);

    pvm_pkint(Random_seed+task, 1, 1);
    pack(Init[task]);

    pvm_pkdouble(&cpu_time_check, 1, 1);
    pvm_pkdouble(&dt_snap, 1, 1);
    pvm_pkdouble(&snap_cube_size, 1, 1);

    // Send the data as message SEND_DATA_MSG.

    pvm_send(tid, SEND_DATA_MSG);

    if (scatter_summary > 0)
	cerr << "task " << task << ":  sent data to " << tid
	     << " (" << host_name(tid) << ")\n" << flush;
    
}

local void closeout_tasks(scatter_profile & prof,
			  sigma_out & out, int rho_zone,
			  real dt_snap, real snap_cube_size,
			  stat_fp acc_stats, int debug,
			  int n_sent, int & n_complete,
			  intermediate_state3 & inter,
			  final_state3 & final)
{
    if (n_sent < out.trials_per_zone) {

	terminate_processors(debug);
	pvm_exit();
	cerr << "Unable to continue -- cannot identify failed task\n";
	err_exit("multiscatter3: calculation aborted.");

    } else {

	// Close out all remaining scatterings and add the
	// appropriate number of errors into the statistics.

	// Note: only the structure elements listed here are
	// necessary in the case of an error.

	inter.descriptor = unknown_intermediate;
	inter.n_osc = 0;

	final.descriptor = error;
	final.time = 0;
	final.error = 0;
	final.n_steps = 0;

	cerr << n_sent - n_complete << " uncompleted tasks:\n";

	for (int task = 0; task < N_task_init; task++)

	    if (Random_seed[task] > 0) {

		summarize_scattering_initial(Init[task],
					     Random_seed[task],
					     dt_snap, snap_cube_size);

		// Close out this task by accumulating a (fake)
		// error status in the statistics.

		single_scatter_stats(prof, Init[task], inter, final,
				     rho_zone,
				     out, acc_stats, 0);

		Random_seed[task] = 0;
		n_complete++;

		// The task still hasn't completed, so it will not be
		// usable for further scatterings.  Kill and restart it here.

		// Keep the old task number, cpu, and nscatt.
		// Modify only the TID.

		// Note that the TID will apparently be garbled after
		// we kill the task, so we can't verify the source of a
		// FAILURE_MSG signal.  Solve this problem by identifying
		// the notification with a different message ID.

		pvm_notify(PvmTaskExit, KILL_CONFIRM, 1, Tid+task);

		if (pvm_kill(Tid[task]) < 0) {

		    cerr << "Unable to kill task #" << task << endl;
		    exit(1);

		} else {

		    cerr << "Task #" << task << " ("
		         << Tid[task] << ", node " << Name[task]
			 << ") killed..." << flush;

		}

		// Wait for confirmation of the kill.

		struct timeval timeout;
		timeout.tv_sec = 600;		// Give up after 10 minutes
		int bufid;

		// Return values for bufid:	< 0 ==> error
		//				= 0 ==> timeout?      <--- ???
		//				> 0 ==> success

		if ( (bufid = pvm_trecv(-1, KILL_CONFIRM, &timeout)) <= 0 ) {

	 	    if (bufid < 0)
		        cerr << "error receiving confirmation.\n";
		    else
		        cerr << "timeout receiving confirmation.\n";

		    exit(1);

		} else {

		    int tid = -1;

		    if (pvm_upkint(&tid, 1, 1) < 0) {
		        cerr << "error unpacking confirmation message.\n";
			exit(1);
		    } else
		        cerr << "confirmation received from TID "
			     << tid << ".\n";

		}

		// Wait for 10 us before continuing.

		starlab_wait(10);

		// Flush spurious messages (apparently coming from
		// deleted tasks).

		int n_flush = 0;
		while ((bufid = pvm_probe(-1, FAILURE_MSG)) > 0) {

		    pvm_recv(-1, FAILURE_MSG);

		    int tid = -1;
		    if (pvm_upkint(&tid, 1, 1) == 0) {
		        if (tid != Tid[task]) {
			    cerr << "..." << tid;
			    n_flush++;
			}
		    }

		}

		if (n_flush > 0)
		    cerr << ":  " << n_flush << " unexpected failure message"
		         << (n_flush > 1 ? "s" : "")
			 << " flushed\n";

		// Start up a new task.

		initialize_task(task, 1, 1);

	    }

	if (n_complete != n_sent)
	    cerr << "warning: " << n_sent - n_complete
		<< " tasks unaccounted for...\n";
    }
}

//----------------------------------------------------------------------------

// Externally visible functions:

void initialize_processors(int ntask, int debug)
{
    // Initialize up to NTASK_MAX scattering processes.
    // Should really determine the number of processors available.

    N_task_init = min(NTASK_MAX, ntask);

    if (getenv("PVM_ROOT") == NULL) err_exit("PVM not available"); // Run time
								   // check...

    strcpy(Exe, getenv("STARLAB_PATH"));	// ASSUME that these
    if (Exe[0] == '\0')				// environment variables
	strcpy(Exe, getenv("PWD"));		// are defined
    else
	strcat(Exe, "/bin");
    strcat(Exe, "/scatter3_slave.pvm");

    // cerr << "slave executable name is " << Exe << endl;

    for (int task = 0; task < N_task_init; task++) {

        initialize_task(task, 0, 0);

	Cpu[task] = 0;
	Nscatt[task] = 0;

    }

    // On exit, we have established N_task_init processes, which are now
    // waiting for data.

    cerr << N_task_init << " PVM slave processes initialized\n";

    if (debug) {
	cerr << "\nhost names:";
	for (int i = 0; i < N_task_init; i += 4) {
	    for (int task = i; task < min(i+4, N_task_init); task++)
		cerr << "  " << Name[task];
	    cerr << endl;
	    if (i < N_task_init - 4) cerr << "           ";
	}

	cerr << "\nhost TIDs:";
	for (int i = 0; i < N_task_init; i += 4) {
	    for (int task = i; task < min(i+4, N_task_init); task++)
		cerr << "  " << Tid[task];
	    cerr << endl;
	    if (i < N_task_init - 4) cerr << "          ";
	}
    }

}

void terminate_processors(int debug)
{
    for (int task = 0; task < N_task_init; task++) {

	pvm_initsend(PvmDataRaw);

	// Send message TERMINATE_MSG to terminate.

	if (pvm_send(Tid[task], TERMINATE_MSG) < 0)
	    cerr << "Error sending termination message to task = "
		 << Tid[task] << "  (" << Name[task] << ")\n";

    }

    cerr << N_task_init << " PVM slave processes terminated\n";

    if (debug) {

	// Print use and load information.

	cerr << "\nscatterings per task:  ";
	for (int i = 0; i < N_task_init; i += 4) {
	    for (int task = i; task < min(i+4, N_task_init); task++)
		cerr << "  " << Nscatt[task];
	    cerr << endl;
	    if (i < N_task_init-4) cerr << "                       ";
	}

	cerr << "\nCPU %load distribution:";
	int p = cerr.precision(4);		// nonstandard precision

	real total = 0;
	for (int task = 0; task < N_task_init; task++) total += Cpu[task];
	for (int i = 0; i < N_task_init; i += 4) {
	    for (int task = i; task < min(i+4, N_task_init); task++)
		cerr << "  " << 100*Cpu[task]/total;
	    cerr << endl;
	    if (i < N_task_init-4) cerr << "                       ";
	}
	cerr.precision(p);
    }

    N_task_init = 0;
    pvm_exit();
}

// multiscatter3_pvm: Perform a specified number of scattering experiments in
//		      a given impact parameter range, accumulating statistics
//		      as we go.  Return the total number of hits in this range.
//
//		      Distribute the task over several processors using PVM.
//
// NOTE that this function is called multiscatter3, the same as the non-PVM
// version.  It is designed to replace that version in its entirety...

int multiscatter3(scatter_profile & prof, sigma_out & out,
		  real rho_sq_min, real rho_sq_max, int rho_zone,
		  real dt_snap, real snap_cube_size,
		  real cpu_time_check, real cpu_init, real &cpu_save,
		  int& scatt_total, real& cpu_total, stat_fp acc_stats,
		  int debug, int scatter_summary_flag)
{
    int n_sent = 0, n_complete = 0;
    int total_hits = 0;

    if (N_task_init <= 0) err_exit("PVM not initialized.");

    // Start N_task_init scatterings.

    for (int task = 0; task < N_task_init; task++) {

	// Get the next scattering setup.

	// The function single_scatter_init will randomize the angles
	// in the same way as the standalone tool scatter3.

	initialize_scattering(prof, rho_sq_min, rho_sq_max, Tid[task],
			      scatter_summary_flag, cpu_time_check,
			      dt_snap, snap_cube_size);
	n_sent++;
    }

    // Now wait for data to come back and send new data as needed.

    int single_result;
    real cpu_scatter;

    while (n_complete < out.trials_per_zone) {

	int n_rand;
	initial_state3 init;
	intermediate_state3 inter;
	final_state3 final;

	int scatter_summary = scatter_summary_flag;

	struct timeval timeout;
	timeout.tv_sec = RECV_TIMEOUT;
	int bufid, length, msgid, tid_r;

	int task;

	// Return values for bufid:	< 0 ==> error
	//				= 0 ==> timeout?	<--- ???
	//				> 0 ==> success

	if ( (bufid = pvm_trecv(-1, -1, &timeout)) <= 0 ) {

	    // A timeout or a system error has occurred.  Go straight to
	    // the output phase if all data have been sent (i.e. don't wait
	    // any longer for the remaining data to return).  Otherwise,
	    // quit with a error.

	    if (bufid < 0)
		cerr << "\nError receiving data...\n";
	    else
		cerr << "\nTimeout receiving data...\n";

	    // We cannot identify the failed task -- do the best we can.

	    closeout_tasks(prof, out, rho_zone,
			   dt_snap, snap_cube_size,
			   acc_stats, debug,
			   n_sent, n_complete, inter, final);
	    break;

	}

	length = msgid = tid_r = -1;
	int status = pvm_bufinfo(bufid, &length, &msgid, &tid_r);
	bool pvm_error = (status < 0 || length <= 0 || msgid < 0 || tid_r < 0);

	if (pvm_error) {

	    cerr << "\nUnknown or illegal message received...\n";

	    // Basic diagnostics:

	    PRL(bufid);
	    PRL(status);
	    PRL(length);
	    PRL(msgid);
	    PRL(tid_r);

	    PRL(n_sent);
	    PRL(n_complete);

	    closeout_tasks(prof, out, rho_zone,
			   dt_snap, snap_cube_size,
			   acc_stats, debug,
			   n_sent, n_complete, inter, final);
	    break;

	}

	if (msgid == RETURN_DATA_MSG || msgid == FAILURE_MSG) {

	    // Check for the failure of a child process.

	    if (msgid == FAILURE_MSG) {

		// A slave process has unexpectedly exited.  Print a message
		// and start a new process.

		pvm_upkint(&tid_r, 1, 1);

		cerr << "\nTask " << tid_r << " (" << host_name(tid_r)
		     << ") has terminated prematurely.\n";

		task = task_num(tid_r);

		// Recover the initial conditions for use below.

		n_rand = Random_seed[task];
		init = Init[task];

		// Start a new slave process.  (Keep old task number,
		// cpu, and nscatt.  Modify only TID.)

		initialize_task(task, 1, 1);

		// Do NOT start a new scattering.  Instead, create a fake
		// "error" return state.  (Only the structure elements
		// listed here are necessary in the case of an error.)

		inter.descriptor = unknown_intermediate;
		inter.n_osc = 0;

		final.descriptor = error;
		final.time = 0;
		final.error = 0;
		final.n_steps = 0;

		single_result = 0;
		cpu_scatter = 0;

		pvm_error = true;

	    } else {

		// Normal completion of a slave subtask.

		task = task_num(tid_r);

		if (scatter_summary > 0)
		    cerr << "\nreceived data from TID = " << tid_r
			 << " (" << Name[task] << ")\n" << flush;

		pvm_upkint(&n_rand, 1, 1);

		// Note: n_rand should be Random_seed[task],
		//	 init should be Init[task].

		unpack(init);
		unpack(inter);
		unpack(final);

		pvm_upkint(&single_result, 1, 1);
		pvm_upkdouble(&cpu_scatter, 1, 1);

	    }

	    n_complete++;
	    total_hits += single_result;
	    scatt_total++;
	    cpu_total += cpu_scatter;

	    Nscatt[task]++;
	    Cpu[task] += cpu_scatter;

	    // Diagnostics:

	    if (abs(debug) > 2) {
		int p = cerr.precision(STD_PRECISION);
		cerr << "single_scatter: rho_max^2 = " << rho_sq_max
		     << " ; returning with n_hit = " << single_result
		     << " and n_hit_tot = " << out.n_hit_tot + single_result
		     << endl;
		cerr.precision(p);
	    }

	    if (final.descriptor == error) {

		// Print out enough information to repeat the calculation.

		summarize_scattering_initial(init, n_rand,
					     dt_snap, snap_cube_size);
		scatter_summary = 2;
		if (pvm_error) cerr << "final state unknown\n";

	    }

	    if (!pvm_error && scatter_summary > 0)
		summarize_scattering_final(inter, final,
					   scatter_summary, cpu_scatter);

	    // Start a new process, if necessary.

	    if (n_sent < out.trials_per_zone) {

		initialize_scattering(prof, rho_sq_min, rho_sq_max, tid_r,
				      scatter_summary_flag, cpu_time_check,
				      dt_snap, snap_cube_size);
		n_sent++;

	    } else {

		// Don't terminate any processes here -- handled elsewhere.
		// Do set Random_seed = 0 to signify that this task is done.

		Random_seed[task] = 0;
	    }

	} else {

	    cerr << "unexpected PVM message ID = "
		 << msgid << " received" << endl;

	    err_exit("multiscatter3: fatal error");

	}

	// Accumulate statistics on the result.

	single_scatter_stats(prof, init, inter, final, rho_zone,
			     out, acc_stats, single_result);

    }

    return total_hits;

}

#endif
