
/* Duffing oscillator in ONE dimension -- Leapfrog integrator. */

#include <stdio.h>
#include <math.h>

/*======================================================================*/

/* Global constants: */

static double ALPHA = 0.1;	/* Damping */
static double BETA = 0.0;

static double AMPL = 9.8;	/* Forcing amplitude */
static double OMEGA_D = 1.0;	/* Forcing frequency */
static double PHASE = 0.0;	/* Forcing phase */


/* Define accelerations and potential (1-D ONLY).
   Allow t as an argument in case of time-dependent forces. */

double accx_cons(double x, double v, double t)
{
    return -x*x*x + AMPL*cos(OMEGA_D*t + PHASE);
}

double accx_diss(double x, double v, double t)
{
    return -v * (ALPHA + BETA * abs(v));
}

double accx(double x, double v, double t)
{
    return accx_cons(x, v, t) + accx_diss(x, v, t);
}

double potential(double x)	/* Must be CONSISTENT with accx_cons()! */
{
    return -0.25*x*x*x*x;
}

/*======================================================================*/

/* Take a single time step, NEGLECTING the work integral (for now). */

static int init = 0;	/* Initialization flag. */

void step(double* t, double* x, double* v, double* w, double dt)
{
    double vv;

    if (init == 0) {

	/* Offset the velocity by half a step the first time through. */

	*v += 0.5*dt*accx(*x, *v, *t);
	init = 1;
    }

    /* Update all quantities using the Leapfrog scheme. */

    *x += (*v)*dt;
    *t += dt;

    /* Temporarily synchronize v with x (call it vv) for use with accx. */

    vv = *v + 0.5*dt*accx(*x, *v, *t);
    *v += accx(*x, vv, *t)*dt;
}

double energy(double x, double v, double t, double dt)
{
    /* Must synchronize v in this case (==> dt and t needed)... */

    v -= -0.5*dt*accx(x, v, t);		/* (Using a local copy of v, note.) */

    return 0.5*v*v + potential(x);
}

/*======================================================================*/

double interp(double t0, double t1, double x0, double x1, double t)
{
    /* Linearly interpolate x to its value at time t */

    return x0 + (x1 - x0) * (t - t0) / (t1 - t0);
}

main()
{
    /* Declare and initialize variables. */

    double x, v;
    double t, dt = 0.025, t_max = 500.0;
    char command[40] = "                    ";
    double w;
    char poincare = 0;
    double t_next;
    double t_start = 100.0;
    double x0 = 1.0, v0 = 0.0;
    char clear = 1;

    /* Use dialog boxes to set essential parameters. */

    create_double_dialog_entry("x0", &x0), set_dialog_width(16);
    create_double_dialog_entry("v0", &v0), set_dialog_width(16);
    create_double_dialog_entry("ALPHA", &ALPHA);
    create_double_dialog_entry("AMPL", &AMPL);
    create_double_dialog_entry("OMEGA_D", &OMEGA_D);
    create_double_dialog_entry("PHASE", &PHASE);
    create_double_dialog_entry("dt", &dt);
    create_double_dialog_entry("t_max", &t_max);
    create_double_dialog_entry("t_start", &t_start);
    create_button_dialog_entry("Poincare", &poincare);
    create_string_dialog_entry("command", command);
    create_button_dialog_entry("clear command", &clear);

    set_up_dialog("Duffing Oscillator", 0, 0);

    while (read_dialog_window()) {

	/* Output graphics control command(s), if any. */

	printf("%s\n", command);

	/* Treat "-e" as a special case */

	if (strcmp(command, "-e") && (!poincare || t_max > t_start)) {

	    if (poincare) printf("-pp\n");

	    /* Set initial conditions (start at rest at x = 0): */

	    init = 0;
	    t = 0.0;
	    x = x0;
	    if (AMPL == 0 && x == 0.0) x = 1.0;
	    v = v0;

	    t_next = 0.0;

	    /* Note that w is never updated in this version of the program. */

	    w = 0.0;

	    /* Initial output */

	    printf("%f %f %f\n", t, x, v);

	    /* Calculate and print out the trajectory. */

	    while (t < t_max) {

		double t_last = t;
		double x_last = x;
		double v_last = v;

		step(&t, &x, &v, &w, dt);

		if (!poincare)

		    printf("%f %f %f\n", t, x, v);

		else {

		    if (t > t_next) {

			if (t >= t_start) 
			    printf("%f %f %f\n", t_next,
				   interp(t_last, t, x_last, x, t_next),
				   interp(t_last, t, v_last, v, t_next));

			t_next += 2*M_PI/OMEGA_D;    /* Sample at the forcing
							frequency */

		    }
		}		
	    }

	    /* Signify end-of-data to force plot_data to plot the points. */

	    printf("end of data\n");
	}

	fflush(stdout);

	if (clear) command[0] = '\0';
	reset_dialog_window();

    }
}
