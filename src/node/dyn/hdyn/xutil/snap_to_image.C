
//// snap_to_image:  Construct images (in Sun rasterfile format) of
////                 a series of snapshots.
////
//// Options:   -f filename  specify root name of image files [snap, - = stdout]
////            -l scale     specify size of field of view (+/- scale) [2]
////            -m           use mass to determine star size [no]
////            -n nmax      specify maximum number of images to produce [-]
////            -p psize     specify star size, in pixels [3]
////            -P axis      specify projection axis [z]
////            -s size      specify image size, in pixels [256]
////            -S nskip     specify snaps to skip between images [0]

//.............................................................................
//
//    version 1:  Nov 1998   Steve McMillan	 email: steve@zonker.drexel.edu
//			     Drexel University, Philadelphia, PA, USA
//.............................................................................

#include "hdyn.h"

void write_image(float* a, int m, int n, char* filename, int scale);

#define IOFF	10
#define CMAX	256

#define L	2.0
#define NX	256
#define NY	256

void add_point(float* a, int nx, int ny, int i, int j, real r, float color)
{
    // Stars are square for r < 2, have corners removed for r = 2,
    // are approximate circles for r > 2.

    int ir = (int)r;

    for (int ii = max(-ir, -i); ii <= min(ir, nx-i); ii++)
	for (int jj = max(-ir, -j); jj <= min(ir, ny-j); jj++)
	    if (ir < 2
//		|| (ir == 2 && (abs(ii) < ir || abs(jj) < ir))
		|| (ii*ii + jj*jj <= (ir+0.5)*(ir+0.5)))
		*(a+(j+jj)*nx+(i+ii)) = color;
}

main(int argc, char** argv)
{
    int count = 0, count1 = 0;
    char filename[64], command[64];
    char* fn;

    real l = L;
    int nx = NX, ny = NY;
    int n = 0, nskip = 0;
    int psize = 3;
    int axis = 3;		// 1 = x, 2 = y, 3 = z

    char file[64];
    strcpy(file, "snap");

    bool mass = false;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "f:l:mn:p:P:s:S:";

    while ((c = pgetopt(argc, argv, param_string)) != -1) {
	switch (c) {
	    case 'f':	strncpy(file, poptarg, 63);
			file[63] = '\0';	// just in case
			break;
	    case 'l':	l = atof(poptarg);
			break;
	    case 'm':	mass = true;
			break;
	    case 'n':	n = atoi(poptarg);
	    		break;
	    case 'p':   psize = atoi(poptarg);
	    		if (psize < 1) psize = 1;
		        break;
	    case 'P':	if (poptarg[0] == 'x' || atoi(poptarg) == 1)
			    axis = 1;
	    		else if (poptarg[0] == 'y' || atoi(poptarg) == 2)
			    axis = 2;
	    		else if (poptarg[0] == 'z' || atoi(poptarg) == 3)
			    axis = 3;
			break;
	    case 's':   nx = ny = atoi(poptarg);
		        break;
	    case 'S':   nskip = atoi(poptarg);
		        break;
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }

    float* a = new float[nx*ny];
    if (!a) exit(1);

    int iax, jax;

    if (axis == 1) {
	iax = 1;
	jax = 2;
    } else if (axis == 2) {
	iax = 2;
	jax = 0;
    } else {
	iax = 0;
	jax = 1;
    }

    real p2 = 0.5*psize;

    int cmin = 1000000, cmax = -1000000;
    real color_scale;

    real mmin = VERY_LARGE_NUMBER, mmax = -VERY_LARGE_NUMBER;
    real logfac = 1;

    // Loop over input snapshots.

    hdyn* b;
    while (b = get_hdyn(cin)) {

	if (count == 0) {

	    // Determine overall color scaling (and mass range, if
	    // relevant) from first snap.

	    for_all_daughters(hdyn, b, bb) {
		if (bb->get_index() >= 0) {
		    cmin = min(cmin, bb->get_index()+IOFF);
		    cmax = max(cmax, bb->get_index()+IOFF);
		}
		if (mass && bb->get_mass() > 0) {
		    mmin = min(mmin, bb->get_mass());
		    mmax = max(mmax, bb->get_mass());
		}
	    }

	    if (cmax >= CMAX) cmax = CMAX - 1;
	    color_scale = 1.0 / cmax;

	    if (mass) logfac = 1.0/log10(mmax/mmin);
	}

	if (nskip <= 0 || count % (nskip+1) == 0) {

	    // Make this snapshot into an image.

	    for (int i = 0; i < nx*ny; i++) a[i] = 0;

	    for_all_daughters(hdyn, b, bb) {

		// Should probably sort stars by coordinate along the
		// projection axis, as in xstarplot.  NOT done yet.

		vector x = bb->get_pos();

		if (x[iax] > -l && x[iax] < l
		    && x[jax] > -l && x[jax] < l) {

		    // Coordinates:

		    int i = (int) ((l+x[iax]) * 0.5 * nx / l);
		    int j = (int) ((l+x[jax]) * 0.5 * ny / l);

		    // Color (by index):

		    float color = 0.9999;		// default is white

		    if (bb->get_index() > 0)
			color = (bb->get_index() % (cmax - cmin) + IOFF)
			    		* color_scale;

		    // Radius of point representing star:

		    real r = p2;

		    if (mass && bb->get_mass() > 0) {

			// Star size depends on its mass.  Minimum size is
			// 1 pixel, maximum is p2.  Scaling is logarithmic.

			r *= log10(bb->get_mass()/mmin) * logfac;
		    }

		    add_point(a, nx, ny, i, j, r, color);
		}
	    }

	    // Image file name counts output images, not input snaps.

	    if (streq(file, "-"))
		fn = NULL;
	    else {
		sprintf(filename, "%s.%3.3d.sun", file, count1++);
		fn = filename;
	    }

	    write_image(a, nx, ny, fn, 0);	// 0 ==> don't scale

	    if (fn) {

		// Runtime compression:

		sprintf(command, "gzip -f -q %s &", filename);
		system(command);
	    }
	}

	rmtree(b);

	count++;
	if (count%10 == 0) cerr << count1 << "/" << count << " ";

	if (n > 0 && count1 > n) break;
    }
    cerr << endl;
}
