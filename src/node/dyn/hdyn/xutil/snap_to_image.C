
//// snap_to_image:  Construct images (in Sun rasterfile format) of
////                 a series of snapshots.
////
//// Options:   -1           toggle combine all frames in a single image   [yes]
////            -c           compress the image file(s) using gzip          [no]
////            -C colormap  specify a colormap file name                   [no]
////            -f filename  specify root name of image files   ["-" --> stdout]
////            -F format    specify image file format
////                             (0 = PNG, 1 = SUN, 2 = GIF)                 [0]
////            -g           write GIF files (currently uses convert and
////                             gifsicle; slow and messy, and may become
////                             uninterruptible, but saves disk space)
////                             -- same as "-F 2"                          [no]
////            -G           toggle forcing particles to grid (nicer single
////                             frames, but jerkier movies)              [true]
////            -H           toggle Herzsprung-Russel diagram or positional
////                             plot                                 [position]
////            -i index     specify (real) color index for all stars
////                                                        [use internal index]
////            -l scale     specify size of field of view (+/- scale)       [3]
////            -m           use mass to determine star color and/or size   [no]
////            -n nmax      specify maximum number of images to produce   [Inf]
////            -N nbody     color using a (small-N) colormap               [no]
////            -p psize     specify star radius, in pixels
////                                           [0 (single image), 1 (animation)]
////            -P axis      specify projection axis                         [z]
////            -q           toggle suppression of diagnostic output
////                                                            [don't suppress]
////            -r           use stellar radius to set point size           [no]
////            -s size      specify image size, in pixels                 [256]
////            -S nskip     specify snaps to skip between images            [0]
////            -t           test the color map [don't test]
////            -x           specify right (log effective temparature) edge
////                             of HRD (-H only)                            [3]
////            -X           specify left (log effective temparature) edge
////                             of HRD (-H only)                            [5]
////            -y           specify minimum (log luminosity/Lsun) limit
////                             of HRD (-H only)                           [-3]
////            -Y           specify maximum (log luminosity/Lsun) limit
////                             of HRD (-H only)                            [3]
////
//// Note: If PNG output is requested and the PNG libraries are unavailable,
////       GIF output is used instead.
//.............................................................................
//
//    version 1:  Nov 1998   Steve McMillan	 email: steve@zonker.drexel.edu
//			     Drexel University, Philadelphia, PA, USA
//    version 2:  Aug 2002   Steve McMillan
//.............................................................................

#include "hdyn.h"
#include "star/single_star.h"

// From gfx/util/write_image.C:

void write_image(float* a, int m, int n, char* filename, int scale,
		 char *colormap_file);
void write_image(float* a, int m, int n, char* filename, int scale,
		 unsigned char *red,
		 unsigned char *green,
		 unsigned char *blue);

#include "write_png.h"

// From gfx/util/make_header.C:

void make_standard_colormap(unsigned char* red,
			    unsigned char* green,
			    unsigned char* blue);

void make_greymap(unsigned char* red,
		  unsigned char* green,
		  unsigned char* blue);

// Our three-color concention now requires us to define our own standards,
// basically the same map replicated three times at different intensities.
// Thus, index, index-85, index-170, and 0 are the 4 color levels (including
// black) associated with a given color.

#define REDUCE 1

void extend_local_colormap(unsigned char* red,
			   unsigned char* green,
			   unsigned char* blue,
			   int psize)
{
    // Fill in the remaining entries.

    red[0] = green[0] = blue[0] = 0;

    for (int i = 1; i <= 85; i++) {

	real fac = 0.65;
	if (psize == 2) fac = 0.7;	// a kludge!
	fac *= REDUCE;

	red[i] = (unsigned char)(fac*red[i+170]);
	green[i] = (unsigned char)(fac*green[i+170]);
	blue[i] = (unsigned char)(fac*blue[i+170]);

	fac = 0.85;
	fac *= REDUCE;

	red[i+85] = (unsigned char)(fac*red[i+170]);
	green[i+85] = (unsigned char)(fac*green[i+170]);
	blue[i+85] = (unsigned char)(fac*blue[i+170]);
    }
}

local void make_local_standard_colormap(unsigned char* red,
					unsigned char* green,
					unsigned char* blue,
					int psize)
{
    int i;

    // Indices 171-255 are the same as the standard 256-color map.
    // Do them first, then derive the others.  Tricky...!

    for (i = 0; i < 10; i++) {
        red[170+i] = 0;
        green[170+i] = 0;
        blue[170+i] = 28*i;
    }
    for (i = 10; i < 32; i++) {
        red[170+i] = 0;
        green[170+i] = 12*(i-10);
        blue[170+i] = 255;
    }
    for (i = 32; i < 53; i++) {
        red[170+i] = 12*(i-32);
        green[170+i] = 255;
        blue[170+i] = 255 - red[170+i];
    }
    for (i = 53; i < 75; i++) {
        red[170+i] = 255;
        green[170+i] = 255 - 12*(i-53);
        blue[170+i] = 0;
    }
    for (i = 75; i < 85; i++) {
        red[170+i] = 255;
        blue[170+i] = green[170+i] = 28*(i-75);
    }

    red[255] = green[255] = blue[255] = 255;
    extend_local_colormap(red, green, blue, psize);
}

void make_alternate_colormap(unsigned char* red,
			     unsigned char* green,
			     unsigned char* blue);

local void make_local_alternate_colormap(unsigned char* red,
					 unsigned char* green,
					 unsigned char* blue,
					 int psize)
{
    make_alternate_colormap(red, green, blue);

    // Compress:

    for (int i = 1; i < 85; i++) {
	red[255-i] = red[255-3*i];
	green[255-i] = green[255-3*i];
	blue[255-i] = blue[255-3*i];
    }

    // Extend:

    extend_local_colormap(red, green, blue, psize);
}

local void make_local_greymap(unsigned char* red,
			      unsigned char* green,
			      unsigned char* blue,
			      int psize)
{
    for (int i = 1; i <= 85; i++)
        red[170+i] = green[170+i] = blue[170+i] = 3*i;

    extend_local_colormap(red, green, blue, psize);
}

local void make_local_small_n_colormap(unsigned char* red,
				       unsigned char* green,
				       unsigned char* blue,
				       int ncolor,
				       int psize)
{
    if (ncolor != 3 && ncolor != 4) return;

    if (ncolor == 3) {		// RGB

	for (int i = 1; i <= 28; i++) {
	    red[170+i] = 255;
	    green[170+i] = blue[170+i] = 0;
	}
	for (int i = 29; i <= 56; i++) {
	    green[170+i] = 255;
	    red[170+i] = blue[170+i] = 0;
	}
	for (int i = 57; i <= 85; i++) {
	    blue[170+i] = 255;
	    red[170+i] = green[170+i] = 0;
	}
	
    } else {			// RGBW

	for (int i = 1; i <= 21; i++) {
	    red[170+i] = 255;
	    green[170+i] = blue[170+i] = 0;
	}
	for (int i = 22; i <= 42; i++) {
	    green[170+i] = 255;
	    red[170+i] = blue[170+i] = 0;
	}
	for (int i = 43; i <= 64; i++) {
	    blue[170+i] = 255;
	    red[170+i] = green[170+i] = 0;
	}
	for (int i = 65; i <= 85; i++) {
	    red[170+i] = green[170+i] = blue[170+i] = 255;
	}
	
    }

    extend_local_colormap(red, green, blue, psize);
}

local void add_point(float *a, int nx, int ny,
		     real x, real y, int i, int j,
		     bool grid,
		     real r, float color,
		     real z, float *zarray)
{
    // The actual position of the particle in pixel coordinates
    // is (x,y).  The central pixel is (i,j).  We expect i <= x,
    // j <= y.  The desired radius is r.

    // Despite the fact that almost everything to do with colors
    // deals with integer indices, the "color" used here is real,
    // with the range 0-1, mapping linearly to indices 0-255...
    // With our new scheme, the color specified should actually be
    // at the bright end of the colormap, in the range 0.667-1.0.
    // May get around to cleaning up the color specification...

    real xref = x;
    real yref = y;
    if (grid) {
	xref = i+0.5;
	yref = j+0.5;
    }

    int ir = (int)(r+0.0001)+1;
    real r2a = pow(r-0.1,2), r2b = pow(r+0.2,2), r2c = pow(r+0.5,2);

    // PRC(x); PRC(i); PRC(y); PRL(j);

    for (int ii = max(-ir, -i); ii <= min(ir, nx-i-1); ii++)
	for (int jj = max(-ir, -j); jj <= min(ir, ny-j-1); jj++) {

	    // Get distance from pixel center to reference point.

	    real xi = i+ii+0.5;
	    real yi = j+jj+0.5;

	    real dx = abs(xi-xref);
	    real dy = abs(yi-yref);

	    real dx2 = pow(xi-xref,2) + pow(yi-yref,2);

	    // Soften the edges of the stars.  (Details fine-tuned by eye.)
	    // Use 3 color levels, in addition to black.

	    real c = color;

	    if ((ii == 0 && jj == 0) || dx2 <= r2a)
		; 				// "color 1" -- do nothing
	    else if (dx2 < r2b)
		c -= 0.3333;			// "color 2" (index - 85)
	    else if (dx2 < r2c)
		c -= 0.6667;			// "color 3" (index - 170)
	    else
		c = 0;				// black

	    if (c > 0) {

//		cerr << "adding " << c << " to ("
//		     << i+ii << "," << j+jj << ")" << endl;

		real zcur = *(zarray+(j+jj)*nx+(i+ii));
		if (z > zcur) {
		    *(a+(j+jj)*nx+(i+ii)) = c;
		    *(zarray+(j+jj)*nx+(i+ii)) = z;
		}
	    }
	}
}

local void initialize_arrays(float *a, float *zarray, int n)
{
    for (int i = 0; i < n; i++) {
	a[i] = 0;
	zarray[i] = -1.e10;
    }
}

local void write_png_image(float *a, int nx, int ny,
			   char *file_name,
			   unsigned char *red,
			   unsigned char *green,
			   unsigned char *blue)
{
    // Input a is in the range [0, 1).  Convert to unsigned char.

    unsigned char *image = (unsigned char *)malloc(nx*ny*sizeof(unsigned char));
    for (int j = 0; j < ny; j++)
	for (int i = 0; i < nx; i++)
	    *(image + nx*j + i) = (unsigned char) (255.9 * (*(a + nx*j + i)));
#ifdef HAS_PNG
    write_png(image, nx, ny, file_name, red, green, blue);
#endif
}

local void write_image_file(float *a, int nx, int ny,
			    char *fn, char *filename,
			    bool compress, int format,
			    bool colormap_set, char *colormap,
			    unsigned char *red,
			    unsigned char *green,
			    unsigned char *blue)
{
    char command[1024];

    // Note that write_image will write the image left to right,
    // top to bottom, but expects "normal" ordering (bottom to top)
    // in the input array... 

//      for (int j = ny-1; j >= 0; j--) {
//  	for (int i = 0; i < nx; i++)
//  	    if (a[i+nx*j] <= 0)
//  		cerr << " ";
//  	    else
//  		cerr << "*";
//  	cerr << endl;
//      }

    if (colormap_set)
	write_image(a, nx, ny, fn, 0, colormap);	// 0 ==> no scaling
    else {
	if (format == 0)
	    write_png_image(a, nx, ny, fn, red, green, blue);	// PNG format
	else
	    write_image(a, nx, ny, fn, 0, red, green, blue);	// SUN format
    }

    if (fn) {
	if (format == 1 && compress) {

	    // Runtime compression:

	    sprintf(command, "gzip -f -q %s.sun &", filename);
	    system(command);

	} else if (format == 2) {

	    // Create a GIF file (messy!).
	    // Note that this system command (without the &) is
	    // uninterruptible (known bug!), but running in the
	    // background causes obvious errors.

	    sprintf(command,
"convert %s.sun tmp1.gif && gifsicle tmp1.gif -o tmp2.gif && mv tmp2.gif %s.gif && rm %s.sun tmp1.gif",
		    filename, filename, filename);
	    system(command);

	}
    }
}

// Some defaults:

#define L	3.0
#define NX	256
#define NY	256

main(int argc, char** argv)
{
    int count = 0, count1 = 0;
    char filename[64], im_filename[64], command[1024];
    char* fn;

    char file[64];
    strcpy(file, "-");

    bool combine = true;
    bool compress = false;
    bool gif = false;
    int format = 0;

    real l = L;
    real xleft = -L;
    real xright = L;
    real ybot = -L;
    real ytop = L;

    bool HRD = false;
    // real xleft = 5;		// suitable for HRDs only...
    // real xright = 3;
    // real ybot = -3;
    // real ytop = 3;

    bool xlim_set = false;
    bool ylim_set = false;

    int nx = NX, ny = NY;
    int n = 0, nskip = 0;

    int axis = 3;		// 1 = x, 2 = y, 3 = z

    char colormap[64];
    bool colormap_set = false;

    bool testmap = false;

    bool grid = true;
    bool mass = false;
    bool radius = false;
    int psize = 1;
    bool psize_set = false;
    int ncolor = 0;
    real index_all = -1;

    bool quiet = true;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "1cC:f:F:gGi:Hl:X:x:Y:y:mn:N:p:P:qrs:S:t";

    while ((c = pgetopt(argc, argv, param_string)) != -1) {
	switch (c) {
	    case '1':	combine = !combine;
			break;
	    case 'c':	compress = true;
			break;
	    case 'C':	strncpy(colormap, poptarg, 63);
			colormap[63] = '\0';	// just in case
	    		colormap_set = true;
			break;
	    case 'f':	strncpy(file, poptarg, 63);
			file[63] = '\0';	// just in case
			break;
	    case 'F':	format = atoi(poptarg);
			gif = (format == 2);
			break;
	    case 'g':	gif = !gif;
			format = 2 - format;
			break;
	    case 'G':	grid = !grid;
			break;
	    case 'H':	HRD = !HRD;
			break;
	    case 'l':	l = atof(poptarg);
			break;
	    case 'i':	index_all = atof(poptarg);
	    		if (index_all < 0)
			    index_all = 0;
	    		else if (index_all > 1)
			    index_all = 1;
			break;
	    case 'X':	xright = atof(poptarg);
			xlim_set = true;
			break;
	    case 'x':	xleft = atof(poptarg);
			xlim_set = true;
			break;
	    case 'Y':	ytop = atof(poptarg);
			ylim_set = true;
			break;
	    case 'y':	ybot = atof(poptarg);
			ylim_set = true;
			break;
	    case 'm':	mass = true;
			break;
	    case 'n':	n = atoi(poptarg);
	    		break;
	    case 'N':	ncolor = atoi(poptarg);
			index_all = -1;
	    		break;
	    case 'p':   psize = atoi(poptarg);
	    		// if (psize < 1) psize = 1;
	    		if (psize < 0) psize = 0;
	    		psize_set = true;
		        break;
	    case 'P':	if (poptarg[0] == 'x' || atoi(poptarg) == 1)
			    axis = 1;
	    		else if (poptarg[0] == 'y' || atoi(poptarg) == 2)
			    axis = 2;
	    		else if (poptarg[0] == 'z' || atoi(poptarg) == 3)
			    axis = 3;
			break;
	    case 'q':	quiet = !quiet;
			break;
	    case 'r':	radius = true;
			break;
	    case 's':   nx = ny = atoi(poptarg);
		        break;
	    case 'S':   nskip = atoi(poptarg);
		        break;
	    case 't':	testmap = true;
			break;
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }

    if (!psize_set) {
	if (combine && !radius)
	    psize = 0;
	else
	    psize = 1;
    }

#ifndef HAS_PNG
    if (format == 0) {
	format = 2;
	gif = true;
    }
#endif

    if (gif) compress = false;

    if (HRD) {
	if (!xlim_set) {
	    xleft = 5;
	    xright = 3;
	}
	if (!ylim_set) {
	    ybot = -3;
	    ytop = 3;
	}
    } else {
	xleft = ybot = -l;
	xright = ytop = l; 
    }

    real lx = xright - xleft;
    real ly = ytop - ybot;

    real xmin = min(xleft, xright);
    real xmax = max(xleft, xright);
    real ymin = min(ybot, ytop);
    real ymax = max(ybot, ytop);

    // PRC(xleft); PRC(xright); PRL(lx);
    // PRC(ybot); PRC(ytop); PRL(ly);

    // Note on color maps and conventions (Steve, 8/02):
    //
    // If a colormap file is specified, use it and simply use the index
    // (internal or global) as a pointer into the file.  Not recommended.
    //
    // If no colormap file is specified and we are told to use:
    //
    //     - the internal particle index (index_all < 0), use a (local)
    // 	     standard colormap, or the small-N colormap, if specified
    //	     (ncolor = 3 or 4: -N)
    //
    //     - a single index (index_all > 0: -i), use the grey map for now
    //
    //	   - the particle mass (mass: -m), use the (local) alternate colormap
    //       (note that mass supercedes index for color, if specified)
    //
    // Particle radii may be determined by
    //
    //	   - the value of psize (command-line option -p)
    //
    //	   - the mass, if mass is set (-m) and radius (-r) is not
    //
    //     - the radius (scaled by psize) if radius is set (-r)
    //
    // Adding HR diagrams raises even more options.  Color is temperature,
    // but could also be mass.  Radius is radius if available, but should
    // probably be log radius, and could use mass if radius is unknown.
    // Not fully coded yet...				(Steve, 8/02)

    unsigned char red[256], green[256], blue[256];

    if (!colormap_set) {

	// Need to clean up these options, as they can be mutually
	// incompatible.

	if (mass || HRD)
	    make_local_alternate_colormap(red, green, blue, psize);
	else if (ncolor > 0)
	    make_local_small_n_colormap(red, green, blue, ncolor, psize);
	else {
	    if (index_all < 0)
		make_local_standard_colormap(red, green, blue, psize);
//		make_standard_colormap(red, green, blue);
	    else
		make_local_greymap(red, green, blue, psize);
//		make_greymap(red, green, blue);
	}
    }

    float* a = new float[nx*ny], *zarray = new float[nx*ny];
    if (!a || !zarray) exit(1);

    if (testmap) {

	// Just draw the colormap and exit.

	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx; i++)
		*(a+j*nx+i) = (i%256)/256.;
	write_image(a, nx, ny, NULL, 0, red, green, blue);
	exit(0);
    }

    int iax, jax, kax;

    if (axis == 1) {
	iax = 1;
	jax = 2;
	kax = 0;
    } else if (axis == 2) {
	iax = 2;
	jax = 0;
	kax = 1;
    } else {
	iax = 0;
	jax = 1;
	kax = 2;
    }

    int cmin = 1000000, cmax = -1000000;
    real color_scale;

    real mmin = VERY_LARGE_NUMBER, mmax = -VERY_LARGE_NUMBER;
    real rmin = VERY_LARGE_NUMBER, rmax = -VERY_LARGE_NUMBER;
    real logfac = 1, rfac = nx/(2*l);

    // Loop over input snapshots.

    hdyn* b;
    while (b = get_hdyn(cin)) {

	float color_all = 1;
	if (index_all >= 0) color_all = 0.6667 + 0.3333*index_all; // !!

	if (nskip <= 0 || count % (nskip+1) == 0) {

	    // Make this snapshot into an image.

	    if (count1 == 0) {

		// This is the first frame to be converted into an image.
		// Determine overall scalings (and the mass range, if
		// relevant) from *this* snapshot.

		if (HRD || mass || radius || index_all < 0) {

		    // We are coloring by mass or index and sizing by mass
		    // or radius.

		    for_all_daughters(hdyn, b, bb) {
			if (index_all < 0 && bb->get_index() >= 0) {
			    cmin = min(cmin, bb->get_index());
			    cmax = max(cmax, bb->get_index());
			}
			if (mass && bb->get_mass() > 0) {
			    mmin = min(mmin, bb->get_mass());
			    mmax = max(mmax, bb->get_mass());
			}
			if (radius) {
			    rmin = min(mmin, bb->get_radius());
			    rmax = max(mmax, bb->get_radius());
			}
		    }

		    color_scale = 1.0 / (cmax-cmin);
		    if (mass && mmax > mmin) logfac = 1.0/log10(mmax/mmin);
		}

		// For a combined image, initialize here.

		if (combine) {

		    // Initialize work arrays.

		    initialize_arrays(a, zarray, nx*ny);

		    // Open the output file.  For now, we don't know how
		    // to write a PNG file to stdout.

		    if (streq(file, "-"))
			fn = NULL;
		    else {
			sprintf(filename, "%s", file);
			if (format == 0)
			    sprintf(im_filename, "%s.png", file);
			else
			    sprintf(im_filename, "%s.sun", file);
			fn = im_filename;
		    }
		}
	    }

	    if (!combine)
		initialize_arrays(a, zarray, nx*ny);

	    for_all_daughters(hdyn, b, bb) {

		real x, y, z;
		if(!HRD) {
		    vector pos = bb->get_pos();
		    x = pos[iax];
		    y = pos[jax];
		    z = pos[kax];
		}
		else {

		    // Neet to be careful where we search for star data...

		    real T_eff, L_sun, stp = -1;
		    story *st = bb->get_star_story();

		    if (st) {
			T_eff = getrq(st, "T_eff"); 
			L_sun = getrq(st, "L_eff");
			// stp = getrq(st, "Type");	// need to convert...

			if (T_eff == -VERY_LARGE_NUMBER) st = NULL;
		    }

		    if (!st) {
			st = bb->get_dyn_story();
			T_eff = getrq(st, "T"); 
			L_sun = getrq(st, "L");
			stp = getrq(st, "S");
		    }

		    x = log10(T_eff);
		    y = log10(L_sun);
		    z = 0;
		}

		// PRC(x); PRC(xmin); PRL(xmax);
		// PRC(y); PRC(ymin); PRL(ymax);

		if (x > xmin && x < xmax && y > ymin && y < ymax) {

		    // Set color (by mass or index) and radius (fixed, by mass,
		    // or by radius):

		    float color = 1;
		    real r = psize;

		    if (HRD) {

			real fac = max(0.0, min(1.0, x - 3.5));	// arbitrary
			color = 0.7 + 0.3 * fac;

		    } else if (mass && bb->get_mass() > 0) {

			// Mass scaling is logarithmic.

			real fac = log10(bb->get_mass()/mmin) * logfac;

			color = 0.6667 + 0.3333*fac;

			// Minimum size is 1 pixel (r=0), maximum radius
			// is psize.

			r *= fac;

		    } else if (index_all < 0 && bb->get_index() > 0)
			color = 0.3*(bb->get_index() - cmin) * color_scale
					+ 0.7;	// note offset !!!
		    else if (index_all >= 0)
			color = color_all;

		    if (radius || (HRD && bb->get_radius() > 0))
			r = psize * bb->get_radius() * rfac;

		    // Single pixels are too small for an animation.

		    if (!combine) r = max(r, 1.0);

		    // Coordinates in the frame:

		    x = ((x - xleft) * 1.0 * nx / lx);
		    int i = (int) x;
		    y = ((y - ybot) * 1.0 * ny / ly);
		    int j = (int) y;

		    add_point(a, nx, ny, x, y, i, j, grid, r, color, z, zarray);
		}
	    }

	    // For separate frames, create and write the image file.

	    if (!combine) {

		// Image file name counts output images, not input snaps.

		if (streq(file, "-"))
		    fn = NULL;
		else {
		    sprintf(filename, "%s.%3.3d", file, count1);
		    if (format == 0)
			sprintf(im_filename, "%s.%3.3d.png", file, count1);
		    else
			sprintf(im_filename, "%s.%3.3d.sun", file, count1);
		    fn = im_filename;
		}

		write_image_file(a, nx, ny, fn, filename, compress, format,
				 colormap_set, colormap, red, green, blue);
	    }
	}

	rmtree(b);

	count1++;
	count++;
	if (count%10 == 0 && !quiet)
	    cerr << count1 << "/" << count << " ";

	if (n > 0 && count1 >= n) break;
    }

    if (!quiet) cerr << endl;

    if (combine) {

	// Write the output file.

	write_image_file(a, nx, ny, fn, filename, compress, format,
			 colormap_set, colormap, red, green, blue);

    }
}
