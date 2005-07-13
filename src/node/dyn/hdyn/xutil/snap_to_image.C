
//// snap_to_image:  Construct images of a series of snapshots.
////
//// Options:
////           -1           combine all frames in a simgle image          [yes]
////           -a           produce a series of frames for animation       [no]
////           -c           compress the image file(s) using gzip          [no]
////           -C colormap  specify a colormap file name                   [no]
////           -d           delete frames once the animation is made       [no]
////           -f filename  specify root name of image files   ["-" --> stdout]
////           -F format    specify image file format
////                            (0 = PNG, 1 = SUN, 2 = GIF)                 [2]
////           -g           write GIF files -- same as "-F 2"             [yes]
////           -G           toggle forcing particles to grid (nicer single
////                            frames, but jerkier movies)              [true]
////           -H           toggle Herzsprung-Russel diagram or positional
////                            plot                                 [position]
////           -i index     specify (real) color index (0-1) for all stars
////                                                       [use internal index]
////           -l scale     specify width of field of view (+/- scale)      [3]
////           -L loop      specify number of loops in animation           [10]
////           -m           use mass to determine star color and/or size   [no]
////           -n nmax      specify maximum number of images to produce   [Inf]
////           -N nbody     color using a (small-N) colormap               [no]
////           -o filename  same as -f (more standard name)
////           -O option    specify how to choose the plot center           [0]
////                            0:  as is (don't adjust)
////                            1:  initial center of mass
////                            2:  modified center of mass of each frame
////                            (x,y,z):  specify coordinates (single argument)
////                                - "(" can also be "[" or "{", but note
////                                  that all must be escaped or quoted,
////                                  since they have special meaning to
////                                  the shell
////                                - separate components by commas or
////                                  spaces (again quoted or escaped, so
////                                  the vector is a single argument)
////                                - the "{" version specifies an incremental
////                                  offest between successive frames.
////           -p psize     specify (maximum) star radius/scale, in pixels
////                        (psize < 0 ==> lower limit on pixel size = 0,
////                        otherwise, limit = 1)
////                                          [0 (single image), 1 (animation)]
////           -P axis      specify projection axis                         [z]
////           -q           toggle suppression of diagnostic output
////                                                           [don't suppress]
////           -r           use stellar radius to set point size           [no]
////           -R           animate in reverse                             [no]
////           -s nx ny     specify image size, in pixels                 [256]
////           -S nskip     specify snaps to skip between images            [0]
////           -t           test the color map [don't test]
////           -T           specify precedence scheme for points in the image
////                            (c = color, r = radius, z = depth)          [z]
////           -x           specify right (log effective temparature) edge
////                            of HRD (-H only)                            [3]
////           -X           specify left (log effective temparature) edge
////                            of HRD (-H only)                            [5]
////           -y           specify minimum (log luminosity/Lsun) limit
////                            of HRD (-H only)                           [-3]
////           -Y           specify maximum (log luminosity/Lsun) limit
////                            of HRD (-H only)                            [3]
////
//// In the case of GIF output, the command line responsible for creation of
//// the image is encoded in the text segmant of the output file.  Note that
//// the default output format has been changed to GIF now that the LZW patent
//// issues seem to have gone away.
////
//// Notes: 1. If animations are specified, an MNG or animated GIF file will
////           be created.  The individual frames will be retained unless
////           the "-d" option is set.
////        2. If PNG output is requested and the PNG libraries are unavailable,
////           then GIF output is produced instead.
////        3. The "T" option determines which attribute is favored in
////           displaying a particle in the image.  The default is depth,
////           but we can also raise particles based on color or radius.
//.............................................................................
//
//    version 1:  Nov 1998   Steve McMillan	 email: steve@zonker.drexel.edu
//			     Drexel University, Philadelphia, PA, USA
//    version 2:  Aug 2002   Steve McMillan
//.............................................................................



#include "hdyn.h"
#include "star/single_star.h"
#include "../../../../gfx/util/image_fmt.h"	// should clean this up!

// Our three-color convention now requires us to define our own standards,
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




local void initialize_arrays(unsigned char *a,
			     float *zarray, float *rarray,
			     int n)
{
    for (int i = 0; i < n; i++) {
	a[i] = 0;
	zarray[i] = -1.e10;
	rarray[i] = -1.e10;
    }
}

int precedence = 0;

local void add_point(unsigned char *a, int nx, int ny,
		     real x, real y, int i, int j,
		     bool grid,
		     real r, float color, real z,
		     float *zarray, float *rarray)
{
    // Do the actual work of adding a particle to the image grid (a).

    // The true position of the particle in pixel coordinates
    // is (x,y).  The central pixel is (i,j).  We expect i <= x,
    // j <= y.  The desired radius is r pixels.  On input, color is
    // a number between 0 and 1 -- it will be rescaled to fit the
    // color map when added to the array.

    // With our new scheme, the color used should actually be at
    // the bright end of the colormap, in the range 0.667-1.0.  The
    // calling function need not know that, so expect color in the
    // range 0 to 1 and rescale it here before continuing.

    color = 0.68 + 0.31*color;

    real xref = x;
    real yref = y;
    if (grid) {
	xref = i+0.5;
	yref = j+0.5;
    }

    int ir = (int)(r+0.0001)+1;
    real r2a = pow(r-0.1,2), r2b = pow(r+0.2,2), r2c = pow(r+0.5,2);

    // PRC(x); PRC(i); PRC(y); PRL(j);

    for (int ii = Starlab::max(-ir, -i);
	 ii <= Starlab::min(ir, nx-i-1); ii++)
	for (int jj = Starlab::max(-ir, -j);
	     jj <= Starlab::min(ir, ny-j-1); jj++) {

	    // Get distance from pixel center to reference point.

	    real xi = i+ii+0.5;
	    real yi = j+jj+0.5;

	    real dx = abs(xi-xref);
	    real dy = abs(yi-yref);

	    real dx2 = pow(xi-xref,2) + pow(yi-yref,2);

	    // Soften the edges of the stars.  (Details fine-tuned by eye.)
	    // Use 3 color levels, in addition to black.

	    real c = color;			// convenient to use real here

	    if ((ii == 0 && jj == 0) || dx2 <= r2a)
		; 				// "color 1" -- do nothing
	    else if (dx2 < r2b)
		c -= 0.33;			// "color 2" (index - 85)
	    else if (dx2 < r2c)
		c -= 0.67;			// "color 3" (index - 170)
	    else
		c = 0;				// black

	    if (c > 0.0039) {			// 0.0039 = 1/256

#if 0
		cerr << "adding " << c << " to ("
		     << i+ii << "," << j+jj << ")" << endl;
#endif

		if (precedence == 0) {

		    // Placement is determined by z.

		    real zcur = *(zarray+(j+jj)*nx+(i+ii));
		    if (z > zcur) {
			*(zarray+(j+jj)*nx+(i+ii)) = z;
			*(a+(j+jj)*nx+(i+ii)) = (unsigned char)(255.999*c);
		    }
		
		} else if (precedence == 1) {

		    // Placement is determined by color.

		    unsigned char cc = (unsigned char)(255.999*c);

		    if (cc >= *(a+(j+jj)*nx+(i+ii)))
			*(a+(j+jj)*nx+(i+ii)) = cc;

		} else if (precedence == 2) {

		    // Placement is determined by radius.

		    real rcur = *(rarray+(j+jj)*nx+(i+ii));
		    if (r > rcur) {
			*(rarray+(j+jj)*nx+(i+ii)) = r;
			*(a+(j+jj)*nx+(i+ii)) = (unsigned char)(255.999*c);
		    }
		}
	    }

	}
}



local void write_image_file(unsigned char *a, int nx, int ny,
			    char *output_file_name,
			    bool compress, int format,
			    bool colormap_set, char *colormap,
			    unsigned char *red,
			    unsigned char *green,
			    unsigned char *blue,
			    const char *comment = NULL)
{
    char *func = "write_image_file";

    // Note that write_xxx will write the image left to right,
    // top to bottom.

//      for (int j = ny-1; j >= 0; j--) {
//  	for (int i = 0; i < nx; i++)
//  	    if (a[i+nx*j] <= 0)
//  		cerr << " ";
//  	    else
//  		cerr << "*";
//  	cerr << endl;
//      }

    FILE *dst = stdout;			// NULL file name ==> stdout

    if (output_file_name) {
	dst = fopen(output_file_name, "w");
	if (!dst) {
	    cerr << func << " can't open file " << output_file_name
		 << "; using stdout" << endl;
	    dst = stdout;
	    output_file_name = NULL;
	}
    }

    if (colormap_set) {

	if (format == 0)
	    write_png(dst, nx, ny, a, colormap, comment);
	else if (format == 1)
	    write_sun(dst, nx, ny, a, colormap);
	else if (format == 2)
	    write_gif(dst, nx, ny, a, colormap, comment);

    } else {

	if (format == 0)
	    write_png(dst, nx, ny, a, red, green, blue, comment);
	else if (format == 1)
	    write_sun(dst, nx, ny, a, red, green, blue);
	else if (format == 2)
	    write_gif(dst, nx, ny, a, red, green, blue, comment);

    }

    if (output_file_name) {

	fclose(dst);

	if (compress) {

	    // Runtime compression of Sun rasterfile (output_file_name
	    // includes the frame number):

	    char command[256];
	    sprintf(command, "gzip -f -q %s&", output_file_name);
	    system(command);
	}
    }
}



// Some defaults:

#define L	3.0
#define NX	256
#define NY	256

#define FILE_NAME_LEN	128

#include <ctype.h>

#define NARGS 8

local void get_offset(char *str, vec& offset)
{
    int len = strlen(str)+1;
    char *copy = new char[len];
    strncpy(copy, str, len);

    // Convert the string into suitable form.

    for (int i = 0; i < len-1; i++) {
	int c = copy[i];
	if (isdigit(c) || c == 'e' || c == 'E'
	    || c == '-' || c == '+' || c == '.')
	    ;
	else
	    copy[i] = ' ';
    }
    sscanf(copy, "%lf %lf %lf", &offset[0], &offset[1], &offset[2]);
}

#include <string>
#include <sstream>

local string identify_run(int argc, char *argv[])
{
    ostringstream s;

    // Identify the Starlab version.

    s << "Image created by Starlab version " << VERSION;

    // Attempt to identify the user and host.

    s << ", run by ";
    if (getenv("LOGNAME"))
	s << getenv("LOGNAME");
    else
	s << "(unknown)";

    s << " on host ";
    if (getenv("HOST"))
	s << getenv("HOST");
    else if (getenv("HOSTNAME"))
	s << getenv("HOSTNAME");
    else
	s << "(unknown)";

    s << endl;

    if (getenv("HOSTTYPE"))
	s << "    (host type " << getenv("HOSTTYPE") <<")";

    // Also log the time and date.

    const time_t tt = time(NULL);
    s << " on  " << ctime(&tt);

    s << "Command: " << argv[0];
    for (int i = 1; i < argc; i++)
	s << " " << argv[i];

    return s.str();		// no newline at end
}

local string identify_frame(hdyn *b, int count)
{
    ostringstream s;

    // Identify the frame.

    s << "Frame " << count << endl;

    if (b == NULL)

	s << "No snapshot information";

    else {

	s << "History:" << endl;

	// Add history from the root log story.

	story *log = b->get_log_story();
	for (story * d = log->get_first_daughter_node(); d != NULL;
	     d = d->get_next_story_node())
	    if (!d->get_chapter_flag()) {
		char *t = d->get_text();
		if (t) {
		    char *sub = strstr(t, "Starlab");
		    if (sub) {
			char *plus = strstr(t, "+");
			if (!plus || plus > sub) {	// omit merger info
			    if (strstr(sub, "(user") && strstr(sub, ":"))
				s << t << endl;
			}
		    }
		}
	    }
	s << "time = " << b->get_system_time() << "  N = " << b->n_leaves();
    }

    return s.str();		// no newline at end
}

main(int argc, char** argv)
{
    char output_file_id[FILE_NAME_LEN];
    strcpy(output_file_id, "-");

    int count = 0, count1 = 0;

    bool combine = true;
    bool compress = false;
    int format = 2;

    int nx = NX, ny = NY;

    real boxw = L;
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

    int n = 0, nskip = 0;

    int origin = 0;
    int axis = 3;		// 1 = x, 2 = y, 3 = z

    char colormap[FILE_NAME_LEN];
    bool colormap_set = false;

    bool testmap = false;

    bool grid = true;
    bool mass = false;
    bool radius = false;
    int psize = 1;
    bool psize_set = false;
    real minpixel = 1;
    int ncolor = 0;
    real index_all = -1;

    bool quiet = true;
    bool delete_frames = false;

    bool reverse = false;

    int loop = 10;

    vec xoffset, voffset, dxoffset;	// origin = 1, 3, or 4

    check_help();

    extern char *poptarg, *poparr[];
    int c;
    char* param_string = "1acC:df:F:gGi:Hl:L:X:x:Y:y:mn:N:o:O:p:P:qrRs:.S:tT:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case '1':	combine = true;
			break;
	    case 'a':	combine = false;
			break;
	    case 'c':	compress = true;
			break;
	    case 'C':	strncpy(colormap, poptarg, 63);
			colormap[63] = '\0';	// just in case
	    		colormap_set = true;
			break;
	    case 'd':	delete_frames = !delete_frames;
	    		break;
	    case 'f':
	    case 'o':	strncpy(output_file_id, poptarg, 63);
			output_file_id[63] = '\0';	// just in case
			break;
	    case 'F':	format = atoi(poptarg);
			break;
	    case 'g':	format = 2;
			break;
	    case 'G':	grid = !grid;
			break;
	    case 'H':	HRD = !HRD;
			break;
	    case 'i':	index_all = atof(poptarg);
	    		if (index_all < 0)
			    index_all = 0;
	    		else if (index_all > 1)
			    index_all = 1;
			break;
	    case 'l':	boxw = atof(poptarg);
			break;
	    case 'L':	loop = atoi(poptarg);
	    		break;
	    case 'm':	mass = true;
			break;
	    case 'n':	n = atoi(poptarg);
	    		break;
	    case 'N':	ncolor = atoi(poptarg);
			index_all = -1;
	    		break;
	    case 'O':	if (poptarg[0] == '[' || poptarg[0] == '('
			    || poptarg[0] == '{') {
			    get_offset(poptarg, dxoffset);
			    if (poptarg[0] != '{') {
				xoffset = dxoffset;
				origin = 3;
			    } else
				origin = 4;
			} else
			    origin = atoi(poptarg);
	    		break;
	    case 'p':   psize = atoi(poptarg);
	    		if (psize < 0) {
			    psize = -psize;
			    minpixel = 0;
			}
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
	    case 'R':	reverse = true;
			break;
	    case 's':   nx = ny = atoi(poparr[0]);
	    		if (poparr[1][0] != '-') ny = atoi(poparr[1]);
		        break;
	    case 'S':   nskip = atoi(poptarg);
		        break;
	    case 't':	testmap = true;
			break;
	    case 'T':	if (poptarg[0] == 'z')
			    precedence = 0;
			else if (poptarg[0] == 'c')
			    precedence = 1;
			else if (poptarg[0] == 'r')
			    precedence = 2;
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
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }


    PRC(origin); PRL(xoffset); PRL(dxoffset);

    if (!psize_set) {
	if (combine && !radius) {
	    psize = 0;
	    minpixel = 0;
	} else {
	    psize = 1;
	    minpixel = 1;
	}
    }

    string run_id_string = identify_run(argc, argv);

#ifndef HAVE_LIBPNG
    if (format == 0) {
	cerr << "No PNG support; switching to GIF" << endl;
	format = 2;
    }
#endif

    if (format != 1) compress = false;

    // At this point, output_file_id is the "root" name for the image file(s).
    // It is "-" for stdout.  Otherwise, the actual image file will be of the
    // form "output_file_id.png" or "output_file_id.nnn.png".  The string
    // output_file_name passed to write_image_file will be NULL for output
    // to stdout, or else it will contain the actual name of the file currently
    // being written, including frame number (if any) and extension.  Declare
    // it here so the name doesn't vanish!

    char *output_file_name = NULL;

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
	xleft = -boxw;
	xright = boxw;
	ybot = -(boxw*ny)/nx;
	ytop = (boxw*ny)/nx;
    }

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
    //	   - the particle mass (mass: -m), use the (local) alternate colormap,
    //	     which runs from red to blue to better mimic the mass (note that
    //	     mass supercedes index for color, if specified)
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

    unsigned char *a = new unsigned char[nx*ny];
    float *zarray = new float[nx*ny];
    float *rarray = new float[nx*ny];
    if (!a || !zarray || !rarray) exit(1);

    if (testmap) {

	// Just draw the colormap and exit.

	for (int j = 0; j < ny; j++)
	    for (int i = 0; i < nx; i++)
		*(a+j*nx+i) = (i%256);
	string image_comment = run_id_string+"\nTest image";
	write_image_file(a, nx, ny, NULL,		// NULL ==> stdout
			 0, format, false, NULL,
			 red, green, blue,
			 image_comment.c_str());
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

    // Box settings (some may change with the first dataset):

    real logfac = 1, rfac = nx/(2*boxw);

    real lx = xright - xleft;
    real ly = ytop - ybot;

    real xmin = Starlab::min(xleft, xright);
    real xmax = Starlab::max(xleft, xright);
    real ymin = Starlab::min(ybot, ytop);
    real ymax = Starlab::max(ybot, ytop);

    // PRC(xleft); PRC(xright); PRL(lx);
    // PRC(ybot); PRC(ytop); PRL(ly);

    // Loop over input snapshots.

    hdyn* b;
    string frame_id_string;

    while (b = get_hdyn()) {

	frame_id_string = identify_frame(b, count);

	if (origin == 2) {

	    // Move the standard center to (0,0,0).

	    vec cmx, cmv;
	    get_std_center(b, cmx, cmv);

	    b->inc_pos(-cmx);
	    b->inc_vel(-cmx);
	}
	
	float color_all = 1;
	if (index_all >= 0) color_all = index_all;

	if (nskip <= 0 || count % (nskip+1) == 0) {

	    // Make this snapshot into an image.

	    if (origin == 4) {
	        xoffset += dxoffset;
		PRL(xoffset);
	    }

	    if (count1 == 0) {

		// This is the first frame to be converted into an image.
		// Determine overall scalings (and the mass range, if
		// relevant) from *this* snapshot.

		// Start by checking that some points will actually be visible!

		if (!HRD) {
		    int total, count;
		    real boxw_save = boxw;
		    boxw /= 2;
		    do {
			boxw *= 2;
			real boxh = (boxw*ny)/nx;
			total = count = 0;
			for_all_daughters(hdyn, b, bb) {
			    total++;
			    vec pos = b->get_pos() + bb->get_pos();
			    real x = pos[iax];
			    real y = pos[jax];
			    if (abs(x) <= boxw && abs(y) < boxh) count++;
			}
		    }  while (total > 10 && count < (3*total)/4);

		    if (boxw > 1.1*boxw_save) {

			cerr << "snap_to_image: box size increased to "
			     << boxw << endl;

			// Repeats lines above.  Note that we are assuming
			// xleft < xright (and similarly y) in this case.

			xleft = -boxw;
			xright = boxw;
			ybot = -(boxw*ny)/nx;
			ytop = (boxw*ny)/nx;
			rfac = nx/(2*boxw);

			lx = xright - xleft;
			ly = ytop - ybot;

			xmin = Starlab::min(xleft, xright);
			xmax = Starlab::max(xleft, xright);
			ymin = Starlab::min(ybot, ytop);
			ymax = Starlab::max(ybot, ytop);
		    }

		    if (origin == 1) compute_com(b, xoffset, voffset);

		    if (total == count) {

			// Try to center the first frame.

			real x1, x2, y1, y2;
			x1 = y2 = VERY_LARGE_NUMBER;
			x2 = y2 = -VERY_LARGE_NUMBER;

			for_all_daughters(hdyn, b, bb) {

			    vec pos = b->get_pos() + bb->get_pos();
			    if (origin == 1 || origin >= 3) pos -= xoffset;
			    real x = pos[iax];
			    real y = pos[jax];

			    if (x > xmin && x < xmax && y > ymin && y < ymax) {
			
				x1 = Starlab::min(x1, x);
				x2 = Starlab::max(x2, x);
				y1 = Starlab::min(y1, y);
				y2 = Starlab::max(y2, y);
			    }
			}

			// Expand borders to the equal the largest gap.

			if (x1-xmin < xmax-x2)
			    xmin = x1 - (xmax-x2);
			else
			    xmax = x2 + (x1-xmin);

			if (y1-ymin < ymax-y2)
			    ymin = y1 - (ymax-y2);
			else
			    ymax = y2 + (y1-ymin);

			xleft = xmin;
			xright = xmax;
			ybot = ymin;
			ytop = ymax;

			lx = xright - xleft;
			ly = ytop - ybot;
			rfac = nx/(2*boxw);
		    }
		}

		if (HRD || mass || radius || index_all < 0) {

		    // We are coloring by mass or index and sizing by mass
		    // or radius.

		    // Special case:

		    if (ncolor > 0) {
			cmin = 1;
			cmax = ncolor;
		    }

		    for_all_leaves(hdyn, b, bb) {
			if (ncolor == 0 && index_all < 0
			    && bb->get_index() >= 0) {
			    cmin = Starlab::min(cmin, bb->get_index());
			    cmax = Starlab::max(cmax, bb->get_index());
			}
			if (mass && bb->get_mass() > 0) {
			    mmin = Starlab::min(mmin, bb->get_mass());
			    mmax = Starlab::max(mmax, bb->get_mass());
			}
			if (radius) {
			    rmin = Starlab::min(mmin, bb->get_radius());
			    rmax = Starlab::max(mmax, bb->get_radius());
			}
		    }

		    color_scale = 1.0 / (cmax-cmin);
		    if (mass && mmax > mmin) logfac = 1.0/log10(mmax/mmin);
		}

		// For a combined image, initialize here.

		if (combine) {

		    // Initialize work arrays.

		    initialize_arrays(a, zarray, rarray, nx*ny);

		    // Define the output file name.

		    if (strcmp(output_file_id, "-")) {

			if (output_file_name) delete [] output_file_name;
			output_file_name = new char[FILE_NAME_LEN+4];

			if (output_file_name) {
			    if (format == 0)
				sprintf(output_file_name, "%s.png",
					output_file_id);
			    else if (format == 1)
				sprintf(output_file_name, "%s.sun",
					output_file_id);
			    else if (format == 2)
				sprintf(output_file_name, "%s.gif",
					output_file_id);
			} else
			    cerr << "Can't create filename string! "
				 << "Output will go to stdout." << endl;
		    }

		}
	    }

	    if (!combine)
		initialize_arrays(a, zarray, rarray, nx*ny);

	    hdyn *root = b->get_root();
	    for_all_daughters(hdyn, b, bb) {

		real x, y, z;
		if(!HRD) {
		    vec pos = bb->get_pos();

		    hdyn *p = bb->get_parent();
		    while (p) {
			pos += p->get_pos();
			p = p->get_parent();
		    }

		    if (origin == 1 || origin >= 3) pos -= xoffset;

		    x = pos[iax];
		    y = pos[jax];
		    z = pos[kax];

		}
		else {

		    // Need to be careful where we search for star data...

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
			T_eff = getrq(st, "T");		// these are really
			L_sun = getrq(st, "L");		// pdyn data...
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
		    // or by radius).  Color is a float between 0 and 1.

		    float color = 1;
		    real r = psize;

		    if (HRD)

			color = Starlab::max(0.0,		  // ~arbitrary
					     Starlab::min(1.0,
							  x - 3.5));

		    else if (mass && bb->get_mass() > 0) {

			// Mass scaling is logarithmic.

			color = log10(bb->get_mass()/mmin) * logfac;

			// Minimum size is 1 pixel (r = 0), maximum radius
			// is psize.

			r *= color;

			// PRC(bb->get_mass()); PRC(mmin); PRC(color); PRL(r);

		    } else if (index_all < 0 && bb->get_index() > 0)

			color = (bb->get_index() - cmin)
					* color_scale;	    // note offset !!!

		    else if (index_all >= 0)

			color = color_all;

		    if (radius || (HRD && bb->get_radius() > 0))
			r = psize * bb->get_radius() * rfac;

		    // Single pixels may be too small for an animation.
		    // Allow specification of a lower limit.

		    r = Starlab::max(r, minpixel);

		    // Note that the image array runs from top to bottom,
		    // while the y coordinates run from bottom to yop...

		    // Coordinates in the frame:

		    x = ((x - xleft) * 1.0 * nx / lx);
		    int i = (int) x;
		    // y = ((y - ybot) * 1.0 * ny / ly);
		    y = ((ytop - y) * 1.0 * ny / ly);
		    int j = (int) y;

		    add_point(a, nx, ny, x, y, i, j,
			      grid, r, color, z, zarray, rarray);
		}
	    }

	    // For separate frames, create and write the image file.
	    
	    if (!combine) {

		// Image file name counts output images, not input snaps.

		if (strcmp(output_file_id, "-")) {

		    if (output_file_name) delete [] output_file_name;
		    output_file_name = new char[FILE_NAME_LEN+8];

		    if (output_file_name) {
			if (format == 0)
			    sprintf(output_file_name, "%s.%3.3d.png",
				    output_file_id, count1);
			else if (format == 1)
			    sprintf(output_file_name, "%s.%3.3d.sun",
				    output_file_id, count1);
			else if (format == 2)
			    sprintf(output_file_name, "%s.%3.3d.gif",
				    output_file_id, count1);
		    }
		}

		string image_comment = run_id_string + "\n" + frame_id_string;
		write_image_file(a, nx, ny, output_file_name,
				 compress, format,
				 colormap_set, colormap, red, green, blue,
				 image_comment.c_str());
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

	char temp[10];
	sprintf(temp, "%d", count);
	string image_comment = run_id_string
				+ "\nCombined frame of " + temp + " snapshots";
	if (frame_id_string.size() > 0)
	    image_comment += "\n" + frame_id_string;

	write_image_file(a, nx, ny, output_file_name,
			 compress, format,
			 colormap_set, colormap, red, green, blue,
			 image_comment.c_str());

    } else {

	if (format == 0 || format == 2) {

	    // Create the movie.

	    char ext1[8], ext2[8], command[1024];
	    char rev[3] = "  ";
	    if (reverse) strcpy(rev, "-r");

	    if (format == 0) {
		strcpy(ext1, ".png");
		strcpy(ext2, ".mng");
		sprintf(command, "png2mng.pl %s -i %s -s %d %d > %s%s",
			rev, output_file_id, nx, ny, output_file_id, ext2);
	    } else {
		strcpy(ext1, ".gif");
		strcpy(ext2, ".gif");
		sprintf(command,
		    "xargs gifsicle `/bin/ls %s %s.*%s` --loopcount=%d -o %s%s",
		    rev, output_file_id, ext1, loop, output_file_id, ext2);
	    }

	    cerr << endl
		 << "Creating animation in " << output_file_id << ext2 << endl;
	    int status = system(command);

	    // Delete individual frames if desired.

	    if (status == 0 && delete_frames) {
		sprintf(command, "/bin/rm -f %s.*%s", output_file_id, ext1);
		cerr << "Deleting individual frames" << endl;
		system(command);
	    }
	}
    }
}
