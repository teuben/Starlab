
//// snap_to_image:  Construct images of a series of snapshots.
////
//// Options:
////           -1           combine all frames in a simgle image          [yes]
////           -a           produce a series of frames for animation       [no]
//   OLD       -c           compress the image file(s) using gzip          [no]
////           -c variable  specify the variable setting point color,
////                        followed by the limits to use in determining
////                        the color; variables are index, mass, radius,
////                        temperature, and limits are logarithmic except
////                        for index; three arguments are required, but
////                        setting equal limits (e.g. 0 0) does nothing,
////                        and using "." means determine the limit
////                        automatically                           [index . .] 
////                        (additional optional arguments set min and max)
////           -C colormap  specify a colormap file name                   [no]
////           -d           delete frames once the animation is made       [no]
////           -D           delay between movie frames in ms                [0]
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
//   OLD       -m           use mass to determine star color and/or size   [no]
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
////           -p psize     specify (maximum) stellar radius, in pixels (if
////                        psize < 0, use |psize} and set the  lower limit
////                        on pixel size to 0, otherwise, the limit is 1)
////                                          [0 (single image), 1 (animation)]
////           -P axis      specify projection axis                         [z]
////           -q           toggle suppression of diagnostic output
////                                                           [don't suppress]
//   OLD       -r           use stellar radius to set point size           [no]
////           -r variable  specify the variable setting point radius,
////                        followed  by the limits to use in determining
////                        the radius; variables are mass, radius,
////                        luminosity; details are as for -c above      [none]
////                        (additional optional arguments set min and max)
////           -R           animate in reverse                             [no]
////           -s nx ny     specify image size, in pixels                 [256]
////           -S nskip     specify snaps to skip between images            [0]
////           -t           test the color map [don't test]
////           -T           specify precedence scheme for points in the image
////                            (c = color, r = radius, z = depth)          [z]
////           -x           specify left (log effective temperature) edge
////                            of HRD (-H only)                          [4.5]
////           -X           specify right (log effective temperature) edge
////                            of HRD (-H only)                          [3.5]
////           -y           specify minimum (log luminosity/Lsun) limit
////                            of HRD (-H only)                         [-2.5]
////           -Y           specify maximum (log luminosity/Lsun) limit
////                            of HRD (-H only)                         [3.75]
////           -z           compress the image file(s) using gzip          [no]
////
//// In the case of GIF output, the command line responsible for creation of
//// the image is encoded in the text segment of the output file.  Note that
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
////
//// Note: the comand-line option list has changed as of June 2008.
//// The old "-c" (gzip) option has been renamed to "-z".  The old
//// "-m" and "-r" options have been removed.  The new "-c" and "-r"
//// options control point size and color using the variables index,
//// mass, radius, temperature, and luminosity, extending and
//// superceding the old "-m" and "-r" options.  The limits are
//// generally logarithmic (except in index).
//.............................................................................
//
//    version 1:  Nov 1998   Steve McMillan	 email: steve@physics.drexel.edu
//			     Drexel University, Philadelphia, PA, USA
//    version 2:  Aug 2002   Steve McMillan
//    version 2:  Jun 2008   Steve McMillan
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

void make_stellar_colormap(unsigned char* red,
			   unsigned char* green,
			   unsigned char* blue);

local void make_local_stellar_colormap(unsigned char* red,
				       unsigned char* green,
				       unsigned char* blue,
				       int psize)
{
    make_stellar_colormap(red, green, blue);

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
    const char *func = "write_image_file";

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

			// Omit merger info for all but the first frame.

			if (count == 0 || !plus || plus > sub) {
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

real get_logtemp(hdyn *bb)
{
    // Need to be careful where we search for star data.  Look in the
    // star story, then in the dyn story, then quit.

    real T_eff = -1;
    story *st = bb->get_star_story();
    if (st) T_eff = getrq(st, "T_eff"); 
    if (T_eff <= 0) {
      st = bb->get_dyn_story();
      T_eff = getrq(st, "T");	    // really pdyn
    }
    if (T_eff <= 0) {
      cerr << "No temperature data available for " << bb->format_label()
	   << endl;
      exit(1);
    }

    return log10(T_eff);
}

real get_loglum(hdyn *bb)
{
    // Need to be careful where we search for star data.  Look in the
    // star story, then in the dyn story, then quit.

    real L_eff = -1;
    story *st = bb->get_star_story();
    if (st) L_eff = getrq(st, "L_eff");
    if (L_eff <= 0) {
      st = bb->get_dyn_story();
      L_eff = getrq(st, "L");	    // really pdyn
    }
    if (L_eff <= 0) {
      cerr << "No luminosity data available for " << bb->format_label()
	   << endl;
      exit(1);
    }

    return log10(L_eff);
}

real get_color_data(hdyn *bb, int colorvar)
{
    real x;
    if (colorvar == 1)
	x = bb->get_index();
    else if (colorvar == 2)
	x = log10(bb->get_mass());
    else if (colorvar == 3) {
	x = bb->get_radius();
	if (x <= 0) {
	    cerr << "No radius data available"
		 << endl;
	    exit(1);
	}
	x = log10(x);
    } else if (colorvar == 4)
	x = get_logtemp(bb);

    return x;
}

real get_radius_data(hdyn *bb, int radiusvar)
{
    real x;
    if (radiusvar == 1)
	x = log10(bb->get_mass());
    else if (radiusvar == 2) {
	x = bb->get_radius();
	if (x <= 0) {
	    cerr << "No radius information available"
		 << endl;
	    exit(1);
	}
	x = log10(x);
    } else if (radiusvar == 3)
	x = get_loglum(bb);

    return x;
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

    bool xlim_set = false;
    bool ylim_set = false;

    int n = 0, nskip = 0;

    int origin = 0;
    int axis = 3;		// 1 = x, 2 = y, 3 = z

    char colormap[FILE_NAME_LEN];
    bool colormap_set = false;

    bool testmap = false;

    bool grid = true;
//    bool mass = false;
//    bool radius = false;
    int psize = 1;
    bool psize_set = false;
    real minpixel = 1;
    int ncolor = 0;
    real index_all = -1;

    const char *color[5] = {"none", "index", "mass", "radius", "temperature"};
    int colorvar = 0;
    bool colorvar_set = false;
    const char *radius[4] = {"none", "mass", "radius", "luminosity"};
    int radiusvar = 0;
    bool radiusvar_set = false;

    bool quiet = true;
    bool delete_frames = false;

    bool reverse = false;

    int loop = 10;

    vec xoffset, voffset, dxoffset;	// origin = 1, 3, or 4

    int delay = 0;

    // Scaling data for colors and radii.

    real cmin = VERY_LARGE_NUMBER, cmax = -VERY_LARGE_NUMBER;
    bool cmin_set = false, cmax_set = false;
    real color_scale;
    
    real rmin = VERY_LARGE_NUMBER, rmax = -VERY_LARGE_NUMBER;
    bool rmin_set = false, rmax_set = false;
    real radius_scale;

    check_help();

    extern char *poptarg, *poparr[];
    int c;
    const char *param_string = 
      "1ac:::C:dD:f:F:gGi:Hl:L:X:x:Y:y:n:N:o:O:p:P:qr:::Rs:.S:tT:z";

    char c0;
    real temp0, temp1;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision$", _SRC_)) != -1) {
	switch (c) {
	    case '1':	combine = true;
			break;
	    case 'a':	combine = false;
			break;
//	    case 'c':	compress = true;
//			break;
	    case 'c': 	c0 = poptarg[0];
			colorvar_set = true;
			if (c0 == 'i' || c0 == 'I')
			  colorvar = 1;		// index
			else if (c0 == 'm' || c0 == 'M')
			  colorvar = 2;		// mass
			else if (c0 == 'r' || c0 == 'R')
			  colorvar = 3;		// radius
			else if (c0 == 't' || c0 == 'T')
			  colorvar = 4;		// temperature
			else
			  colorvar_set = false;
			if (colorvar_set) index_all = -1;
			if (strcmp(poparr[1], ".")) {
			    cmin_set = true;
			    cmin = atof(poparr[1]);
			}
			if (strcmp(poparr[2], ".")) {
			    cmax_set = true;
			    cmax = atof(poparr[2]);
			}
			if (cmin_set && cmax_set && cmin == cmax) {
			    cmin = VERY_LARGE_NUMBER;
			    cmax = -VERY_LARGE_NUMBER;
			    cmin_set = cmax_set = false;
			}
			break;
	    case 'C':	strncpy(colormap, poptarg, 63);
			colormap[63] = '\0';	// just in case
	    		colormap_set = true;
			break;
	    case 'd':	delete_frames = !delete_frames;
	    		break;
	    case 'D':	delay = atoi(poptarg);
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
			if (HRD) {
			  if (colorvar_set && colorvar != 4)
			    cerr << "warning: -H option overrides color = \""
				 << color[colorvar] << "\"" << endl;
			  if (radiusvar_set && radiusvar != 2)
			    cerr << "warning: -H option overrides radius = \""
				 << radius[radiusvar] << "\"" << endl;
			  colorvar = 4;		// temperature
			  radiusvar = 2;	// radius
			  colorvar_set = radiusvar_set = false;
			}
			break;
	    case 'i':	index_all = atof(poptarg);
	    		if (index_all < 0)
			    index_all = 0;
	    		else if (index_all > 1)
			    index_all = 1;
			colorvar = 0;
			break;
	    case 'l':	boxw = atof(poptarg);
			break;
	    case 'L':	loop = atoi(poptarg);
	    		break;
//	    case 'm':	mass = true;
//			break;
	    case 'n':	n = atoi(poptarg);
	    		break;
	    case 'N':	ncolor = atoi(poptarg);
			colorvar = 1;
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
//	    case 'r':	radius = true;
//			break;
	    case 'r': 	c0 = poptarg[0];
			radiusvar_set = true;
			if (c0 == 'm' || c0 == 'M')
			  radiusvar = 1;	// mass
			else if (c0 == 'r' || c0 == 'R')
			  radiusvar = 2;	// radius
			else if (c0 == 'l' || c0 == 'L')
			  radiusvar = 3;	// luminosity
			else
			  radiusvar_set = false;
			if (strcmp(poparr[1], ".")) {
			    rmin_set = true;
			    rmin = atof(poparr[1]);
			}
			if (strcmp(poparr[2], ".")) {
			    rmax_set = true;
			    rmax = atof(poparr[2]);
			}
			if (rmin_set && rmax_set && rmin == rmax) {
			    rmin = VERY_LARGE_NUMBER;
			    rmax = -VERY_LARGE_NUMBER;
			    rmin_set = rmax_set = false;
			}
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
	    case 'z':	compress = true;
			break;
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			return false;
	}
    }

    PRC(origin); PRL(xoffset); PRL(dxoffset);

    // Check color and radius choices and echo settings.

    if (colorvar < 0 || colorvar > 4) colorvar = 0;
    if (radiusvar < 0 || radiusvar > 3) radiusvar = 0;

    cerr << "color variable choice is \""<< color[colorvar] << "\"" << endl;
    if (index_all >= 0) PRL(index_all);
    if (ncolor > 0) PRL(ncolor);
    cerr << "radius variable choice is \""<< radius[radiusvar] << "\"" << endl;

    if (!psize_set) {
	if (combine && radiusvar == 0) {
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
	if (!xlim_set) {	// also use these for color scaling, if chosen
	    xleft = 4.5;
	    xright = 3.5;
	}
	if (!ylim_set) {
	    ybot = -2.5;
	    ytop = 3.75;
	}
    } else {
	xleft = -boxw;
	xright = boxw;
	ybot = -(boxw*ny)/nx;
	ytop = (boxw*ny)/nx;
    }

    // Note on color maps and conventions (Steve, 8/02):
    //
    // If a colormap file is specified, use it and simply use the
    // index (internal or global) or scaled color ariable as a pointer
    // into the file.  Not recommended.
    //
    // If no colormap file is specified and we are told to use:
    //
    //     - the internal particle index (colorvar = 1), use a (local)
    // 	     standard colormap, or the small-N colormap, if specified
    // 	     (ncolor = 3 or 4: -N option)
    //
    //     - a single index (colorvar = 0; 0 < index_all < 1: -i
    //     - option), use the grey map for now
    //
    //	   - the particle mass or radius (colorvar = 1 or 2): use the
    //	     (local) alternate colormap, which runs from red to blue
    //	     to better mimic the mass
    //
    //	   - the particle luminosity (colorvar = 3): use the
    //	     (local) stellar colormap, which runs from red to voilet
    //	     to better mimic stellar colors
    //
    // Particle radii are determined by
    //
    //	   - the value of psize (radiusvar = 0; command-line option -p)
    //
    //	   - the particle mass, radius, or luminosity (radiusvar = 1,
    //	   - 2, 3), optionally scaled by psize; we actually use the
    //	   - log of the quantity, not the quantity itself

    unsigned char red[256], green[256], blue[256];

    if (!colormap_set) {

	if (colorvar == 4)
	    make_local_stellar_colormap(red, green, blue, psize);
	else if (ncolor > 0)
	    make_local_small_n_colormap(red, green, blue, ncolor, psize);
	else {
	    if (index_all < 0)
		make_local_standard_colormap(red, green, blue, psize);
	    else
		make_local_greymap(red, green, blue, psize);
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

    // Box settings (some may change with the first dataset):

    real logfac = 1, rfac = nx/(2*boxw);

    real lx = xright - xleft;
    real ly = ytop - ybot;

    real xmin = Starlab::min(xleft, xright);
    real xmax = Starlab::max(xleft, xright);
    real ymin = Starlab::min(ybot, ytop);
    real ymax = Starlab::max(ybot, ytop);

    PRC(xleft); PRC(xright); PRL(lx);
    PRC(ybot); PRC(ytop); PRL(ly);

    hdyn* b;
    string frame_id_string;

    // Loop over input snapshots.  Loop is too long -- should be split
    // into more elementary functions.

    while (b = get_hdyn()) {

	if (!combine || count == 0)
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

		// This is the first frame to be converted into an
		// image.  Initialize all dynamic quantities here.
		// Determine overall scalings and the mass and other
		// ranges (as relevant) from *this* snapshot.

		// Start by ensuring that some points will actually be
		// visible!

		if (!HRD) {

		    // Determine spatial box scale and limits.

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
			x1 = y1 = VERY_LARGE_NUMBER;
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

		// Set color scaling.

		if (ncolor > 0) {

		    cmin = 1;
		    cmax = ncolor;

		} else if (colorvar > 0) {

		    // Derive scaling information for point color.

		    if (colorvar < 4) {
			for_all_leaves(hdyn, b, bb) {
			    real x = get_color_data(bb, colorvar);
			    if (!cmin_set) cmin = Starlab::min(cmin, x);
			    if (!cmax_set) cmax = Starlab::max(cmax, x);
			}
		    } else {
			if (!cmin_set) cmin = 3.5;	// hard-wire, for now
			if (!cmax_set) cmax = 4.5;
		    }
		    color_scale = 1.0 / (cmax-cmin);
		    cerr << "color (log " << color[colorvar] << ") limits: ";
		    PRC(cmin); PRL(cmax);
		}

		// Set radius scaling.

		if (radiusvar > 0) {

		    // Derive scaling information for point radii.

		    for_all_leaves(hdyn, b, bb) {
			real x = get_radius_data(bb, radiusvar);
			if (!rmin_set) rmin = Starlab::min(rmin, x);
			if (!rmax_set) rmax = Starlab::max(rmax, x);
		    }
		    radius_scale = 1.0 / (rmax-rmin);
		    cerr << "radius (log " << radius[radiusvar] << ") limits: ";
		    PRC(rmin); PRL(rmax);
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
	    for_all_leaves(hdyn, b, bb) {

		// First determine cooordinates on the display.

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

		} else {
		    x = get_logtemp(bb);
		    y = get_loglum(bb);
		    z = 0;
		}

		// If point is visible, determine its properties.

		if (x > xmin && x < xmax && y > ymin && y < ymax) {

		    // Set color and radius.  Color is real, between 0 and 1.

		    float color = 1;	// default = top end of the colormap

		    if (index_all >= 0)
			color = color_all;
		    else if (colorvar > 0) {
			real x = get_color_data(bb, colorvar);
			color = (x-cmin)*color_scale;
		    }

		    if (color < 0) color = 0;
		    if (color > 1) color = 1;

		    real r = psize;	// default size
		    if (radiusvar > 0) {
			real x = get_radius_data(bb, radiusvar);
			r *= (x-rmin)*radius_scale;
		    }

		    // Single pixels may be too small for an animation.
		    // Allow specification of a lower limit.

		    r = Starlab::max(r, minpixel);

		    // Note that the image array runs from top to bottom,
		    // while the y coordinates run from bottom to yop...

		    // Coordinates in the frame:

		    x = ((x - xleft) * 1.0 * nx / lx);
		    int i = (int) x;
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
			    sprintf(output_file_name, "%s.%4.4d.png",
				    output_file_id, count1);
			else if (format == 1)
			    sprintf(output_file_name, "%s.%4.4d.sun",
				    output_file_id, count1);
			else if (format == 2)
			    sprintf(output_file_name, "%s.%4.4d.gif",
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

    }	// end of (excessively long) loop over input snapshots

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

//		Should be OK, but seems to fail on MacOS (Steve, 7/05).
//
//		sprintf(command,
//		    "xargs gifsicle `/bin/ls %s %s.*%s` --loopcount=%d -o %s%s",
//		    rev, output_file_id, ext1, loop, output_file_id, ext2);

//		Works on linux and MacOS.

		sprintf(command,
	    "/bin/ls %s %s.*%s | xargs gifsicle --loopcount=%d -d %d -o %s%s",
 	       rev, output_file_id, ext1, loop, delay, output_file_id, ext2);
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
