
#define ERROR	0
#define OK	1

void make_standard_colormap(unsigned char *red,
			    unsigned char *green,
			    unsigned char *blue);

void make_alternate_colormap(unsigned char *red,
			     unsigned char *green,
			     unsigned char *blue);

void greymap(unsigned char *red,
	     unsigned char *green,
	     unsigned char *blue);

void get_colormap(unsigned char *red,
		  unsigned char *green,
		  unsigned char *blue,
		  char *colormap_file = NULL);

void write_sun_header(int m, int n, FILE *out_file,
		      unsigned char *red,
		      unsigned char *green,
		      unsigned char *blue);

void write_sun_header(int m, int n, FILE *out_file,
		      char* colormap_file = NULL);

void write_gif_header(int m, int n, FILE *out_file,
		      unsigned char *red,
		      unsigned char *green,
		      unsigned char *blue);

void write_gif_header(int m, int n, FILE *out_file,
		      char* colormap_file = NULL);

void encode_gif(FILE *dst, unsigned char *pixels, int depth, int siz);

int write_sun(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      char *colormap_file = NULL);

int write_sun(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      unsigned char *red = NULL,
	      unsigned char *green = NULL,
	      unsigned char *blue = NULL);

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      char *colormap_file = NULL);

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      unsigned char *red = NULL,
	      unsigned char *green = NULL,
	      unsigned char *blue = NULL);

int write_png(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      char *colormap_file = NULL);

int write_png(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      unsigned char *red = NULL,
	      unsigned char *green = NULL,
	      unsigned char *blue = NULL);
