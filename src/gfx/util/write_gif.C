
#include "stdinc.h"
#include "image_fmt.h"

#define BUFLEN 1000

void write_gif_header(FILE *dst, int cols, int rows, int depth,
		      char *colormap_file)
{
    char *pos, *buffer;
    int i;

    buffer = (char *)malloc((BUFLEN+1)*sizeof(char))+1;
    pos = buffer;

    *pos++ = 'G';
    *pos++ = 'I';
    *pos++ = 'F';
    *pos++ = '8';
    *pos++ = '7';
    *pos++ = 'a';
  
    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols)/0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows)/0x100;
    *pos++ = 0xf0 | (0x7&(depth-1));
    *pos++ = 0xff;
    *pos++ = 0x0;

    /* Color map. */

    unsigned char red[256], green[256], blue[256];

    get_colormap(red, green, blue, colormap_file);

    for(i=0;i<256;i++) {
	*pos++ = 0xff & red[i];
	*pos++ = 0xff & green[i];
	*pos++ = 0xff & blue[i];
    }

    *pos++ = 0x2c;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols)/0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows)/0x100;
    *pos++ = 0x7&(depth-1);

    fwrite(buffer, pos-buffer, 1, dst);
}

// Almost identical to the above (except for arguments and colormap).

void write_gif_header(FILE *dst, int cols, int rows, int depth,
		      unsigned char *red,
		      unsigned char *green,
		      unsigned char *blue)
{
    char *pos, *buffer;
    int i;

    buffer = (char *)malloc((BUFLEN+1)*sizeof(char))+1;
    pos = buffer;

    *pos++ = 'G';
    *pos++ = 'I';
    *pos++ = 'F';
    *pos++ = '8';
    *pos++ = '7';
    *pos++ = 'a';
  
    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols)/0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows)/0x100;
    *pos++ = 0xf0 | (0x7&(depth-1));
    *pos++ = 0xff;
    *pos++ = 0x0;

    /* Color map. */

    for(i=0;i<256;i++) {
	*pos++ = 0xff & red[i];
	*pos++ = 0xff & green[i];
	*pos++ = 0xff & blue[i];
    }

    *pos++ = 0x2c;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0x00;
    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols)/0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows)/0x100;
    *pos++ = 0x7&(depth-1);

    fwrite(buffer, pos-buffer, 1, dst);
}

void write_gif_footer(FILE *dst)
{
    fputc(0, dst);
    fputc(';', dst);	/* end of gif file */
}

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      char *colormap_file)
{
    int depth = 8;
    write_gif_header(dst, cols, rows, depth, colormap_file);
    encode_gif(dst, pixels, depth, cols*rows);
    write_gif_footer(dst);

    return OK;
}

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      unsigned char *red,
	      unsigned char *green,
	      unsigned char *blue)
{
    int depth = 8;
    write_gif_header(dst, cols, rows, depth, red, green, blue);
    encode_gif(dst, pixels, depth, cols*rows);
    write_gif_footer(dst);

    return OK;
}
