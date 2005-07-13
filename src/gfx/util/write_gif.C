
#include "stdinc.h"
#include "image_fmt.h"

#define BUFLEN 1000

local void write_gif_header(FILE *dst, int cols, int rows, int depth,
			    unsigned char *red,
			    unsigned char *green,
			    unsigned char *blue)
{
    char *pos, *buffer;
    int i;

    buffer = (char *)malloc((BUFLEN+1)*sizeof(char))+1;
    pos = buffer;

    // Header.

    *pos++ = 'G';
    *pos++ = 'I';
    *pos++ = 'F';
    *pos++ = '8';
    *pos++ = '9';
    *pos++ = 'a';
    
    // Logical screen descriptor.

    *pos++ = 0xff & cols;
    *pos++ = (0xff00 & cols)/0x100;
    *pos++ = 0xff & rows;
    *pos++ = (0xff00 & rows)/0x100;
    *pos++ = 0xf0 | (0x7&(depth-1));
    *pos++ = 0xff;
    *pos++ = 0x00;

    // Color map.

    for(i=0;i<256;i++) {
	*pos++ = 0xff & red[i];
	*pos++ = 0xff & green[i];
	*pos++ = 0xff & blue[i];
    }

    // Image descriptor.

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

#define BUFSIZE 80
#define INDENT 4

void write_gif_comment(FILE *dst, char *comment)
{
    // Add an optional text comment.  Comments can be up to 255
    // characters long, but limit ours to blocks of length 80
    // for readability.  Indent extra lines by 4 chars, and also
    // take a new line and reset the indentation if a '\n' is
    // encountered.

    int comment_size = strlen(comment);
    // PRL(comment_size);
    // PRL(comment);

    char buffer[BUFSIZE+1];
    int i = 0, indent = 0;
    while (i < comment_size) {

	// Get the next piece of the comment, starting at location i.

	int n = comment_size - i;			// chars remaining
	if (n > BUFSIZE-indent) n = BUFSIZE-indent;	// chars that will fit

	// Look for a newline in the buffered portion.

	bool newline = false;
	for (int k = 0; k < n; k++)
	    if (comment[i+k] == '\n') {
		n = k;
		newline = true;
		break;
	    }

	// Next piece will consist of n chars plus possible indentation.

	for (int k = 0; k < indent; k++) buffer[k] = ' ';
	if (indent > 0) buffer[0] = '+';
	strncpy(buffer+indent, comment+i, n);
	buffer[n+indent] = '\0';

	// PRC(n); PRL(buffer);

	// Write a new text block.

	fputc(0x21, dst);
	fputc(0xfe, dst);
	fputc(n+indent, dst);
	fputs(buffer, dst);
	fputc(0x00, dst);

	i += n;
	if (newline) {
	    i++;
	    indent = 0;
	} else
	    indent = INDENT;

	// PRC(i); PRL(indent);
    }
}

void write_gif_end_image(FILE *dst)
{
    fputc(0x00, dst);
}

void write_gif_footer(FILE *dst)
{
    fputc(0x3b, dst);	// end of gif file
}

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      unsigned char *red,
	      unsigned char *green,
	      unsigned char *blue,
	      char *comment)			// default = NULL
{
    int depth = 8;
    write_gif_header(dst, cols, rows, depth, red, green, blue);
    encode_gif(dst, pixels, depth, cols*rows);
    write_gif_end_image(dst);
    if (comment) write_gif_comment(dst, comment);
    write_gif_footer(dst);

    return OK;
}

int write_gif(FILE *dst, int cols, int rows,
	      unsigned char *pixels,
	      char *colormap_file,
	      char *comment)			// default = NULL
{
    unsigned char red[256], green[256], blue[256];
    get_colormap(red, green, blue, colormap_file);
    return write_gif(dst, cols, rows, pixels, red, green, blue, comment);
}
