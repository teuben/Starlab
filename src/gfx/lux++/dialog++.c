
/* dialog.c: Lux-based routines to permit easy construction of simple
 *	     dialog boxes.
 *
 *	     Steve McMillan (Drexel U.), Jan-May 1996.
 */

#include "c_lux.h"
#include "local.h"

local int mymin(int a, int b)
{
    return (a < b ? a : b);
}

local int mymax(int a, int b)
{
    return (a > b ? a : b);
}

local int count_left_digits(double d)
{
    /* How many digits to the left of the decimal point in d? */

    int n = 1;

    if (d < 0) {
	d = -d;
	n = 2;
    }

    while (d > 10.0) {
	n++;
	d /= 10.0;
    }

    return n;
}

local int format_double_string(double d, int w, char* s)
{
    /* Attempt to format the value d into string s, of maximum width w. */

    char* temp;
    char fmt[40];
    int size;

    size = w;

    /* Try f format first. */

    sprintf(fmt, "%%.%df", mymax(2, w - 1 - count_left_digits(d)));
    sprintf(s, fmt, d);

    /* Strip trailing zeros (and the decimal point, if possible). */

    if (strchr(s, '.')) {

	temp = s + strlen(s) - 1;
	while (*temp == '0') {
	    *temp = '\0';
	    temp--;
	}
	if (*temp == '.') *temp = '\0';
    }

    if (strlen(s) > w) {

	/* Won't fit -- resort to e format. */

	sprintf(fmt, "%%.%de", mymax(2, w - 7));
	sprintf(s, fmt, d);

	if (strlen(s) > w) size = strlen(s);

    }

    return size;
}

local int convert_value(dialog_entry* d, char* s)
{
    int i, size;
    double dd;
    char* temp;

    switch (d->type) {

	case 0:	*s = *(d->value);
		size = 1;
		break;

	case 1:	i = sprintf(s, "%d", *((int*)(d->value)));
		size = 3;
		if (i > 3 && i <= 8) size = i;
		if (i > 8) size = 8;
		break;

	case 2:
	case 3: if (d->type == 2)
	    	    dd = (double)*((float*)(d->value));
		else
		    dd = *((double*)(d->value));

		size = format_double_string(dd, d->width, s);

		break;

	case 4: size = strlen(d->value);

		/* Strip leading blanks */
		temp = d->value;
	 	while(*temp <= ' ' && *temp > '\0') temp++;

		sprintf(s, "%s", temp);
		break;
    }

    if (size < d->width) size = d->width;
    return size;
}

local void set_dialog_item(Window* dialog_win, int i, dialog_entry* d,
			   int x1, int x2, int y)
{
    /* Translate d into a string representing the value. */

    int size, j;
    char temp_buffer[80];

    size = convert_value(d, temp_buffer);

    /* Initialize the appropriate box (twice, just in case...). */

    for (j = 0; j < 2; j++) {

	if (d->type) {
	    lux_set_item(*dialog_win, i, TEXT_WINDOW, NO_TYPE,
			 x1, y, strlen(d->name), d->name);
	    lux_set_item(*dialog_win, i, INPUT_WINDOW, NO_TYPE,
			 x2, y, size, temp_buffer);
	    
	} else {
	    lux_set_item(*dialog_win, i, TEXT_WINDOW, CHECK_BUTTON,
			 x1, y, strlen(d->name), d->name);
	    lux_set_item(*dialog_win, i, BUTTON_WINDOW, CHECK_BUTTON,
			 x2, y, size, temp_buffer);
	}

    }
}

void initialize_dialog(Window* dialog_win,
		       int xorigin, int yorigin, char* title,
		       int* xsize, int* ysize,
		       dialog_entry* list, int n, int split)
{
    int l1, h, lchar, l2, xgap, xmid, dy, x1, x2, ymax, y, uchar;

    /* Get maximum name length. */

    int i_max = -1, max_length = -1;
    int i_maxb = -1, ll, max_box = -1;

    char temp[40];

    int i, l;
    for (i = 0; i < n; i++) {

	if ( (l = strlen(list[i].name)) > max_length) {
	    max_length = l;
	    i_max = i;
	}

	if ((ll = convert_value(list+i, temp)) > max_box) {
	    max_box = ll;
	    i_maxb = i;
	}
    }

    if (i_max < 0 || max_length < 0 || i_maxb < 0 || max_box < 0) {
	fprintf(stderr, "Error in initialize_dialog.\n");
	exit(1);
    }

    lux_string_dimensions(list[i_max].name, &l1, &h);

    lux_string_dimensions("0", &lchar, &h);
    l2 = lchar;
    if (max_box > 0) l2 *= max_box;
    l2 += lchar;			/* Take box itself into account */

    xgap = 10;
    xmid = 15;
    dy = 15;

    if (!split) {
	*xsize = 3*xgap + l1 + l2;
	*ysize = (n+1) * (dy + h) + dy;
    } else {
	*xsize = 5*xgap + xmid + 2*(l1 + l2);
	*ysize = ((split < n-split ? n-split : split)+1) * (dy + h) + dy;
    }
    if (*xsize < 3*xgap + 12*lchar) *xsize = 3*xgap + 12*lchar; /* OK/CANCEL */

    /* Make room for the box title */

    if (title)
	i = strlen(title);
    else
	i = 10;

    i += 8;	/* Take size of buttons. etc., into account */

    lux_string_dimensions("X", &uchar, &h);
    i *= uchar*1.25;	/* Trial and error! */
    if (*xsize < i) *xsize = i;

    /* Initialize dialog box material. */

    *dialog_win = lux_open_dialog(xorigin, yorigin, *xsize, *ysize);

    if (title)
	lux_set_window_name(*dialog_win, title);
    else
	lux_set_window_name(*dialog_win, "Parameters");

    x1 = xgap;
    x2 = x1 + l1 + xgap;

    ymax = *ysize - dy/2;
    y = ymax + dy/2;

    for (i = 0; i < n; i++) {
	set_dialog_item(dialog_win, i, list+i, x1,  x2, (y -= dy+h));
	if (split > 0 && i == split-1) {
	    y = ymax + dy/2;
	    x1 = 3*xgap + xmid + l1 + l2;
	    x2 = x1 + l1 + xgap;
	}
    }

    /* OK/cancel boxes: */

    xgap = (*xsize - 12*lchar)/3;
    lux_set_item(*dialog_win, n, BUTTON_WINDOW, OK_KEEP_BUTTON,
		 xgap, dy, 6, "  OK  ");
    lux_set_item(*dialog_win, n+1, BUTTON_WINDOW, CANCEL_BUTTON,
		 2*xgap + 6*lchar, dy, 6, "CANCEL");
}

/*---------------------------------------------------------------------- */

void update_from_dialog(Window dialog_win, dialog_entry* list, int n)
{
    /* Read all items from the dialog box. */

    int i;

    for (i = 0; i < n; i++) {

	if (list[i].type == 0)
	    lux_get_itemvalue(dialog_win, i, BUTTON_WINDOW, CHECK_BUTTON,
			      list[i].value);
	else {

	    if (list[i].type < 4) {

		char temp[40];
		lux_get_itemvalue(dialog_win, i, INPUT_WINDOW, NO_TYPE, temp);

		if (list[i].type == 1)
		    sscanf(temp, "%d", (int*)list[i].value);
		else if (list[i].type == 2)
		    sscanf(temp, "%f", (float*)list[i].value);
		else
		    sscanf(temp, "%lf", (double*)list[i].value);

	    } else

		lux_get_itemvalue(dialog_win, i, INPUT_WINDOW, NO_TYPE,
				  list[i].value);

	}

    }
}

void reset_dialog_entries(Window dialog_win, dialog_entry* list, int n)
{
    /* Reset all items in the dialog box. */

    int i;

    for (i = 0; i < n; i++) {
	if (list[i].type == 0)

	    lux_update_itemvalue(dialog_win, i, BUTTON_WINDOW, CHECK_BUTTON,
				 list[i].value);

	else {

	    char temp[40];
	    convert_value(list+i, temp);

	    lux_update_itemvalue(dialog_win, i, INPUT_WINDOW, NO_TYPE, temp);
	}
    }
}

/*---------------------------------------------------------------------- */

/* For "easy" dialog options: */

static int n_dialog = 0;
static int split = 0;
static dialog_entry * d = NULL;

void print_dialog_list(char* s)
{
    int i;

    fprintf(stderr, "%s\n", s);
    for (i = 0; i < n_dialog; i++) {
	fprintf(stderr, "Item %d: name = \"%s\", type = %d, ",
		i, d[i].name, d[i].type);
	if (d[i].type == 0)
	    fprintf(stderr, "value = %d\n", *d[i].value);
	else if (d[i].type == 1)
	    fprintf(stderr, "value = %d\n", *(int*)d[i].value);
	else if (d[i].type == 2)
	    fprintf(stderr, "value = %g\n", *(float*)d[i].value);
	else if (d[i].type ==3)
	    fprintf(stderr, "value = %g\n", *(double*)d[i].value);
	else
	    fprintf(stderr, "value = %s\n", d[i].value);
    }
}

local void extend_dialog_list()
{
    /* Untested realloc!! */

    n_dialog++;
    d = realloc(d, n_dialog*sizeof(dialog_entry));
}

/* Note: These functions could be merged into one in C++! */

void create_button_dialog_entry(char* name, char* value)
{
    extend_dialog_list();

    /*
     *	Define buttons:
     *
     *		label,  type,  address of variable.
     *
     *  Types:  0 = char
     *		1 = int
     *		2 = float
     *		3 = double
     *		4 = string
     */

    d[n_dialog-1].name = name;
    d[n_dialog-1].type = 0;
    d[n_dialog-1].value = value;
    d[n_dialog-1].width = 1;
}

void create_int_dialog_entry(char* name, int* value)
{
    extend_dialog_list();

    d[n_dialog-1].name = name;
    d[n_dialog-1].type = 1;
    d[n_dialog-1].value = (char*)value;
    d[n_dialog-1].width = 6;
}

void create_float_dialog_entry(char* name, float* value)
{
    extend_dialog_list();

    d[n_dialog-1].name = name;
    d[n_dialog-1].type = 2;
    d[n_dialog-1].value = (char*)value;
    d[n_dialog-1].width = 8;
}

void create_double_dialog_entry(char* name, double* value)
{
    extend_dialog_list();

    d[n_dialog-1].name = name;
    d[n_dialog-1].type = 3;
    d[n_dialog-1].value = (char*)value;
    d[n_dialog-1].width = 8;
}

void create_string_dialog_entry(char* name, char* value)
{
    extend_dialog_list();

    d[n_dialog-1].name = name;
    d[n_dialog-1].type = 4;
    d[n_dialog-1].value = value;
    d[n_dialog-1].width = 10;
}

void set_dialog_width(int width)
{
    /* Set width of the *most recently defined* dialog entry. */

    if (width > 0) d[n_dialog-1].width = width;
}

void initialize_graphics()
{
    /* Set up LUX. */

    set_default_font("7x13");	/* Better for MacX */
    lux_xinit();
}

void set_dialog_split(int s)
{
    if (s >= 0)
	split = s;
    else
	split = -1;
}

/* "Ez" user probably doesn't care about these,
    so hide them from view: */

static Window ez_dialog_win;
static int ez_xsize, ez_ysize;

void set_up_dialog(char* title, int xorigin, int yorigin)
{
    int spl;

    if (n_dialog > 0) {

	if (split >= 0)
	    spl = split;
	else
	    spl = (n_dialog+1)/2;

	initialize_graphics();
	initialize_dialog(&ez_dialog_win,
			  xorigin, yorigin, title,
			  &ez_xsize, &ez_ysize,
			  d, n_dialog, spl);
    }
}

int read_dialog_window()
{
    /*
     * Display the dialog box; lux_getevent() = 3 on return when the
     * user exits with "OK".  Then read the dialog box.
     */

    int status;

    lux_show_dialog(ez_dialog_win);

    if ((status = lux_getevent()) == 3)
	update_from_dialog(ez_dialog_win, d, n_dialog);

    return (status == 3);
}

void reset_dialog_window()
{
    reset_dialog_entries(ez_dialog_win, d, n_dialog);
}
