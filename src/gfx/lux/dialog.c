/* dialog.c
 * Biao Lu                     biao@eagle.drexel.edu
 */

#include "win.h"

extern lux_wins      *get_currentwin();
extern lux_wins      *lux_allocwin();

extern XFontStruct  *font_info;
extern unsigned long BLACK_COLOR, WHITE_COLOR;
extern unsigned long lux_get_fgcolor();
extern unsigned long lux_get_bgcolor();


lux_highlight(win)
Window win;
{
    register lux_wins *current, *parent;
    XGCValues values;


    current = get_currentwin(win);
    parent  = get_currentwin(current->win.parent);

    XGetGCValues(parent->win.display, parent->win.gc, 
		 GCLineWidth | GCLineStyle | GCForeground | GCBackground |
		 GCCapStyle | GCJoinStyle, 
		 &values);

    XSetLineAttributes(parent->win.display, parent->win.gc, 1, LineSolid, 
		      CapButt, JoinBevel);

    XSetForeground(parent->win.display, parent->win.gc, BLACK_COLOR);

    XDrawRectangle(parent->win.display, parent->win.window, parent->win.gc,
		   current->win.x-1, current->win.y-1, 
		   current->win.width+3, current->win.height+3);

    XSetLineAttributes(parent->win.display, parent->win.gc, 
		      values.line_width, values.line_style, 
		      values.cap_style, values.join_style);

    XSetForeground(parent->win.display, parent->win.gc, values.foreground);

}

lux_unhighlight(win)
Window win;
{
    register lux_wins *current, *parent;
    XGCValues values;


    current = get_currentwin(win);
    parent  = get_currentwin(current->win.parent);

    XGetGCValues(parent->win.display, parent->win.gc, 
		 GCLineWidth | GCLineStyle | GCForeground | GCBackground |
		 GCCapStyle | GCJoinStyle, 
		 &values);

    XSetLineAttributes(parent->win.display, parent->win.gc, 1, LineSolid, 
		      CapButt, JoinBevel);

    XSetForeground(parent->win.display, parent->win.gc, values.background);

    XDrawRectangle(parent->win.display, parent->win.window, parent->win.gc,
		   current->win.x-1, current->win.y-1, 
		   current->win.width+3, current->win.height+3);

    XSetLineAttributes(parent->win.display, parent->win.gc, 
		      values.line_width, values.line_style, 
		      values.cap_style, values.join_style);

    XSetForeground(parent->win.display, parent->win.gc, values.foreground);

}

Window lux_open_dialog(/*win*/x, y, width, height)
unsigned int x, y, width, height;
{
    register lux_wins *mywin;
    static   int       status = 0;

    if ( !status ) {
      lux_xinit();
      status = 1;
    }

    if ((mywin = lux_allocwin()) == (lux_wins *)NULL ) return((Window)NULL);
    
    mywin->win.type = DIALOG_WINDOW;
    mywin->win.window_name = "Dialog Window";

    mywin->win.x = x;
    mywin->win.y = y;
    mywin->win.width  = width;
    mywin->win.height = height;
    mywin->win.old_width   = width;
    mywin->win.old_height  = height;
    mywin->win.user_width  = width;
    mywin->win.user_height = height;
    mywin->win.xresizefactor = 1.0;
    mywin->win.yresizefactor = 1.0;

    mywin->win.fxorg = 0.0;
    mywin->win.fyorg = 0.0;
    mywin->win.fxsize = 10.0;
    mywin->win.fysize = 10.0;
    mywin->win.xfactor = 1.0;
    mywin->win.yfactor = 1.0;
    mywin->win.xmin = 0.0;
    mywin->win.xmax = width;
    mywin->win.ymin = 0.0;
    mywin->win.ymax = height;
    mywin->win.xsize = width;
    mywin->win.ysize = height;
    mywin->win.xorg = 0;
    mywin->win.yorg = 0;
    
    lux_createwin(&mywin->win);

    return(mywin->win.window);
}

lux_set_item(dialog, num, type, subtype, x, y, length, def_value) 
Window dialog;
int num, subtype, type, x, y, length;     /* Only string item allowed */
char *def_value;                          /* length < 255 */
{
    register lux_wins *mywin;
    register lux_wins *diawin;
    static   int       status = 0;



    switch(type) {
      
    case TEXT_WINDOW:
      lux_draw_string(dialog, (float)x, (float)y, 0.0, def_value, -1);
      break;
      
    case BUTTON_WINDOW:
    case INPUT_WINDOW:

      if ((mywin = lux_allocwin()) == (lux_wins *)NULL ) return((Window)NULL);
      mywin->win.type = type;
      mywin->win.subtype = subtype;
      mywin->win.parent = dialog;
      mywin->win.serial = num;
      
      diawin = get_currentwin(dialog);

      mywin->win.x = x;
      mywin->win.y = diawin->win.height - y - 
	             (font_info->ascent + font_info->descent);
      if (subtype == CHECK_BUTTON) 
	mywin->win.width=mywin->win.height=font_info->ascent+font_info->descent;
      else {
	mywin->win.width  = length * XTextWidth(font_info, "X", 1)+4;
	mywin->win.height = font_info->ascent + font_info->descent + 4;
      }
      mywin->win.old_width   = mywin->win.user_width  = mywin->win.width;
      mywin->win.old_height  = mywin->win.user_height = mywin->win.height;
      mywin->win.xresizefactor = 1.0;
      mywin->win.yresizefactor = 1.0;

      mywin->win.fxorg = 0.0;
      mywin->win.fyorg = 0.0;
      mywin->win.fxsize = 10.0;
      mywin->win.fysize = 10.0;
      mywin->win.xfactor = 1.0;
      mywin->win.yfactor = 1.0;
      mywin->win.xmin = 0.0;
      mywin->win.xmax = mywin->win.width;
      mywin->win.ymin = 0.0;
      mywin->win.ymax = mywin->win.height;
      mywin->win.xsize = mywin->win.width;
      mywin->win.ysize = mywin->win.height;
      mywin->win.xorg = 0;
      mywin->win.yorg = 0;

      lux_createwin(&mywin->win);

      if (type == BUTTON_WINDOW) {
	if (subtype == OK_BUTTON || subtype == OK_KEEP_BUTTON
	    || subtype == CANCEL_BUTTON) 
	  lux_draw_string(mywin->win.window, mywin->win.width/2.0, 
			  4.0, 0.0, def_value, 0);
	else if (subtype == CHECK_BUTTON) {
	  if (def_value[0] == 1) {
	    lux_draw_line(mywin->win.window, 0, 0, 
		      mywin->win.width, mywin->win.height);
	    lux_draw_line(mywin->win.window, 0, mywin->win.height, 
		      mywin->win.width, 0);
	  }
	  mywin->win.msg[1] = mywin->win.msg[0] = def_value[0]; 
	                                  /* Save the old one */
	}	  
      }
      else if (type == INPUT_WINDOW) {
	lux_draw_string(mywin->win.window, 2.0, 4.0, 0.0, def_value, -1);
	{
	  char *new;
	  int i;
	  new = (char *)malloc(255); 
	  new[0] = (char)NULL; mywin->win.msg[0] = (char)NULL;
	  for(i=0;i<sizeof(float)*3+1;i++) 
	    new[i] = mywin->win.data->data.b[i];
	  new[i] = (char)NULL;
	  strcat(&new[i],&(mywin->win.data->data.b[i]));
	  strcat(mywin->win.msg,&(mywin->win.data->data.b[i]));
	  free(mywin->win.data->data.b);
	  mywin->win.data->data.b = new;
	}	  
      }
      break;
      
    default:
      fprintf(stderr, "\n Something not expected in create item!\n ");
      break;
    }

    return(1);
}

lux_draw_palette(dialog)  /* Only for 512x512 dialog window */
Window dialog;
{
    register lux_wins *current;
    register unsigned long i;
    unsigned long c;

    current = get_currentwin(dialog);

    c = lux_get_fgcolor(dialog);
    for(i=0;i<256;i++) { 
      lux_set_color(dialog, i); 
      lux_draw_line(dialog, 0, (int)(i*2),   20, (int)(i*2));
      lux_draw_line(dialog, 0, (int)(i*2+1), 20, (int)(i*2+1));
    }
    lux_set_color(dialog,c);
    lux_draw_line(dialog, 20, 0, 20, 512);
}

lux_update_itemvalue(dialog, num, type, subtype, value)
Window dialog;
int num, type, subtype;
char *value;
{
    register lux_wins *current, *old;

    old = current = get_currentwin(dialog);

    while(current->win.parent != dialog || current->win.serial != num ||
	  current->win.type != type || current->win.subtype != subtype ) {
      current = current->next;
      if (old == current) {
	  fprintf(stderr, "lux_update_itemvalue: dialog item not found\n");
	  return(0);
      }
    }

    switch(current->win.type) {
    case INPUT_WINDOW:
      current->win.msg[0] = 0;
      strcat(current->win.msg, value);
      current->win.data->data.b[sizeof(float)*3+1] = 0;
      strcat(&(current->win.data->data.b[sizeof(float)*3+1]), value);
      redraw(current->win.window,1);
      break;
    case BUTTON_WINDOW:
      if (current->win.subtype == CHECK_BUTTON) {
	if (current->win.msg[0] == 0 && value[0] == 1) {
	    lux_draw_line(current->win.window, 0, 0, 
		      current->win.width, current->win.height);
	    lux_draw_line(current->win.window, 0, current->win.height, 
		      current->win.width, 0);
	}
	else if (current->win.msg[0] == 1 && value[0] == 0) 
	    lux_reset_window(current->win.window);
	current->win.msg[1] = current->win.msg[0] = value[0];
      }
      break;
    default:
      break;
    }
    return(1);
}

lux_get_itemvalue(dialog, num, type, subtype, value)
Window dialog;
int num, type, subtype;
char *value;
{
    register lux_wins *current, *old;

    old = current = get_currentwin(dialog);

    value[0] = (char)NULL;

    while(current->win.parent != dialog || current->win.serial != num ||
	  current->win.type != type || current->win.subtype != subtype ) {
      current = current->next;
      if (old == current) {
	  fprintf(stderr, "lux_get_itemvalue: dialog item not found\n");
	  return(0);
      }
    }

    switch(current->win.type) {
    case INPUT_WINDOW:
      strcat(value, current->win.msg);
      break;
    case BUTTON_WINDOW:
      if (current->win.subtype == CHECK_BUTTON) 
	value[0] = current->win.msg[0];
      break;
    default:
      break;
    }
    return(1);
}


lux_show_dialog(win)
Window win;
{
    register lux_wins *current, *old;

    old = current = get_currentwin(win);

    while(current->win.type != DIALOG_WINDOW /*||
	  current->win.serial != win*/) {
      current = current->next;
      if (old == current) {
	fprintf(stderr, "No dialog window opened!\n");
	return(0);
      }

    }

    XMapWindow(current->win.display, current->win.window);
/*    XBell(current->win.display, 100); */
    return(1);
}

lux_hide_dialog(dia)
Window dia;
{
    register lux_wins *current;

    current = get_currentwin(dia);

    XUnmapWindow(current->win.display, current->win.window);
    return(1);
}

lux_ok_data(dialog)
Window dialog;
{

    register lux_wins *current, *old;

    old = current = get_currentwin(dialog);

    current = current->next;
    while(old != current) {
      if (current->win.parent != dialog ||
	  (current->win.type  != INPUT_WINDOW &&
	   current->win.type  != BUTTON_WINDOW)) {
	current = current->next;continue; 
      }
      if (current->win.type   == INPUT_WINDOW ) 
	
	  strcpy(current->win.msg,
		 &(current->win.data->data.b[3*sizeof(float)+1]));
      else 
	switch(current->win.subtype) {
	case CHECK_BUTTON:
	  current->win.msg[1] = current->win.msg[0];
	  break;
	case OK_BUTTON:
	case OK_KEEP_BUTTON:
	case CANCEL_BUTTON:
	default:
	  break;
	}
      current = current->next;
    }  
}

lux_cancel_data(dialog)
Window dialog;
{
    register lux_wins *current, *old;

    old = current = get_currentwin(dialog);

    current = current->next;
    while(old != current) {
      if (current->win.parent != dialog ||
	  (current->win.type  != INPUT_WINDOW &&
	   current->win.type  != BUTTON_WINDOW)) {
	current = current->next; continue; 
      }
      if (current->win.type == INPUT_WINDOW) {
	  strcpy(&(current->win.data->data.b[3*sizeof(float)+1]),
		 current->win.msg);
      }
      else 
	switch(current->win.subtype) {
	case CHECK_BUTTON:
	  if (current->win.msg[1] == 0 && current->win.msg[0] == 1) 
	    lux_reset_window(current->win.window);
	  else if (current->win.msg[1] == 1 && current->win.msg[0] == 0) {
	    lux_draw_line(current->win.window, 0, 0, 
		      current->win.width, current->win.height);
	    lux_draw_line(current->win.window, 0, current->win.height, 
		      current->win.width, 0);
	  }
	  current->win.msg[0] = current->win.msg[1];
	  break;
	case OK_BUTTON:
	case OK_KEEP_BUTTON:
	case CANCEL_BUTTON:
	default:
	  break;
	}
      current = current->next;
    }   
}


