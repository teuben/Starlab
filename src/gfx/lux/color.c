/* color.c
 * Biao Lu                     biao@eagle.drexel.edu
 */

#include "win.h"
#include <string.h>

#define  ColormapFromDefault    0       /* Read Only Colormap */
#define  ColormapFromDefaultC   1       /* Changable Colormap */
#define  ColormapFromStandard   2       
#define  ColormapFromFile       3
#define  ImmutableColormap      4       /* Immutable Colormap */

extern Display           *display;
extern Visual            *visual;
extern int                screen;
extern XStandardColormap  map_info;
extern Colormap           colormap;
extern int                colormap_size;
extern char              *colorfile;
extern lux_wins          *get_currentwin();
unsigned   long           BLACK_COLOR,WHITE_COLOR;

XColor colors[256];              /* This should be associated with win
				    structure. But now all the windows
				    have same colormap  */
extern int    lux_colors;
extern int    lux_colormap;

lux_setup_colorfile(name)
char *name;
{
    if (colorfile != NULL) free(colorfile);
    colorfile = (char *) malloc((strlen(name)+1)*sizeof(char));
    colorfile[0] = 0;
    strcat(colorfile,name);
}

int lux_get_colormap_name(n)
char *n;
{
    if (n == (char *)NULL) return(ColormapFromDefault);
    if (!strcmp(n,"ColormapFromDefault"))  return(ColormapFromDefault);
    if (!strcmp(n,"ColormapFromDefaultC")) return(ColormapFromDefaultC);
    if (!strcmp(n,"ColormapFromStandard")) return(ColormapFromStandard);
    return(ColormapFromFile);
}

Colormap lux_setup_colormap(display, screen, visual)
Display *display;
int      screen;
Visual **visual;
{
    Colormap            colormap;
    Visual             *vis;
    XVisualInfo         visual_info, *vlist, vinfo_template, *v;
    int                 num_vis=0;
    unsigned int        default_depth;

    vlist = XGetVisualInfo(display, VisualNoMask, &vinfo_template, &num_vis);

    for (v = vlist; v < vlist + num_vis; v++) 
      if (v->visual == DefaultVisual(display,screen)) break;
    
    default_depth = DefaultDepth(display, screen);

    if (default_depth == 1) { /* Static Colormap */
	lux_colormap = ImmutableColormap;
	colormap_size = 2;
	BLACK_COLOR = BlackPixel(display,screen);
	WHITE_COLOR = WhitePixel(display,screen);
	*visual = v->visual;
	return(DefaultColormap(display, screen));	  
    }

    if (!XMatchVisualInfo(display, screen, default_depth, 
			  PseudoColor, &visual_info)) {
	if (!XMatchVisualInfo(display, screen, default_depth, 
			      DirectColor, &visual_info)) {
	    /*	  fprintf(stderr, "lux_msg: hardware colormap is immutable: cannot create new  colormap. Using default colormap\n");*/
	    lux_colormap = ImmutableColormap;
	    colormap_size = v->colormap_size;
	    BLACK_COLOR = BlackPixel(display,screen);
	    WHITE_COLOR = WhitePixel(display,screen);
	    *visual = v->visual;
	    return(DefaultColormap(display, screen));
	}
    }
    if (visual_info.visual != DefaultVisual(display,screen)) 
      for (v = vlist; v < vlist + num_vis; v++) 
	if (v->visual == visual_info.visual) break;

    colormap_size = v->colormap_size;
    *visual = v->visual;

    if (colormap_size>256) colormap_size = 256;

    lux_setup_map_info();

    BLACK_COLOR = BlackPixel(display,screen);
    WHITE_COLOR = WhitePixel(display,screen);

/*    fprintf(stderr,"%s %d\n",colorfile,lux_get_colormap_name(colorfile));*/

    switch(lux_colormap = lux_get_colormap_name(colorfile)) {

	case ColormapFromDefault:
	case ColormapFromDefaultC:
/*	  fprintf(stderr,"Using ColormapFromDefault\n");*/
	  return(DefaultColormap(display, screen)); break;
	case ColormapFromStandard:
	  colormap = XCreateColormap(display, RootWindow(display, screen), 
				     *visual, AllocAll);

	  if (colormap == DefaultColormap(display, screen)) {
	      fprintf(stderr,"lux_msg: hardware colormap is immutable: cannot create new colormap. Using default colormap\n");
	      lux_colormap = ColormapFromDefault;
	      return(DefaultColormap(display, screen));
	  }
	  lux_store_default_color(display,colormap,screen);
/*	  fprintf(stderr,"Using ColormapFromStandard\n");*/
	  return(colormap);
	  break;
        case ColormapFromFile:
	  colormap = XCreateColormap(display, RootWindow(display, screen), 
				     *visual, AllocAll);

	  if (colormap == DefaultColormap(display, screen)) {
	      fprintf(stderr,"lux_msg: hardware colormap is immutable: cannot create new colormap. Using default colormap\n");
	      lux_colormap = ColormapFromDefault;
	      return(DefaultColormap(display, screen));
	  }
	  if(lux_store_colorfromfile(display, colormap, screen, colorfile)) {
/*	      fprintf(stderr,"Using ColormapFromFile\n");*/
	      return(colormap);
	  }
	  else {
	      lux_colormap = ColormapFromDefault;
/*	      fprintf(stderr,"Using ColormapFromDefault\n");*/
	      return(DefaultColormap(display, screen));
	  }	          
	  break;
        default:
	  fprintf(stderr,"lux_msg: Unknown case in lux_setup_colormap\n");

      }
}

lux_store_colorfromfile(display, colormap, screen, filename)
Display            *display;
Colormap colormap;
int screen;
char *filename;
{
  int     ncells, i, j, k, l;
  XColor *exact_defs;
  unsigned char red[256],green[256],blue[256];
  FILE   *input;

  input = fopen(filename,"rb");

  if (input == (FILE *)NULL) return(0);
  ncells = colormap_size;
  exact_defs = (XColor *) calloc(sizeof(XColor), ncells); 

  fread(red,  1,256,input);
  fread(green,1,256,input);
  fread(blue, 1,256,input);

  for(i=0;i<ncells;i++) {
    exact_defs[i].red   = red[i]*256;
    exact_defs[i].green = green[i]*256;
    exact_defs[i].blue  = blue[i]*256;
    exact_defs[i].pixel = (unsigned long)i;
    exact_defs[i].flags = DoRed | DoGreen | DoBlue;
  }

  exact_defs[WhitePixel(display,screen)].red   = 65535;
  exact_defs[WhitePixel(display,screen)].green = 65535;
  exact_defs[WhitePixel(display,screen)].blue  = 65535;

  exact_defs[BlackPixel(display,screen)].red   = 0;
  exact_defs[BlackPixel(display,screen)].green = 0;
  exact_defs[BlackPixel(display,screen)].blue  = 0;
/*
  BLACK_COLOR = BlackPixel(display,screen);
  WHITE_COLOR = WhitePixel(display,screen);
*/
  XStoreColors (display, colormap, exact_defs, ncells);
  fclose(input);
  return(1);
}

lux_store_default_color(display,colormap,screen)
Display     *display;
Colormap     colormap;
int          screen;
{
  int     ncells, i, j, k, l;
  XColor *exact_defs;

  ncells = colormap_size;

  exact_defs = (XColor *)calloc(sizeof(XColor), ncells); 

  /* permute the levels of red, green and blue */
  l = map_info.base_pixel;
  for (k = 0; k <= map_info.red_max; k++) {
    for (j = 0; j <= map_info.green_max; j++) {
      for (i = 0; i <= map_info.blue_max; i++) {
	exact_defs[l].blue  = 65535 * i / map_info.blue_max;
	exact_defs[l].green = 65535 * j / map_info.green_max;
	exact_defs[l].red   = 65535 * k / map_info.red_max;
	exact_defs[l].pixel = (unsigned long)l;
	exact_defs[l].flags = DoRed | DoGreen | DoBlue;
	l++;
      }
    }
  }

  exact_defs[WhitePixel(display,screen)].red   = 65535;
  exact_defs[WhitePixel(display,screen)].green = 65535;
  exact_defs[WhitePixel(display,screen)].blue  = 65535;

  exact_defs[BlackPixel(display,screen)].red   = 0;
  exact_defs[BlackPixel(display,screen)].green = 0;
  exact_defs[BlackPixel(display,screen)].blue  = 0;
/*  
  BLACK_COLOR = BlackPixel(display,screen);
  WHITE_COLOR = WhitePixel(display,screen);
*/
  XStoreColors(display, colormap, exact_defs, ncells);
  free(exact_defs);
}

lux_setup_map_info()
{
        if (colormap_size == 256) {
          map_info.red_max   = (unsigned long)7;
          map_info.green_max = (unsigned long)7;
          map_info.blue_max  = (unsigned long)3;
 
          map_info.red_mult   = (unsigned long)32;
          map_info.green_mult = (unsigned long)4;
          map_info.blue_mult  = (unsigned long)1;
        }
        else if (colormap_size == 128) {
          map_info.red_max   = (unsigned long)7;
          map_info.green_max = (unsigned long)3;
          map_info.blue_max  = (unsigned long)3;
 
          map_info.red_mult   = (unsigned long)16;
          map_info.green_mult = (unsigned long)4;
          map_info.blue_mult  = (unsigned long)1;
        }
        else if (colormap_size == 16) {
          map_info.red_max   = (unsigned long)3;
          map_info.green_max = (unsigned long)1;
          map_info.blue_max  = (unsigned long)1;
 
          map_info.red_mult   = (unsigned long)4;
          map_info.green_mult = (unsigned long)2;
          map_info.blue_mult  = (unsigned long)1;
        }
        else if (colormap_size == 4) {
          map_info.red_max   = (unsigned long)3;
          map_info.green_max = (unsigned long)3;
          map_info.blue_max  = (unsigned long)3;
 
          map_info.red_mult   = (unsigned long)1;
          map_info.green_mult = (unsigned long)0;
          map_info.blue_mult  = (unsigned long)0;
        }
	else return(0);
        map_info.base_pixel = 0;
	return(1);
}

unsigned long lux_rgb_pixel(win,red, green, blue)
Window win;
float red,green,blue;
{

    register lux_wins *current;
    register int      i;
    unsigned long     plane_mask[1], pixels[1];
    int      redi, greeni, bluei;

    current = get_currentwin(win);


  if (red<0.000001 && green < 0.000001 && blue < 0.000001) return(BLACK_COLOR);
  if (red>0.999999 && green > 0.999999 && blue > 0.999999) return(WHITE_COLOR);

  switch(current->win.lux_colormap) {
      case ImmutableColormap:
        fprintf(stderr,"Colormap is immutable, use pixel value directly\n");
	return(BLACK_COLOR);
	break;
      case ColormapFromDefault:
        redi   = 65535L * (unsigned long)(0.5+red*map_info.red_max)/map_info.red_max;
        greeni = 65535L * (unsigned long)(0.5+green*map_info.green_max)/map_info.green_max;
        bluei  = 65535L * (unsigned long)(0.5+blue*map_info.blue_max)/map_info.blue_max;

	for(i=0;i<lux_colors;i++) 
	  if (redi   == colors[i].red   &&
	      greeni == colors[i].green &&
	      bluei  == colors[i].blue  )
	    return(colors[i].pixel);

	colors[lux_colors].red   = redi;
	colors[lux_colors].green = greeni;
	colors[lux_colors].blue  = bluei;
	colors[lux_colors].flags = DoRed | DoGreen | DoBlue;

	if(XAllocColor(current->win.display, current->win.colormap, &colors[lux_colors]))  {
	    lux_colors ++;
	    return(colors[lux_colors-1].pixel);    
	}   
	else {
	    fprintf(stderr,"Can't set read only color!\n");
	    return(BLACK_COLOR);
	}
	break;
      case ColormapFromDefaultC:
        redi   = 65535L * (unsigned long)(0.5+red*map_info.red_max)/map_info.red_max;
        greeni = 65535L * (unsigned long)(0.5+green*map_info.green_max)/map_info.green_max;
        bluei  = 65535L * (unsigned long)(0.5+blue*map_info.blue_max)/map_info.blue_max;

	for(i=0;i<lux_colors;i++) 
	  if (redi   == colors[i].red   &&
	      greeni == colors[i].green &&
	      bluei  == colors[i].blue  )
	    return(colors[i].pixel);

	if(XAllocColorCells(current->win.display, current->win.colormap, False, plane_mask, 0, pixels, 1)) {
	    colors[lux_colors].red   = redi;
	    colors[lux_colors].green = greeni;
	    colors[lux_colors].blue  = bluei;
	    colors[lux_colors].pixel = pixels[0];
	    colors[lux_colors].flags = DoRed | DoGreen | DoBlue;
	    XStoreColors(current->win.display, current->win.colormap, &colors[lux_colors], 1);
	    lux_colors++;
	    return(colors[lux_colors-1].pixel);
	}
	else {
	    fprintf(stderr,"Can't set red/write color!\n");
	    return(BLACK_COLOR);
	}
	break;
      case ColormapFromStandard:
	colors[lux_colors].pixel = map_info.base_pixel + 
	     (unsigned long)(0.5+red*map_info.red_max)*map_info.red_mult+
	     (unsigned long)(0.5+green*map_info.green_max)*map_info.green_mult+
	     (unsigned long)(0.5+blue*map_info.blue_max)*map_info.blue_mult;
        colors[lux_colors].red = 65535L * (unsigned long)(0.5+red*map_info.red_max)/map_info.red_max;
        colors[lux_colors].green = 65535L * (unsigned long)(0.5+green*map_info.green_max)/map_info.green_max;
        colors[lux_colors].blue = 65535L * (unsigned long)(0.5+blue*map_info.blue_max)/map_info.blue_max;
	colors[lux_colors].flags = DoRed | DoGreen | DoBlue;

	for(i=0;i<lux_colors;i++) 
	  if (colors[lux_colors].red   == colors[i].red   &&
	      colors[lux_colors].green == colors[i].green &&
	      colors[lux_colors].blue  == colors[i].blue  )
	    return(colors[i].pixel);

	lux_colors++;
        return(colors[lux_colors-1].pixel);
        break;
      case ColormapFromFile:
        fprintf(stderr,"Using colormap from file!!\n");
        return(BLACK_COLOR);
        break;
      default:
        fprintf(stderr,"Unknown case in rgb_pixel\n");
        return(BLACK_COLOR);
  }
}

unsigned long lux_lookup_color(win, name) 
Window    win;
char     *name;
{
    register lux_wins *current;
    register int      i;
    unsigned long     plane_mask[1], pixels[1];
    float    red, green, blue;
    XColor   db_def,  h_def;
    int      status;

    current = get_currentwin(win);

    if (!XLookupColor(current->win.display, current->win.colormap, name,
		      &db_def, &h_def)) 
      {fprintf(stderr,"Color %s not found.\n",name);return(BLACK_COLOR);}

    if (h_def.red   == 256 &&
	h_def.green == 256 &&
	h_def.blue  == 256) return(WHITE_COLOR);
    else if (h_def.red   == 0 &&
	     h_def.green == 0 &&
	     h_def.blue  == 0) return(BLACK_COLOR);


    switch(current->win.lux_colormap) {

	case ImmutableColormap:

	    /* Special treatment of immutable CM -- use pixel value directly. */

	    status = XAllocColor(current->win.display,
				 current->win.colormap,
				 &h_def);
	    return h_def.pixel;	/* ignore returned status for now! */
	    break;

	case ColormapFromDefault:

	    for(i = 0; i < lux_colors; i++) 
		if (h_def.red   == colors[i].red   &&
		    h_def.green == colors[i].green &&
		    h_def.blue  == colors[i].blue  )
		    return(colors[i].pixel);

	    h_def.flags = DoRed | DoGreen | DoBlue;
	    if(XAllocColor(current->win.display, current->win.colormap,
			   &h_def)) {
		colors[lux_colors].red   = h_def.red;
		colors[lux_colors].green = h_def.green;
		colors[lux_colors].blue  = h_def.blue;
		colors[lux_colors].flags = DoRed | DoGreen | DoBlue;
		colors[lux_colors].pixel = h_def.pixel;
		lux_colors++;
		return(colors[lux_colors-1].pixel);  
	    }
	    else {
		fprintf(stderr,"Can't set read only color %s!\n",name);
		return(BLACK_COLOR);
	    }
	    break;

	case ColormapFromDefaultC:

	    for(i=0;i<lux_colors;i++) 
		if (h_def.red   == colors[i].red   &&
		    h_def.green == colors[i].green &&
		    h_def.blue  == colors[i].blue  )
		    return(colors[i].pixel);

	    h_def.flags = DoRed | DoGreen | DoBlue;
	    if(XAllocColorCells(current->win.display, current->win.colormap,
				False, plane_mask, 0, pixels, 1)) {
		colors[lux_colors].red   = h_def.red;
		colors[lux_colors].green = h_def.green;
		colors[lux_colors].blue  = h_def.blue;
		colors[lux_colors].pixel = pixels[0];
		colors[lux_colors].flags = DoRed | DoGreen | DoBlue;
		XStoreColors(current->win.display, current->win.colormap,
			     &colors[lux_colors], 1);
		lux_colors++;
		return(colors[lux_colors-1].pixel);
	    }
	    else {
		fprintf(stderr,"Can't set red/write color!\n");
		return(BLACK_COLOR);
	    }
	    break;

      case ColormapFromStandard:

	  red   = ((float)h_def.red  )/65535.0;
	  green = ((float)h_def.green)/65535.0;
	  blue  = ((float)h_def.blue )/65535.0;
	  colors[lux_colors].pixel = map_info.base_pixel + 
	      (unsigned long)(0.5+red*map_info.red_max)*map_info.red_mult+
	      (unsigned long)(0.5+green*map_info.green_max)*map_info.green_mult+
	      (unsigned long)(0.5+blue*map_info.blue_max)*map_info.blue_mult;
	  colors[lux_colors].red
	      = 65535L * (unsigned long)(0.5+red*map_info.red_max)
						/map_info.red_max;
	  colors[lux_colors].green
	      = 65535L * (unsigned long)(0.5+green*map_info.green_max)
						/map_info.green_max;
	  colors[lux_colors].blue
	      = 65535L * (unsigned long)(0.5+blue*map_info.blue_max)
						/map_info.blue_max;
	  colors[lux_colors].flags = DoRed | DoGreen | DoBlue;

	  for(i = 0 ; i < lux_colors ; i++) 
	      if (colors[lux_colors].red   == colors[i].red   &&
		  colors[lux_colors].green == colors[i].green &&
		  colors[lux_colors].blue  == colors[i].blue  )
		  return(colors[i].pixel);

	  lux_colors++;
	  return(colors[lux_colors-1].pixel);
	  break;

      case ColormapFromFile:

	  fprintf(stderr, "Using colormap from file!!\n");
	  return(BLACK_COLOR);
	  break;

      default:

	  fprintf(stderr, "Unknown case in rgb_pixel\n");
/*        return(BLACK_COLOR);	*/

    }
    return(BLACK_COLOR);
}


lux_rotate_colormap(win,flag)
Window win;
int    flag;
{
    register lux_wins *current;
    register unsigned long pixel;
    XColor  *exact_defs;
    int      i;

    if (flag == 0) return(0);

    current = get_currentwin(win);

    switch(current->win.lux_colormap) {
	case ImmutableColormap:
	case ColormapFromDefault:
	fprintf(stderr,"Colormap can't be rotate\n");
	break;
        case ColormapFromFile:
	exact_defs = (XColor *) calloc(sizeof(XColor), 258); 
	for(i=1;i<257;i++) {
	    exact_defs[i].pixel = i-1;
	    exact_defs[i].flags = DoRed | DoGreen | DoBlue;
	}
	XQueryColors(current->win.display, current->win.colormap, 
		      &(exact_defs[1]), 256);
	exact_defs[ 0 ].pixel = exact_defs[256].pixel;
	exact_defs[257].pixel = exact_defs[ 1 ].pixel;
	if (flag > 0) {
	    for(i=0;i<257;i++)	exact_defs[i].pixel = exact_defs[i+1].pixel;
/*
	    exact_defs[ 0 ].pixel = exact_defs[256].pixel;
	    exact_defs[257].pixel = exact_defs[ 1 ].pixel;
	    if ( BLACK_COLOR > WHITE_COLOR ) {
		exact_defs[BLACK_COLOR].pixel = exact_defs[BLACK_COLOR+1].pixel;
		exact_defs[BLACK_COLOR+1].pixel = BLACK_COLOR;
		exact_defs[WHITE_COLOR].pixel = exact_defs[WHITE_COLOR+1].pixel;
		exact_defs[WHITE_COLOR+1].pixel = WHITE_COLOR;		
	    }
	    else {
		exact_defs[WHITE_COLOR].pixel = exact_defs[WHITE_COLOR+1].pixel;
		exact_defs[WHITE_COLOR+1].pixel = WHITE_COLOR;		
		exact_defs[BLACK_COLOR].pixel = exact_defs[BLACK_COLOR+1].pixel;
		exact_defs[BLACK_COLOR+1].pixel = BLACK_COLOR;
	    }	    
	    exact_defs[256].pixel = exact_defs[ 0 ].pixel;
	    exact_defs[ 1 ].pixel = exact_defs[257].pixel;
*/
	}
	else {
	    for(i=257;i>0;i--) exact_defs[i].pixel = exact_defs[i-1].pixel;
/*
	    exact_defs[ 0 ].pixel = exact_defs[256].pixel;
	    exact_defs[257].pixel = exact_defs[ 1 ].pixel;
	    if ( BLACK_COLOR > WHITE_COLOR ) {
		exact_defs[BLACK_COLOR+2].pixel = exact_defs[BLACK_COLOR].pixel;
		exact_defs[BLACK_COLOR].pixel = BLACK_COLOR;
		exact_defs[WHITE_COLOR].pixel = exact_defs[WHITE_COLOR+1].pixel;
		exact_defs[WHITE_COLOR+1].pixel = WHITE_COLOR;
		if ( WHITE_COLOR == 0 ) 
		  
	    }
	    else {
		exact_defs[WHITE_COLOR].pixel = exact_defs[WHITE_COLOR+1].pixel;
		exact_defs[WHITE_COLOR+1].pixel = WHITE_COLOR;		
		exact_defs[BLACK_COLOR].pixel = exact_defs[BLACK_COLOR+1].pixel;
		exact_defs[BLACK_COLOR+1].pixel = BLACK_COLOR;
	    }	    
	    exact_defs[256].pixel = exact_defs[ 0 ].pixel;
	    exact_defs[ 1 ].pixel = exact_defs[257].pixel;
*/
	}
	
	XStoreColors(current->win.display, current->win.colormap, 
		      &(exact_defs[1]), 256);
	free(exact_defs);
/*	fprintf(stderr,"Colormap rotated\n");*/
/*	XFlush(current->win.display);*/
	break;
	case ColormapFromDefaultC:
	case ColormapFromStandard:
	if (flag > 0) {
	    pixel = colors[0].pixel;
	    for(i=0;i<lux_colors-2;i++) {
		colors[i].pixel = colors[i+1].pixel;
	    }
	    colors[lux_colors-1].pixel = pixel;
	}
	else {
	    pixel = colors[lux_colors-1].pixel;
	    for(i=lux_colors-2;i>=0;i--) {
		colors[i+1].pixel = colors[i].pixel;
	    }
	    colors[0].pixel = pixel;
	}
	XStoreColors(current->win.display,current->win.colormap,colors,lux_colors);
/*	fprintf(stderr,"Colormap rotated\n");*/
/*	XFlush(current->win.display);*/
	return(1);
	default:
	fprintf(stderr,"Unknown case in rotate colormap\n");
	break;
    }
    return(0);
}

lux_get_colorinfo(size,im_colormap) /* if im_colormap is 1, the colormap
				       is immutable */
int *size, *im_colormap;
{
    *size = colormap_size;
    if (lux_colormap == ImmutableColormap) *im_colormap = 1;
    else *im_colormap = 0;
}

lux_set_window_colormap(win, name)	/* Return the number of colors! */
Window win;
char  *name;
{
    register lux_wins *current;
    Colormap cm;

    current = get_currentwin(win);

    /* Reset to default colormap */
    if (!strcmp(name,"LuxDefaultColormap")) {
	if (current->win.colormap != colormap) {
	    XFreeColormap(current->win.display,current->win.colormap);
	    if (current->win.colormapfile != NULL) 
	      free(current->win.colormapfile);
	    if (colorfile != NULL) {
		current->win.colormapfile = (char *)malloc(strlen(colorfile)+1);
		current->win.colormapfile[0] = 0;
		strcat(current->win.colormapfile,colorfile);
	    }
	    else current->win.colormapfile = NULL;
	    
	    current->win.lux_colormap = lux_colormap;
	    current->win.colormap = colormap;
	    XSetWindowColormap(current->win.display,win,colormap);
	}
	return 16;
    }
    


    if (lux_colormap == ImmutableColormap) {
	fprintf(stderr,"Error: The colormap is immutable.\n");
	return(0);
    }

    cm = XCreateColormap(current->win.display, 
			 RootWindow(current->win.display, screen), 
			 visual, AllocAll);

    if (lux_store_colorfromfile(current->win.display, cm, screen, name) == 0) {
	fprintf(stderr,"Error: File not found.\n");
	XFreeColormap(current->win.display, cm);
	return(0);
    }

    if (current->win.colormap != colormap) /* Not default colormap */
      XFreeColormap(current->win.display,current->win.colormap);

    if (current->win.colormapfile != NULL) free(current->win.colormapfile);
	
    current->win.colormapfile = (char *)malloc(strlen(name)+1);
    current->win.colormapfile[0] = 0;
    strcat(current->win.colormapfile,name);

    current->win.lux_colormap = ColormapFromFile;
    current->win.colormap = cm;

    XSetWindowColormap(current->win.display,win,cm);
    return 256;
}




