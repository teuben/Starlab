/*
 * interface.c
 * Biao Lu                     biao@eagle.drexel.edu
 *
 *
 * (Not entirely clear why this is needed... (SLWM))
 */

extern unsigned long lux_openwin();
extern unsigned long lux_rgb_pixel();
extern unsigned long lux_open_dialog();

unsigned long LUX_openwin(int a, int b, int c, int d)
{
    return(lux_openwin(a, b, c, d));
}
LUX_set_window_name(unsigned long a, char *s)
{
    return(lux_set_window_name(a,s));
}
LUX_setup_region(unsigned long a, float b,float c,float d,float e)
{
    return(lux_setup_region(a, b,c, d, e));
}
LUX_clear_current_region(unsigned long a)
{
    return(lux_clear_current_region(a));
}
LUX_setup_axis(unsigned long a, float b ,float c ,float d ,float e)
{
    return(lux_setup_axis(a, b ,c , d , e));
}
LUX_draw_linef(unsigned long a,float b ,float c ,float d,float e)
{
    return(lux_draw_linef(a, b , c , d, e));
}
LUX_draw_pointf(unsigned long a,float b,float c)
{
    return(lux_draw_pointf(a,b,c));
}
LUX_draw_rectanglef(unsigned long a,float b,float c,float d,float e)
{
    return(lux_draw_rectanglef(a, b, c, d,e));
}
LUX_draw_arcf(unsigned long a,float b,float c,float d,float e,float f,float g)
{
    return(lux_draw_arcf(a,b,c,d,e,f,g));
}
LUX_fill_arcf(unsigned long a,float b,float c,float d,float e,float f,float g)
{
    return(lux_fill_arcf(a,b,c,d,e,f,g));
}
LUX_draw_axis(unsigned long a)
{
    return(lux_draw_axis(a));
}
LUX_set_noupdate(unsigned long a)
{
    return(lux_set_noupdate(a));
}
LUX_next_keypress(unsigned long w, char* a, char* b, char* c, char* d)
{
    return(lux_next_keypress(w, a, b, c, d));
}
LUX_set_discard(unsigned long a)
{
    return(lux_set_discard(a));
}
LUX_set_save(unsigned long a)
{
    return(lux_set_save(a));
}
LUX_getevent()
{
    return(lux_getevent());
}
LUX_exit()
{
    return(lux_exit());
}
LUX_quick_exit()
{
    return(lux_quick_exit());
}
LUX_set_color(unsigned long a, long b)
{
    return(lux_set_color( a,  b));
}
LUX_set_bgcolor(unsigned long a, long b)
{
    return(lux_set_bgcolor( a,  b));
}
LUX_set_window_bgcolor(unsigned long a, long b)
{
    return(lux_set_window_bgcolor(a,  b));
}
unsigned long LUX_rgb_pixel(unsigned long a,float red,float green, float blue)
{
    return(lux_rgb_pixel(a,red , green,  blue));
}
unsigned long LUX_lookup_color(unsigned long a, char *b)
{
    return(lux_lookup_color(a,b));
}
LUX_draw_string(unsigned long  a, float b, float c, float d, char *e, char f)
{
    return(lux_draw_string(a, b,  c,  d, e, f));
}
LUX_draw_vstring(unsigned long a, float b, float c, float d, char *e, char f)
{
    return(lux_draw_vstring(a,  b,  c,  d, e,  f));
}
LUX_draw_image_string(unsigned long a, float b, float c, float d, char *e, char f)
{
    return(lux_draw_image_string(a, b, c,  d, e,  f));
}
LUX_check_keypress(unsigned long a,char b)
{
    return(lux_check_keypress(a,b));
}
LUX_check_buttonpress(unsigned long a)
{
    return(lux_check_buttonpress(a));
}
unsigned long LUX_open_dialog(int a, int b, int c, int d)
{
    return(lux_open_dialog(a, b, c,  d));
}
LUX_set_item(unsigned long a, int b, int c, int d, int e, int f, int g, char *h)
{
    return(lux_set_item(a, b, c, d, e, f, g, h));
}
LUX_draw_palette(unsigned long a)
{
    return(lux_draw_palette(a));
}
LUX_get_itemvalue(unsigned long a, int b, int c, int d, char *e)
{
    return(lux_get_itemvalue(a,  b,  c,  d, e));
}
LUX_update_itemvalue(unsigned long a, int b, int c, int d, char *e)
{
    return(lux_update_itemvalue(a,  b,  c, d, e));
}
LUX_clear_window(unsigned long a)
{
    return(lux_clear_window(a));
}
LUX_reset_window(unsigned long a)
{
    return(lux_reset_window(a));
}
LUX_update_fg(unsigned long a)
{
    return(lux_update_fg(a));
}
LUX_show_dialog(unsigned long a)
{
    return(lux_show_dialog(a));
}
LUX_reconvert_rcoord(unsigned long a, int b, int c, float *d, float *e)
{
    return(lux_reconvert_rcoord(a, b, c,d, e));
}
LUX_set_linestyle(unsigned long a, int b)
{
    lux_set_linestyle(a,b);
}
void LUX_set_nobuffer()
{
    lux_set_nobuffer();
}
