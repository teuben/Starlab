
/* Declarations for simple lux graphics. */

void set_box(int xo, int yo, int xs, int ys);
void nodraw_box();
void draw_box(float xmin_in, float xmax_in, char *xlabel,
              float ymin_in, float ymax_in, char *ylabel);
void add_label(char *label);
void set_color(char *color);
void plot(float *x, float *y, int n);
void move_to(float x, float y);
void draw_to(float x, float y);
void point(float x, float y, float size);
void get_mouse(float *x, float *y);
void pause(int time);
void clear_graphics();
void exit_graphics();

void crop_gfx();
void nocrop_gfx();
void wait_for_mouse();
int check_mouse();
