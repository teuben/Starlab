
#include "sdyn3.h"
#include "xhrdplot.h"
#include "single_star.h"

// Open windows for plotting and instructions, set up color schemes, etc.
static float linespacing;    // Set in xstarplot.  Too many references


void initialize_hrd_graphics(float r_factor,         // (input) scale factor
                         bool b_flag,                // (input) color scheme
                         unsigned long& hrd_win,     // ID of plot window
                         unsigned long& hrd_instr,   // ID of status window
                         unsigned long& hrd_colwin,  // ID of color window
                         unsigned long* c_index,     // color by index 
			 unsigned long* c_energy,    // color by energy
                         int& win_size,              // width of plot window
                         int& hrd_xorigin, int& hrd_yorigin) // origin for dialog box
{
    // This is a kludge to try to get the font size right...
    // Note that the same font is used everywhere.

    // The layout works best for 0.6 < r_factor < 1.3.

    if (r_factor >= 1.0) {
        set_default_font("9x15");
        linespacing = 1.2;
    } else if (r_factor >= 0.8) {
        set_default_font("7x13");
        if (r_factor >= 0.9)
            linespacing = 1.15;
        else
            linespacing = 1.1;
    } else if (r_factor >= 0.65) {
        set_default_font("6x12");
        if (r_factor >= 0.725)
            linespacing = 1.05;
        else
            linespacing = 1.0;
    } else if (r_factor >= 0.5) {
        set_default_font("6x10");
        if (r_factor >= 0.575)
            linespacing = 1.05;
        else
            linespacing = 1.0;
    } else if (r_factor >= 0.4) {
        set_default_font("6x9");
        linespacing = 1.0;
    } else {
        set_default_font("5x8");
        linespacing = 1.0;
    }

    // The twm title bar is ~27-30 pixels (*independent* of scaling),
    // but the bar doesn't count in the location of the window!
    // Mac windows need more space (frame at bottom too).
    // Assume that scale factors of 0.75 or less mean we are
    // using a Mac or PC display.  Also, for r_factor > 1.1, move
    // the  "colwin" box up and over.

    // Main window:

    int title = 27;
    int bottom = 0;
    if (r_factor <= 0.75) {
        title = 20;
        bottom = 20;
    }

    int xsize = _R_(400);
    int ysize = _R_(400);
    hrd_xorigin = xsize + 6;
    hrd_yorigin = 0;
    hrd_win = lux_openwin(hrd_xorigin, hrd_yorigin, xsize, ysize);
    lux_set_window_name(hrd_win, "Hertzsprung-Russel");

    win_size = xsize;

    // Status window:

    if (r_factor <= 1)
        hrd_yorigin += ysize + title + bottom;
    else
        hrd_yorigin += ysize + 2*title + bottom;

    // (Extra factor of 2 seems to improve layout...)

    ysize = _R_(260);
    if (r_factor <= 0.6) {
        xsize = (int) (xsize * 1.25);
        ysize = (int) (ysize * 1.1);
    }
    hrd_instr = lux_openwin(hrd_xorigin, hrd_yorigin, xsize, ysize);
    lux_set_window_name(hrd_instr, "Status & Commands");

    // Color window:

    if (r_factor <= 1.1 && r_factor > 0.75) {
        hrd_yorigin += ysize + title + bottom;
    } else {
        int edge = 3;
        if (r_factor <= 0.6) xsize = (int) (xsize / 1.25);
        hrd_xorigin += xsize + edge;
        hrd_yorigin = 0;
    }
    if (r_factor <= 0.6) {
        xsize = (int) (xsize * 1.25);
        ysize = (int) (ysize * 1.1);
    }
    if (r_factor <= 0.6) ysize =  (int) (ysize * 1.1);
    hrd_colwin = lux_openwin(hrd_xorigin, hrd_yorigin, xsize, ysize);
    lux_set_window_name(hrd_colwin, "Color coding");

    init_colors(hrd_win, c_energy, c_index, b_flag);

    lux_set_noupdate(hrd_win);
    lux_set_bgcolor(hrd_win,c_energy[background_color]);
    lux_set_window_bgcolor(hrd_win,c_energy[background_color]);
    lux_set_color(hrd_win,c_energy[default_color]);

    lux_set_noupdate(hrd_colwin);
    lux_set_bgcolor(hrd_colwin,c_energy[background_color]);
    lux_set_window_bgcolor(hrd_colwin,c_energy[background_color]);
    lux_set_color(hrd_colwin,c_energy[default_color]);

    lux_set_noupdate(hrd_instr);
    lux_set_bgcolor(hrd_instr,c_energy[background_color]);
    lux_set_window_bgcolor(hrd_instr,c_energy[background_color]);
    lux_set_color(hrd_instr,c_energy[default_color]);

    lux_setup_region(hrd_win, 1.0, 1.0, 8.0, 8.0);

    // Set up origin of dialog box.

    if (r_factor <= 1.1 && r_factor > 0.75) {
        int edge = 3;
        hrd_xorigin += xsize + edge;
        hrd_yorigin = 0;
    } else
        hrd_yorigin += ysize + title + bottom;
}

void set_hrd_limits(float& xmin, float& xmax,
                    float& ymin, float& ymax)
{
    xmin = -1.0;
    xmax = 7.0;
    ymin = -1;
    ymax = 7;
//    xmin = -3.0;
//    xmax = 10.0;
//    ymin = -5;
//    ymax = 8;
}

void derive_stellar_type(dyn* bj, dyn* bi, char& c)
{
    c = default_color;
    real rmin = VERY_LARGE_NUMBER;
    dynptr bmin = bi;

    stellar_type type = NAS;
    real t_cur, t_rel, m_rel, m_env, m_core, T_eff, L_eff, p_rot, b_fld;
    real co_core=0;
    t_cur=t_rel=m_rel=m_env=m_core=T_eff=L_eff=p_rot=b_fld=0;
    story *s = bi->get_starbase()->get_star_story();
    extract_story_chapter(type, t_cur, t_rel,
                          m_rel, m_env, m_core, co_core,
			  T_eff, L_eff, p_rot, b_fld, *s);
//cerr<<"star: "<<type<<" " <<t_rel<<" " <<m_rel<<" " <<m_env<<" "
//<<m_core<<" "<<T_eff<< " "<< L_eff<<endl;
    switch(summarize_stellar_type(type)) {
       case ZAMS          : c = bound_binary;
                            break;
       case Early_Giant   : c = unbound_binary;
                            break;
       case Late_Giant    : c = unbound_single;
                            break;
       case Helium_Remnant: c = default_color;
                            break;
       case Inert_Remnant : c = bound_single;
                            break;
       default            : c = default_color;
       break;
    }
}

void show_hrd_color_scheme(unsigned long co, unsigned long *c_e, unsigned long *c_i,
                       float r, char ce, bool b_flag, int u)
{
    lux_clear_window(co);

    if (ce==1) {
        show_instructions(co, r, "Color by stellar type:\n",0,u);
        lux_set_color(co, c_e[bound_binary]);
        show_instructions(co, r, "  Green  = main sequence star\n",3,u);
        lux_set_color(co, c_e[unbound_binary]);
        show_instructions(co, r, "  Gold   = sub giant\n",5,u);
        lux_set_color(co, c_e[unbound_single]);
        show_instructions(co, r, "  Red    = red giant\n",4,u);
        lux_set_color(co, c_e[default_color]);
        if (b_flag)
            show_instructions(co, r,
                              "  White  = white dwarf of helium star\n",2,u);
        else
            show_instructions(co, r,
                              "  Black  = white dwarf of helium star\n",2,u);
        lux_set_color(co, c_e[bound_single]);
        show_instructions(co, r,
                              "  Blue   = neutron star or black hole\n",1,u);
        lux_set_color(co, c_e[default_color]);

    } else if(ce==2) {
	show_instructions(co, r, "Color by energy:\n",0,u);
	lux_set_color(co, c_e[default_color]);
	if (b_flag) 
	    show_instructions(co, r,
			      "  White  = default (bound to cluster)\n",1,u);
	else 
	    show_instructions(co, r,
			      "  Black  = default (bound to cluster)\n",1,u);
	lux_set_color(co, c_e[bound_single]);
	show_instructions(co, r, "  Blue   = bound to nearest neighbor\n",2,u);
	lux_set_color(co, c_e[bound_binary]);
	show_instructions(co, r, "  Green  = bound binary\n",3,u);
	lux_set_color(co, c_e[unbound_single]);
	show_instructions(co, r, "  Red    = unbound single star\n",4,u);
	lux_set_color(co, c_e[unbound_binary]);
	show_instructions(co, r, "  Gold   = unbound binary\n",5,u);
	lux_set_color(co, c_e[default_color]);
     }
     else {

        show_instructions(co, r, "Color by index:\n",0,u);
        lux_set_color(co, c_i[1]);

        if (b_flag) {
            format_and_show_instructions(co, r, c_i,  1, 0, "white", 1, u);
            format_and_show_instructions(co, r, c_i,  2, 1, RV_COLOR2, 1, u);
            format_and_show_instructions(co, r, c_i,  3, 2, RV_COLOR3, 1, u);
            format_and_show_instructions(co, r, c_i,  4, 0, RV_COLOR4, 2, u);
            format_and_show_instructions(co, r, c_i,  5, 1, RV_COLOR5, 2, u);
            format_and_show_instructions(co, r, c_i,  6, 2, RV_COLOR6, 2, u);
            format_and_show_instructions(co, r, c_i,  7, 0, RV_COLOR7, 3, u);
            format_and_show_instructions(co, r, c_i,  8, 1, RV_COLOR8, 3, u);
            format_and_show_instructions(co, r, c_i,  9, 2, RV_COLOR9, 3, u);
            format_and_show_instructions(co, r, c_i, 10, 0, RV_COLORa, 4, u);
            format_and_show_instructions(co, r, c_i, 11, 1, RV_COLORb, 4, u);
            format_and_show_instructions(co, r, c_i, 12, 2, RV_COLORc, 4, u);
            format_and_show_instructions(co, r, c_i, 13, 0, RV_COLORd, 5, u);
            format_and_show_instructions(co, r, c_i, 14, 1, RV_COLORe, 5, u);
            format_and_show_instructions(co, r, c_i, 15, 2, RV_COLORf, 5, u);
            format_and_show_instructions(co, r, c_i, 16, 0, RV_COLORg, 6, u);
        } else {
            format_and_show_instructions(co, r, c_i,  1, 0, "black", 1, u);
            format_and_show_instructions(co, r, c_i,  2, 1, NV_COLOR2, 1, u);
            format_and_show_instructions(co, r, c_i,  3, 2, NV_COLOR3, 1, u);
            format_and_show_instructions(co, r, c_i,  4, 0, NV_COLOR4, 2, u);
            format_and_show_instructions(co, r, c_i,  5, 1, NV_COLOR5, 2, u);
            format_and_show_instructions(co, r, c_i,  6, 2, NV_COLOR6, 2, u);
            format_and_show_instructions(co, r, c_i,  7, 0, NV_COLOR7, 3, u);
            format_and_show_instructions(co, r, c_i,  8, 1, NV_COLOR8, 3, u);
            format_and_show_instructions(co, r, c_i,  9, 2, NV_COLOR9, 3, u);
            format_and_show_instructions(co, r, c_i, 10, 0, NV_COLORa, 4, u);
            format_and_show_instructions(co, r, c_i, 11, 1, NV_COLORb, 4, u);
            format_and_show_instructions(co, r, c_i, 12, 2, NV_COLORc, 4, u);
            format_and_show_instructions(co, r, c_i, 13, 0, NV_COLORd, 5, u);
            format_and_show_instructions(co, r, c_i, 14, 1, NV_COLORe, 5, u);
            format_and_show_instructions(co, r, c_i, 15, 2, NV_COLORf, 5, u);
            format_and_show_instructions(co, r, c_i, 16, 0, NV_COLORg, 6, u);
        } 
    lux_set_color(co, c_i[default_color]);
    }
}

void show_hrd_information(unsigned long instr,float r,int d, int u, dyn* b) {

//              Methos for plotting HRD.
      real total_mass=0;
      real t_cur, t_rel, m_rel, m_env, m_core, T_eff, L_eff, p_rot, b_fld;
      real co_core=0;
      t_cur=t_rel=m_rel=m_env=m_core=T_eff=L_eff=p_rot=b_fld=0;
      stellar_type type=NAS;
      story *st = NULL;
      int su[6] = {0, 0, 0, 0, 0, 0};
      real mu[6] = {0, 0, 0, 0, 0, 0};
      for_all_leaves(dyn, b, bi) {
         st = bi->get_starbase()->get_star_story();
         extract_story_chapter(type, t_cur, t_rel, m_rel, m_env, m_core, 
			       co_core,
			       T_eff, L_eff, p_rot, b_fld, *st);
         total_mass += m_env+m_core;
         if (type<= Main_Sequence) {
            su[0]++;
            mu[0] += (m_env + m_core);
         }
         else 
         if (type<= Sub_Giant) {
           su[1]++; 
           mu[1] += (m_env + m_core);
         }
         else 
         if (type<= Super_Giant) {
           su[2]++; 
           mu[2] += (m_env + m_core);
         }
         else 
         if (type== Carbon_Star || type== Helium_Star || type== Helium_Giant ||
	     type== Carbon_Dwarf || type== Helium_Dwarf || type== Oxygen_Dwarf) {
           su[3]++; 
           mu[3] += (m_env + m_core);
         }
         else 
         if (type== Neutron_Star || type== Black_Hole) {
           su[4]++; 
           mu[4] += (m_env + m_core);
         }
         else {
           su[5]++; 
           mu[5] += (m_env + m_core);
         }
      // extract binary information should go here.
      // compute_energies(b, bi, c);
      }

      mu[0] /= total_mass;
      mu[1] /= total_mass;
      mu[2] /= total_mass;
      mu[3] /= total_mass;
      mu[4] /= total_mass;
      mu[5] /= total_mass;

      char *centence = new char[500];
      sprintf(centence,  
"Hertzsprung-Russel Diagram information
  total mass       : %6.1f Mo\n\
  Stellar contents : N       mass fr.\n\
   Main sequence   : %5d   %5.3f\n\
   Early giant     : %5d   %5.3f\n\
   Late type giant : %5d   %5.3f\n\
   Helium rich star: %5d   %5.3f\n\
   inert remnant   : %5d   %5.3f\n\
   other           : %5d   %5.3f\n\n\
   Binary information should go here:\n\
",total_mass, su[0], mu[0], su[1], mu[1], su[2], mu[2], su[3], mu[3],
              su[4], mu[4], su[5], mu[5]); 

    show_instructions(instr, r, centence, u);
    delete centence;
}



