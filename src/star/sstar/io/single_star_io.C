//
// single_star_io.C
//

//#include "hdyn.h"
#include "single_star.h"
#include "util_io.h"

istream & single_star::scan_star_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    cerr << "\n\nUsing single_star_io version of scan_star_story()...\n\n";

    while(get_line(s,input_line), strcmp(END_STAR, input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
	sscanf(input_line,"%s%s",keyword,should_be_equal_sign);
	if(strcmp("=",should_be_equal_sign)){
	    cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
	    exit(1);
	}

        real number;
    	if(0){   // trick to keep the if() statement, for the else below
        }else if(!strcmp("Type",keyword)){
            char str_tpe[MAX_INPUT_LINE_LENGTH];
            sscanf(input_line,"%*s%*s%s",str_tpe);
        }else if(!strcmp("Class",keyword)){
            char str_cls[MAX_INPUT_LINE_LENGTH];
            sscanf(input_line,"%*s%*s%s",str_cls);
        }else if(!strcmp("T_cur",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_current_time(number);
        }else if(!strcmp("T_rel",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_relative_age(number);
        }else if(!strcmp("M_rel",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_relative_mass(number);
        }else if(!strcmp("M_env",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_envelope_mass(number);
        }else if(!strcmp("M_core",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_core_mass(number);
        }else if(!strcmp("M_COcore",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_COcore_mass(number);
        }else if(!strcmp("T_eff",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
        }else if(!strcmp("L_eff",keyword)){
            sscanf(input_line,"%*s%*s%lf",&number);
            set_luminosity(number);
        }else if(!strcmp("P_rot",keyword)){       // Experimental extra 
	  sscanf(input_line,"%*s%*s%lf",&number); // information for 
	  set_rotation_period(number);            // radio pulsars.
        }else if(!strcmp("B_fld",keyword)){
	  sscanf(input_line,"%*s%*s%lf",&number);
	  set_magnetic_field(number);
	}else{
	    add_story_line(star_story, input_line);
	}
    }

    //if(current_time<=0) current_time = relative_age;
    return s;
}

ostream& single_star::print_star_story(ostream& s,
				       int short_output)  // default = 0
{
    put_story_header(s, STAR_ID);

    put_string(s,      "  Type   =  ", type_string(get_element_type()));
    if (get_spec_type(Accreting)==Accreting)
      put_string(s,      "  Class   =  ", type_string(get_spec_type(Accreting)));
    put_real_number(s, "  T_cur  =  ", get_current_time());
    if (get_current_time()!=get_relative_age())
      put_real_number(s, "  T_rel  =  ", get_relative_age());
    put_real_number(s, "  M_rel  =  ", get_relative_mass());
    put_real_number(s, "  M_env  =  ", get_envelope_mass());
    put_real_number(s, "  M_core =  ", get_core_mass());
    put_real_number(s, "  T_eff  =  ", temperature());
    put_real_number(s, "  L_eff  =  ", get_luminosity());

    // Extra output for stars with CO cores
    if (star_with_COcore()) {
      put_real_number(s, "  M_COcore  =  ", get_COcore_mass());
    }

    // Extra output for radio- and X-ray pulsars.
    if (get_element_type()==Xray_Pulsar ||
	get_element_type()==Radio_Pulsar ||
	get_element_type()==Neutron_Star) {
      put_real_number(s, "  P_rot  =  ", get_rotation_period());
      put_real_number(s, "  B_fld  =  ", get_magnetic_field());
    }

    if (star_story)
      put_story_contents(s, *star_story);

    put_story_footer(s, STAR_ID);

    return s;
}

void extract_line_text(stellar_type& type, real& t_cur, real& t_rel, 
		       real& m_rel, real& m_env, real& m_core, real& co_core,
		       real& t_eff, real& l_eff,
		       real& p_rot, real& b_fld,
		       story& s)
    {

        char keyword[MAX_INPUT_LINE_LENGTH];
        char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
        char line [MAX_INPUT_LINE_LENGTH];
	char type_string[MAX_INPUT_LINE_LENGTH];
        strcpy(line, s.get_text());
        sscanf(line,"%s%s",keyword,should_be_equal_sign);
        if(strcmp("=",should_be_equal_sign)){
            cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
            exit(1);
        }

        real number;
        if(0){   // trick to keep the if() statement, for the else below
        }else if(!strcmp("Type",keyword)){
            char str_tpe[MAX_INPUT_LINE_LENGTH];
            sscanf(line,"%*s%*s%s",type_string); 
	    type = extract_stellar_type_string(type_string);
        }else if(!strcmp("T_cur",keyword)){
            sscanf(line,"%*s%*s%lf",&t_cur);
        }else if(!strcmp("T_rel",keyword)){
            sscanf(line,"%*s%*s%lf",&t_rel);
        }else if(!strcmp("M_rel",keyword)){
            sscanf(line,"%*s%*s%lf",&m_rel);
        }else if(!strcmp("M_env",keyword)){
            sscanf(line,"%*s%*s%lf",&m_env);
        }else if(!strcmp("M_core",keyword)){
            sscanf(line,"%*s%*s%lf",&m_core);
        }else if(!strcmp("M_COcore",keyword)){
            sscanf(line,"%*s%*s%lf",&co_core);
        }else if(!strcmp("T_eff",keyword)){
            sscanf(line,"%*s%*s%lf",&t_eff);
        }else if(!strcmp("L_eff",keyword)){
            sscanf(line,"%*s%*s%lf",&l_eff);
        }else if(!strcmp("P_rot",keyword)){       // Experimental extra 
	    sscanf(line,"%*s%*s%lf",&p_rot);     // information for 
        }else if(!strcmp("B_fld",keyword)){       // Only for neutron stars
  	    sscanf(line,"%*s%*s%lf",&b_fld);
        }

    }

void extract_story_chapter(stellar_type& type, real& t_cur, real& t_rel, 
                           real& m_rel, real& m_env, real& m_core, 
			   real& co_core,
                           real& T_eff, real& L_eff,
			   real& p_rot, real& b_fld, story& s)
{
    
    if (!s.get_chapter_flag())
        {
        cerr << "extract_story_chapter: not a story\n";
        exit(1);
        }

    t_rel=-1;
      
    for (story * d = s.get_first_daughter_node();
	 d != NULL;
	 d = d->get_next_story_node()) {

      if (d->get_chapter_flag()) 
	extract_story_chapter(type, t_cur, t_rel,
			      m_rel, m_env, m_core, co_core,
			      T_eff, L_eff, p_rot, b_fld, *d);
      else 
	extract_line_text(type, t_cur, t_rel,
			  m_rel, m_env, m_core, co_core,
			  T_eff, L_eff, p_rot, b_fld, *d);
    }
    
    if (t_rel<=0)
      t_rel = t_cur;
}

void single_star::star_transformation_story(stellar_type new_type)
{
    char info_line[MAX_STORY_LINE_LENGTH];
    stellar_type old_type = get_element_type();
    real time = get_current_time();

    sprintf(info_line,
	    "%s_to_%s_at_time = %6.2f",
	    type_string(old_type), type_string(new_type), time);

    add_story_line(get_node()->get_log_story(), info_line);

    // cerr output for debugging!

    cerr << endl << get_node()->format_label() << " " << info_line
	 << " Myr (old mass = " << get_total_mass() << ")" << endl;

    if(is_binary_component()) {
	cerr << "binary companion: " 
	     << type_string(get_companion()->get_element_type()) << endl;
	cerr << "parent: " << get_node()->get_parent()->format_label() 
	     << endl;
    }
}

void single_star::post_supernova_story()
{
    char info_line[MAX_STORY_LINE_LENGTH];

    sprintf(info_line, "v_kick = %6.2f, %6.2f, %6.2f",
	    anomal_velocity[0], anomal_velocity[1], anomal_velocity[2]);

    add_story_line(get_node()->get_log_story(), info_line);
}

void single_star::first_roche_lobe_contact_story(stellar_type accretor)
{
    char info_line[MAX_STORY_LINE_LENGTH];
    real time = get_current_time();
     
    sprintf(info_line,
	    "%s_%s_%s_RLOF_to_%s_at_time = %6.2f",
	    type_string(get_element_type()),
	    type_string(get_spectral_class(temperature())),
	    type_string(get_luminosity_class(temperature(), luminosity)),
	    type_string(accretor), time);

    add_story_line(get_node()->get_log_story(), info_line);

}
