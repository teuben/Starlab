//
// chydro_io.C
//

#include "chydro.h"
#include "util_io.h"

istream & chydro::scan_hydro_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), strcmp(END_HYDRO, input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	char *val = getequals(input_line, keyword);

    	if(0){   // trick to keep the if() statement, for the else below
	}else if(!strcmp("R_eff",keyword)){
	    effective_radius = strtod(val, &val);
    	}else if(!strcmp("M_core",keyword)){
	    core_mass = strtod(val, &val);
	}else if(!strcmp("R_core",keyword)){
	    core_radius = strtod(val, &val);
	}else if(!strcmp("mf",keyword)){
            m_conv_hydro_to_dyn = strtod(val, &val);
	}else if(!strcmp("rf",keyword)){
            r_conv_hydro_to_dyn = strtod(val, &val);
	}else if(!strcmp("tf",keyword)){
            t_conv_hydro_to_dyn = strtod(val, &val);
	}else{
	    add_story_line(hydro_story, input_line);
	}
    }
    return s;
}

ostream& chydro::print_hydro_story(ostream& s)
{
    put_story_header(s, HYDRO_ID);

    put_real_number(s, "  R_eff =  ", effective_radius);
    put_real_number(s, "  M_core =  ", core_mass);
    put_real_number(s, "  R_core =  ", core_radius);

    put_real_number(s, "  mf =  ", m_conv_hydro_to_dyn);
    put_real_number(s, "  rf =  ", r_conv_hydro_to_dyn);
    put_real_number(s, "  tf =  ", t_conv_hydro_to_dyn);

    if (hydro_story)
        put_story_contents(s, *hydro_story);
    put_story_footer(s, HYDRO_ID);
    
    return s;
}

