//
// hydro_io.C
//

#include "hydro.h"
#include "util_io.h"

istream & hydro::scan_hydro_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), strcmp(END_HYDRO, input_line)){

	char keyword[MAX_INPUT_LINE_LENGTH];
	char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
	sscanf(input_line,"%s%s",keyword,should_be_equal_sign);
	if(strcmp("=",should_be_equal_sign)){
	    cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
	    exit(1);
	}

    	if(0){   // trick to keep the if() statement, for the else below
	}else if(!strcmp("R_eff",keyword)){
	    sscanf(input_line,"%*s%*s%lf",&effective_radius);
	}else if(!strcmp("mf",keyword)){
            sscanf(input_line,"%*s%*s%lf",&m_conv_hydro_to_dyn);
	}else if(!strcmp("rf",keyword)){
            sscanf(input_line,"%*s%*s%lf",&r_conv_hydro_to_dyn);
	}else if(!strcmp("tf",keyword)){
            sscanf(input_line,"%*s%*s%lf",&t_conv_hydro_to_dyn);
	}else{
	    add_story_line(hydro_story, input_line);
	}
    }
    return s;
}

ostream& hydro::print_hydro_story(ostream& s)
{
    put_story_header(s, HYDRO_ID);

    put_real_number(s, "  R_eff =  ", effective_radius);

    put_real_number(s, "  mf =  ", m_conv_hydro_to_dyn);
    put_real_number(s, "  rf =  ", r_conv_hydro_to_dyn);
    put_real_number(s, "  tf =  ", t_conv_hydro_to_dyn);

    if (hydro_story)
        put_story_contents(s, *hydro_story);
    put_story_footer(s, HYDRO_ID);
    
    return s;
}

