//
// hydrobase_io.C
//

#include "hydrobase.h"
#include "util_io.h"

istream & hydrobase::scan_hydro_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while(get_line(s,input_line), strcmp(END_HYDRO, input_line)) {
	char keyword[MAX_INPUT_LINE_LENGTH];
	char should_be_equal_sign[MAX_INPUT_LINE_LENGTH];
	sscanf(input_line,"%s%s",keyword,should_be_equal_sign);
	if(strcmp("=",should_be_equal_sign)){
	    cerr << "Expected '=', but got '"<< should_be_equal_sign <<"'\n";
	    exit(1);
	}

    	add_story_line(hydro_story, input_line);
    }
    return s;
}

ostream& hydrobase::print_hydro_story(ostream& s)
{
    put_story_header(s, HYDRO_ID);

    if (hydro_story)
        put_story_contents(s, *hydro_story);

    put_story_footer(s, HYDRO_ID);
    
    return s;
}
