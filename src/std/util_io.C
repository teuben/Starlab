
//
// util_io.C
//

#include "starlab_vector.h"
#include "util_io.h"
#include "story.h"

int get_line(istream & s, char * line)
{
    s.get(line,MAX_INPUT_LINE_LENGTH,'\n');
    char c;
    if(s.get(c) && c!='\n'){
	cerr << "get_line : input line too long :'"<<line<<"'\n";
	exit(1);
    }
    return strlen(line);
}

int check_input_line(istream &s, char* reference_string)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    get_line(s,input_line);
    if(s.eof()){
	cerr << "check_input_line : unexpected EOF, expected '";
	cerr << reference_string <<"'\n";
	exit(1);
    }
    return !strcmp(input_line,reference_string);
}    

int check_and_skip_input_line(istream &s, char* reference_string)
{
//  char dummy;
    char input_line[MAX_INPUT_LINE_LENGTH];
    while (! get_line(s,input_line) && !s.eof())
	{
//	if(s.eof())
//	    return 0;
//	cin >> dummy;           // to get the g++ compiler to the EOF
// cerr << "  dummy = " << dummy << endl;
//	if(s.eof())
//	    return 0;
//	if (dummy == '\0')

	int shhh = 1;           // to make the SG compiler shut up
	if (shhh)               // to make the SG compiler shut up
	    return 0;
	}
    if(strcmp(input_line,reference_string) != 0 ){
	if(s.eof()){
	    return 0;
	}else{
	    cerr << "Input line must be '"<<reference_string;
	    cerr <<"', I got '" << input_line << "'\n";
	    exit(1);
	}
    }
    return 1;
}    

int get_data_line(istream & s,char * input_line)
{
    get_line(s,input_line);
    return strcmp(input_line,")");
}


void set_vector_from_input_line(vector & v, char * input_line)
{
    real vector_component[3];
    sscanf(input_line,"%*s%*s%lf%lf%lf",vector_component,
	   vector_component+1,vector_component+2);
    v=vector(vector_component[0],vector_component[1],vector_component[2]);
}

static bool print = true;

xreal get_xreal_from_input_line(char * input_line)
{
#if defined USE_XREAL

    // "True" xreal:

    long long i;
    unsigned long long f;
    int n = sscanf(input_line, "%*s%*s%Ld%Lu",&i, &f);

    if (n < 2) {

	// Hmmm... most likely we have real input data.  Try just reading
	// a real number.  (Steve, 6/00)

	if (print) {
	    cerr << "get_xreal_from_input_line: error reading xreal input "
		 << "from line" << endl
		 << "    " << input_line << endl
		 << "Assuming real data." << endl << endl;
	    print = false;
	}

	real x;
	sscanf(input_line, "%*s%*s%lf", &x);
	return (xreal)x;
    }

    return xreal(i, f);

#else

    // xreal is really just real:

    xreal x;
    sscanf(input_line, "%*s%*s%lf", &x);
    return x;

#endif
}

//
//---------------------------------------------------------------------
//
// Horrible kludges:
// ----------------
//

void put_story_header(ostream & s, char * id)
{
#ifndef BAD_GNU_IO
    s << '(' << id << endl;
#else
    if (s == cout)
	fprintf(stdout, "(%s\n", id);
    else
	s << '(' << id << endl;
#endif
}

void put_story_footer(ostream & s, char * id)
{
#ifndef BAD_GNU_IO
    s << ')' << id << endl;
#else
    if (s == cout)
	fprintf(stdout, ")%s\n", id);
    else
	s << ')' << id << endl;
#endif
}

#ifdef USE_XREAL
void put_real_number(ostream & s, char * label, xreal x)
{
    s << label << x.get_i() << " " << x.get_f() << endl;
}
#endif

// NOTE use of precision here.  If the STARLAB_PRECISION environment
// variable is set, the first call to set_starlab_precision will use
// its value (whatever it may be).  Subsequent calls will return the
// same value originally read from the environment.  If no environment
// variable is set, a default value (currently 18) is used.

void put_real_number(ostream & s, char * label, real x)
{
    static char format[40];
    static int  local_precision = -1;

    int old_precision = set_starlab_precision(s);

#ifndef BAD_GNU_IO
    s << label << x << endl;
#else
    if (s == cout) {
	int precision = get_starlab_precision();
	if (local_precision != precision) {
	    local_precision = precision;
	    int p = precision;
	    if (p < 0) p = 5;
	    sprintf(format, "%%s%%.%dg\n", p);
	}
	fprintf(stdout, format, label, x);
    } else
	s << label << x << endl;
#endif

    // Restore the current precision.

    if (get_starlab_precision() != old_precision)
	s.precision(old_precision);
}

void put_real_vector(ostream & s, char * label, vector v)
{
    static char format[40];
    static int  local_precision = -1;

    int old_precision = set_starlab_precision(s);

#ifndef BAD_GNU_IO
    s << label << v << endl;
#else
    if (s == cout){
	int precision = get_starlab_precision();
	if (local_precision != precision) {
	    local_precision = precision;
	    int p = precision;
	    if (p < 0) p = 5;
	    sprintf(format,"%%s%%.%dg %%.%dg %%.%dg\n", p, p, p);
	}
	fprintf(stdout, format,	label, v[0], v[1], v[2]);
    } else
	s << label << v << endl;
#endif

    // Restore the current precision.

    if (get_starlab_precision() != old_precision)
	s.precision(old_precision);
}

void put_integer(ostream & s, char * label, int i)
{
#ifndef BAD_GNU_IO
    s << label << i << endl;
#else
    if (s == cout)
	fprintf(stdout, "%s%d\n", label, i);
    else
	s << label << i << endl;
#endif
}

void put_string(ostream & s, char * label, char * str)
{
#ifndef BAD_GNU_IO
    s << label << str << endl;
#else
    if (s == cout) 
	fprintf(stdout, "%s%s\n", label, str);
    else
	s << label << str << endl;
#endif
}
