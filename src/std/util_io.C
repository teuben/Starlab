
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// util_io.C:  functions for formatted I/O.

#include "starlab_vector.h"
#include "util_io.h"
#include "story.h"
#include <ctype.h>

#undef isalnum		/* Hacks for Irix 6.5 <ctype.h> backward compatibility */
#undef isspace

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
    return matchbracket(reference_string, input_line);
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
    if(!matchbracket(reference_string, input_line)){
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

/* Returns 1 if either token matches line, or first two chars of token
 * matches line.  So matchbracket("(Particle", line) matches "(Particle" or "(P".
 */
int matchbracket(const char *token, const char *line) {
  while(*line == ' ' || *line == '\t')
      line++;
  if(token[0] != line[0] || token[1] != line[1])
    return 0;
  return (line[2] == '\0') || (0 == strcmp(token+2, line+2));
}

const char *getequals(const char *input_line, char *keyword)
{
    const char *cp = input_line;

    /* Grab first token from line, like sscanf %s */
    while(isspace(*cp)) cp++;
    int i;
    for(i = 0; isalnum(*cp) || *cp == '_'; )
	keyword[i++] = *cp++;
    keyword[i] = '\0';

    cp = strchr(cp, '=');
    if(cp == NULL) {
	cerr << "Expected keyword = value, but got '"<< input_line <<"'\n";
	exit(1);
    }
    cp++;
    while(isspace(*cp)) cp++;
    return cp;
}

void set_vector_from_input_line(vec & v, char * input_line)
{
    char *eq = strchr(input_line, '=');
    if(eq)
	set_vector_from_string( v, eq+1 );
}

void set_vector_from_string(vec & v, char *val)
{
    real component[3];
    char *cp, *ep;
    component[0] = strtod(val, &cp);
    component[1] = strtod(cp, &cp);
    component[2] = strtod(cp, &ep);
    if(cp == ep) {
	cerr << "Expected three reals, got: " << val << endl;
	exit(1);
    }
    v = vec(component[0],component[1],component[2]);
}
    

static bool print = true;

xreal get_xreal_from_input_line(char * input_line)
{
    char *val = strchr(input_line, '=');
    if(val == NULL) return (xreal)0;
    val++;

#if defined USE_XREAL

    // "True" xreal:

    char *sp, *ep;
    long long i = STRTOL(val, &sp, 10);		  // signed integer part
    unsigned long long f = STRTOUL(sp, &ep, 0);   // unsigned fractional part
						  // "0" here means that we
						  // can read hex or integer

    if (sp == ep) {				  // if we didn't get both
						  // of above,

	// Hmmm... most likely we have real input data.  Try just reading
	// a real number.  (Steve, 6/00)

	if (print) {
	    cerr << "get_xreal_from_input_line: error reading xreal input "
		 << "from line" << endl
		 << "    " << input_line << endl
		 << "Assuming real data." << endl << endl;
	    print = false;
	}

	return (xreal)strtod( val, NULL );
    }

    return xreal(i, f);

#else

    // xreal is really just real:

    return (xreal)strtod( val, NULL );

#endif
}

//
//---------------------------------------------------------------------
//
// Horrible kludges:
// ----------------
//

static bool short_story_keywords = false;

bool use_short_story_keywords( bool useshort ) {
  bool was = short_story_keywords;
  short_story_keywords = useshort;
  return was;
}


void put_story_header(ostream & s, char * id)
{
#ifdef BAD_GNU_IO
    if (s == cout) {
	if(short_story_keywords) {
	    fprintf(stdout, "(%c\n", id[0]);
	} else {
	    fprintf(stdout, "(%s\n", id);
	}
    }
    else
#endif
    if(short_story_keywords) {
	s << '(' << id[0] << endl;
    } else {
	s << '(' << id << endl;
    }
}

void put_story_footer(ostream & s, char * id)
{
#ifdef BAD_GNU_IO
    if (s == cout) {
	if(short_story_keywords) {
	    fprintf(stdout, ")%c\n", id[0]);
	} else {
	    fprintf(stdout, ")%s\n", id);
	}
    }
    else
#endif
    if(short_story_keywords) {
	s << ')' << id[0] << endl;
    } else {
	s << ')' << id << endl;
    }
}


#ifdef USE_XREAL
void put_real_number(ostream & s, char * label, xreal x)
{
    // Simplest version:

    // s << label << x.get_i() << " " << x.get_f() << endl;

    // Better: Use hex for the fractional part...

    xfrac_t f = x.get_f();
    char tmp[128];

    if (f == 0)
	sprintf(tmp, "0");			// handy
    else
	sprintf(tmp, "%#16.16llx", f);

    s << label << x.get_i() << " " << tmp << endl;
}
#endif

// NOTE use of precision here.  If the STARLAB_PRECISION environment
// variable is set, the first call to set_starlab_precision will use
// its value (whatever it may be).  Subsequent calls will return the
// same value originally read from the environment.  If no environment
// variable is set, a default value (currently 18) is used.

void put_real_number(ostream & s, char * label, real x)
{
    int old_precision = set_starlab_precision(s);

#ifndef BAD_GNU_IO
    s << label << x << endl;
#else
    if (s == cout) {
	static char format[40];
	static int local_precision = -1;
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

void put_real_vector(ostream & s, char * label, vec v)
{
    int old_precision = set_starlab_precision(s);

#ifndef BAD_GNU_IO
    s << label << v << endl;
#else
    if (s == cout){
	static char format[40];
	static int local_precision = -1;
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
